from math import erf
from types import MappingProxyType
from typing import Union, Optional

import numpy as np
from brian2 import units, Equations, NeuronGroup, check_units, BrianObject
from scipy.integrate import quad

from meanfield.parameters import PP, IP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import lazyproperty


class MFLinearPopulation(MFPopulation):
    """pop: similar neurons"""

    arguments = MappingProxyType({
        PP.GM: units.siemens,
        PP.VL: units.volt,
        PP.CM: units.farad,
        PP.VTHR: units.volt,
        PP.VRES: units.volt,
        PP.TAU_RP: units.second
    })

    defaults = MappingProxyType({})

    def __init__(self, n: int, parameters: Union[dict, MFParams], **kwargs):
        super().__init__(n, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    @property
    @check_units(result=units.siemens)
    def total_cond(self) -> units.siemens:
        """
        Gm * SE in [1]
        Units of S
        """
        return self[PP.GM] + np.sum(s.conductance for s in self.inputs + self.noises)

    @property
    @check_units(result=units.second)
    def tau_eff(self) -> units.second:
        """
        Seconds
        """
        return self[PP.CM] / self.total_cond

    @property
    @check_units(result=units.volt)
    def mu(self) -> units.volt:
        """
        Volt
        """
        return np.sum(s.voltage_conductance for s in self.inputs + self.noises) / self.total_cond

    @property
    @check_units(result=units.volt ** 2)
    def sigma_square(self) -> units.volt ** 2:
        """
        Volt^2
        """
        if not len(self.noises):
            return 0. * units.volt ** 2

        # only single source of noise supported
        noise = self.noises[0]
        return (noise[IP.GM] / self[PP.CM] * (self.v_mean - noise[IP.VREV])) ** 2 * self.tau_eff * noise.g_dyn() * noise[IP.TAU]

    @property
    @check_units(result=units.Hz)
    def rate_prediction(self) -> units.Hz:
        '''
        phi_firing_func
        '''
        # FIXME no noise
        if not len(self.noises):
            return 0 * units.Hz

        sigma = np.sqrt(self.sigma_square)
        tau_eff = self.tau_eff
        noise = self.noises[0]

        # Brunel Wang 2001 / Brunel Sergi 1998
        beta = (self[PP.VRES] - self[PP.VL] - self.mu) / sigma
        alpha = -0.5 * noise[IP.TAU] / tau_eff \
                + 1.03 * np.sqrt(noise[IP.TAU] / tau_eff) \
                + (- self.mu - self[PP.VL] + self[PP.VTHR]) * (
                        1. + (0.5 * noise[IP.TAU] / tau_eff)) / sigma

        # Fourcauld Brunel 2002
        # beta = (self[NP.VRES] - self[NP.VL] - self.mu) / sigma + 1.03 * np.sqrt(noise[SP.TAU] / tau_eff)
        # alpha = 1.03 * np.sqrt(noise[SP.TAU] / tau_eff) \
        #        + (- self.mu - self[NP.VL] + self[NP.VTHR])/ sigma

        # FIXME very called
        def integrand2(x, max_exp=25):
            if x < -max_exp:
                return np.exp(-max_exp ** 2) * (1. + erf(-max_exp))
            if x > max_exp:
                return 0.
            return np.exp(x ** 2) * (1. + erf(x))

        import scipy.special

        def integrand(x):
            return np.exp(x ** 2) * (1 + scipy.special.erf(x))

        # import time
        # s = time.time()
        q = quad(integrand, beta, alpha, limit=100) # FIXME quadpad too costly
        # print(q)
        # print('->', time.time() - s)

        return 1. / (self[PP.TAU_RP] + tau_eff * np.sqrt(np.pi) * q[0])

    @property
    @check_units(result=units.volt)
    def v_mean_prediction(self) -> units.volt:
        """
        Volt
        """
        return self[PP.VL] + self.mu - (self[PP.VTHR] - self[PP.VRES]) * self.rate * self.tau_eff

    # Simulation

    @lazyproperty
    def brian2(self) -> BrianObject:
        pop = NeuronGroup(
            self.n,
            self.brian2_model,
            method='euler',
            threshold=f'v > {self[PP.VTHR] / units.mV} * mV',
            reset=f'v = {self[PP.VRES] / units.mV} * mV',
            refractory=self[PP.TAU_RP],
            name=self.ref
        )
        pop.v = self[PP.VRES]
        return pop

    @property
    def brian2_model(self) -> Optional[Equations]:
        eqs = Equations(
            'dv / dt = (- g * (v - vl) - I) / cm : volt (unless refractory)',
            g=self[PP.GM],
            vl=self[PP.VL],
            cm=self[PP.CM]
        )

        all_currents = []
        for s in self.inputs + self.noises:
            eqs += s.brian2_model
            all_currents.append(s.current_name)

        if len(all_currents):
            sum_currents = ' + '.join(all_currents)
            eqs += f'I = {sum_currents} : amp'
        else:
            eqs += 'I = 0 : amp'

        return eqs
