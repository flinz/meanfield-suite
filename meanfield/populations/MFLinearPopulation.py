from math import erf
from types import MappingProxyType

import numpy as np
from brian2 import units, Equations, NeuronGroup, check_units
from scipy.integrate import quad

from meanfield.parameters import PP, IP
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

    def __init__(self, name, n, params):
        super().__init__(name, n, params)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    @property
    @check_units(result=units.siemens)
    def total_cond(self):
        """
        Gm * SE in [1]
        Units of S
        """
        return self[PP.GM] + np.sum(s.conductance for s in self.inputs + self.noises)

    @property
    @check_units(result=units.second)
    def tau_eff(self):
        """
        Seconds
        """
        return self[PP.CM] / self.total_cond

    @property
    @check_units(result=units.volt)
    def mu(self):
        """
        Volt
        """
        return np.sum(s.voltage_conductance for s in self.inputs + self.noises) / self.total_cond

    @property
    @check_units(result=units.volt ** 2)
    def sigma_square(self):
        """
        Volt^2
        """
        if not len(self.noises):
            return 0. * units.volt ** 2

        # only single source of noise supported
        noise = self.noises[0]
        return (noise.g_base / self[PP.CM] * (self.v_mean - noise.params[IP.VREV])) ** 2 * self.tau_eff * noise.g_dyn() * noise.params[IP.TAU]

    @property
    @check_units(result=units.Hz)
    def rate_prediction(self):
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
        alpha = -0.5 * noise.params[IP.TAU] / tau_eff \
                + 1.03 * np.sqrt(noise.params[IP.TAU] / tau_eff) \
                + (- self.mu - self[PP.VL] + self[PP.VTHR]) * (
                        1. + (0.5 * noise.params[IP.TAU] / tau_eff)) / sigma

        # Fourcauld Brunel 2002
        # beta = (self.params[NP.VRES] - self.params[NP.VL] - self.mu) / sigma + 1.03 * np.sqrt(noise.params[SP.TAU] / tau_eff)
        # alpha = 1.03 * np.sqrt(noise.params[SP.TAU] / tau_eff) \
        #        + (- self.mu - self.params[NP.VL] + self.params[NP.VTHR])/ sigma

        def integrand(x, max_exp=20):
            if x < -max_exp:
                return np.exp(-max_exp ** 2) * (1. + erf(-max_exp))
            if x > max_exp:
                return 0.
            return np.exp(x ** 2) * (1. + erf(x))

        # import time
        # s = time.time()
        q = quad(integrand, beta, alpha, limit=10000)
        # print(q)
        # print('->', time.time() - s)

        return 1. / (self[PP.TAU_RP] + tau_eff * np.sqrt(np.pi) * q[0])

    @property
    @check_units(result=units.volt)
    def v_mean_prediction(self):
        """
        Volt
        """
        return self[PP.VL] + self.mu - (self[PP.VTHR] - self[PP.VRES]) * self.rate * self.tau_eff

    # Simulation

    @lazyproperty
    def brian2(self):
        pop = NeuronGroup(
            self.n,
            self.brian2_model(),
            method='euler',
            threshold='v > {} * mV'.format(self[PP.VTHR] / units.mV),
            reset='v = {} * mV'.format(self[PP.VRES] / units.mV),
            refractory=self[PP.TAU_RP],
            name=self.ref
        )
        pop.v = self[PP.VRES]
        return pop

    def brian2_model(self):
        eqs = Equations(
            'dv / dt = (- g * (v - vl) - I) / cm : volt (unless refractory)',
            g=self[PP.GM],
            vl=self[PP.VL],
            cm=self[PP.CM]
        )

        all_currents = []
        for s in self.inputs + self.noises:
            eqs += s.brian2_model()
            all_currents.append(s.current_name)

        if len(all_currents):
            eqs += 'I = {} : amp'.format(' + '.join(all_currents))
        else:
            eqs += 'I = 0 : amp'

        return eqs
