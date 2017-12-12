from math import erf

import numpy as np
from brian2 import units, Equations, NeuronGroup, check_units
from scipy.integrate import quad

from meanfield.populations.MFPop import MFPop
from meanfield.utils import lazyproperty
from meanfield.parameters import NP, SP


class MFLinearPop(MFPop):
    """pop: similar neurons"""

    def __init__(self, name, n, params):
        super().__init__(name, n, params)

        defaults = {
        }
        expectations = {
            NP.GM: units.siemens,
            NP.VL: units.volt,
            NP.CM: units.farad,
            NP.VTHR: units.volt,
            NP.VRES: units.volt,
            NP.TAU_RP: units.second
        }

        self.params.fill(defaults)
        self.params.verify(expectations)

    def brian2_model(self):
        eqs = Equations(
            'dv / dt = (- g * (v - vl) - I) / cm : volt (unless refractory)',
            g=self.params[NP.GM],
            vl=self.params[NP.VL],
            cm=self.params[NP.CM]
        )

        total = []
        for s in self.sources + self.noises:
            eqs += s.brian2_model()
            total.append(s.current_name)

        if len(total):
            eqs += 'I = {} : amp'.format(' + '.join(total))
        else:
            eqs += 'I = 0 : amp'

        return eqs

    def brian2_threshold(self):
        return 'v > {} * mV'.format(self.params[NP.VTHR] / units.mV)

    def brian2_reset(self):
        return 'v = {} * mV'.format(self.params[NP.VRES] / units.mV)

    @lazyproperty
    def brian2(self):
        P = NeuronGroup(
            self.n,
            self.brian2_model(),
            method='euler',
            threshold=self.brian2_threshold(),
            reset=self.brian2_reset(),
            refractory=self.params[NP.TAU_RP]
        )
        P.v = self.params[NP.VRES]
        return P


    @property
    @check_units(result=units.siemens)
    def total_cond(self):
        """
        Gm * SE in [1]
        Units of S
        """
        return self.params[NP.GM] + np.sum(s.conductance for s in self.sources + self.noises)

    @property
    @check_units(result=units.second)
    def tau_eff(self):
        """
        Seconds
        """
        return self.params[NP.CM] / self.total_cond

    @property
    @check_units(result=units.volt)
    def mu(self):
        """
        Volt
        """
        return np.sum(s.voltage_conductance for s in self.sources + self.noises) / self.total_cond

    @property
    @check_units(result=units.volt ** 2)
    def sigma_square(self):
        """
        Volt^2
        """
        if not len(self.noises):
            return 0. * units.volt ** 2

        noise = self.noises[0]
        return (noise.g_base / self.params[NP.CM] * (self.v_mean - noise.params[SP.VREV])) ** 2 * self.tau_eff * noise.g_dyn() * noise.params[SP.TAU]

    @check_units(result=units.Hz)
    def phi_firing_func(self):
        # FIXME no noise
        if not len(self.noises):
            return 0 * units.Hz

        sigma = np.sqrt(self.sigma_square)
        tau_eff = self.tau_eff
        noise = self.noises[0]

        # Brunel Wang 2001 / Brunel Sergi 1998
        beta = (self.params[NP.VRES] - self.params[NP.VL] - self.mu) / sigma
        alpha = -0.5 * noise.params[SP.TAU] / tau_eff \
                + 1.03 * np.sqrt(noise.params[SP.TAU] / tau_eff) \
                + (- self.mu - self.params[NP.VL] + self.params[NP.VTHR]) * (
            1. + (0.5 * noise.params[SP.TAU] / tau_eff)) / sigma

        # Fourcauld Brunel 2002
        #beta = (self.params[NP.VRES] - self.params[NP.VL] - self.mu) / sigma + 1.03 * np.sqrt(noise.params[SP.TAU] / tau_eff)
        #alpha = 1.03 * np.sqrt(noise.params[SP.TAU] / tau_eff) \
        #        + (- self.mu - self.params[NP.VL] + self.params[NP.VTHR])/ sigma

        def integrand(x, max_exp=50):
            if x < -max_exp:
                return np.exp(-max_exp ** 2) * (1. + erf(-max_exp))
            if x > max_exp:
                return 0.
            return np.exp(x ** 2) * (1. + erf(x))


        #import time
        #s = time.time()
        q = quad(integrand, beta, alpha, limit=10000)
        #print(q)
        #print('->', time.time() - s)

        return 1. / (self.params[NP.TAU_RP] + tau_eff * np.sqrt(np.pi) * q[0])

    @property
    @check_units(result=units.Hz)
    def rate_prediction(self):
        return self.phi_firing_func()

    @property
    @check_units(result=units.volt)
    def v_mean_prediction(self):
        """
        Volt
        """
        return self.params[NP.VL] + self.mu - (self.params[NP.VTHR] - self.params[NP.VRES]) * self.rate * self.tau_eff

