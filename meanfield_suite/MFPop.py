from abc import abstractproperty
from math import erf

import numpy as np
from brian2 import units, Equations, NeuronGroup
from scipy.integrate import quad

from MFParams import MFParams
from params import NP


class MFPop(object):

    def __init__(self, name, n, params=None):
        self.name = name
        self.n = n
        self.params = MFParams({}) if params is None else params

        self.sources = []
        self.noise = None
        # base estimation
        self.rate = 0.
        self.v_mean = -60.

    def add_noise(self, noise):
        self.noise = noise

    def add_source(self, source):
        self.sources.append(source)

    @property
    def rate_hz(self):
        return self.rate

    @rate_hz.setter
    def rate_hz(self, value):
        self.rate = value

    @property
    def rate_ms(self):
        return self.rate / 1e3

    @rate_ms.setter
    def rate_ms(self, value):
        self.rate = value * 1e3

    @abstractproperty
    def rate_prediction(self):
        pass

    @abstractproperty
    def v_mean_prediction(self):
        pass


class MFLinearPop(MFPop):
    """pop: similar neurons"""

    def __init__(self, name, n, params):
        super().__init__(name, n, params)
        expectation = {
            NP.GM: units.siemens,
            NP.VL: units.volt,
            NP.CM: units.farad,
            NP.VTHR: units.volt,
            NP.VRES: units.volt,
            NP.TAU_RP: units.second
        }
        params.verify(expectation)
        self.params = params

    def brian_v(self):
        eqs = Equations(
            'dv / dt = (- g * (v - vl) - I) / cm : volt (unless refractory)',
            g=self.params[NP.GM],
            vl=self.params[NP.VL],
            cm=self.params[NP.CM]
        )

        total = []
        for i, s in enumerate(self.sources):
            id = 'i' + str(i)
            eqs += s.brian2_model(str(i), id)
            total.append(id)

        eqs += 'I = ' + '+'.join(total) + ': amp'
        return eqs

    @property
    def brian2(self):
        thres = 'v > {}'.format(self.params[NP.VTHR])
        reset = 'v = {}'.format(self.params[NP.VRES])
        P = NeuronGroup(self.n, self.brian_v(), method='euler', threshold=thres, reset=reset, refractory=self.params[NP.TAU_RP])
        P.v = self.params[NP.VRES]
        return P

    # TODO : when parameter, source, modelling, when to construct

    @property
    def total_cond(self):
        """
        Gm * SE in [1]
        Units of S
        """
        return self.params[NP.GM] / units.siemens +\
               np.sum(s.conductance for s in self.sources)

    @property
    def tau_eff(self):
        """
        Seconds
        """
        return self.params[NP.CM] / units.farad / self.total_cond

    @property
    def mu(self):
        """
        Volt
        """
        return 1. / self.total_cond * np.sum(s.voltage_conductance for s in self.sources)

    @property
    def sigma_square(self):
        """
        Volt^2
        """
        if not self.noise:
            return 0.

        return (self.noise.g_base / self.params[NP.CM] * units.farad * (
            self.v_mean - self.noise.E_rev)) ** 2 * self.tau_eff * self.noise.g_dyn() * self.noise.noise_tau

    def phi_firing_func(self):

        sigma = np.sqrt(self.sigma_square)
        tau_eff = self.tau_eff

        beta = (self.params[NP.VRES] / units.volt - self.params[NP.VL] / units.volt - self.mu) / sigma
        alpha = -0.5 * self.noise.noise_tau / tau_eff \
                + 1.03 * np.sqrt(self.noise.noise_tau / tau_eff) \
                + (-self.mu - self.params[NP.VL] / units.volt + self.params[NP.VTHR] / units.volt) * (
            1. + (0.5 * self.noise.noise_tau / tau_eff)) / sigma

        def integrand(x):
            if x < -10.:
                return np.exp(10. ** 2) * (1. + erf(10.))
            if x > 10.:
                return 0.
            return np.exp(x ** 2) * (1. + erf(x))

        return 1. / (self.params[NP.TAU_RP] / units.second + tau_eff * np.sqrt(np.pi) * quad(integrand, beta, alpha)[0])

    @property
    def rate_prediction(self):
        return 1e3 * self.phi_firing_func()

    @property
    def v_mean_prediction(self):
        """
        Volt
        """
        return self.params[NP.VL] / units.volt + self.mu - (self.params[NP.VTHR] - self.params[NP.VRES]) / units.volt * self.rate_ms * self.tau_eff

    def __repr__(self):
        return "MFpop [%s] <%s (%i sources, n: %i, rate: %.4f, v_mean: %.4f)>" % \
               (id(self), self.name, len(self.sources), self.n, self.rate_hz, self.v_mean)

    def print_sys(self, mf=False):
        print("\t%s - tau_eff: %.1fms, mu: %.4f, sig^2: %.4f, rate_pred: %.4f, v_mean_pred: %.4f" % (
            self, self.tau_eff, self.mu, self.sigma_square, self.rate_prediction, self.v_mean_prediction))
        for s in self.sources:
            print("\t\t", s.print_sys())























class MFOldPop(object):
    """pop: similar neurons"""

    def __init__(self, name, params):
        self.name = name
        self.sources = []
        self.params = params
        self.noise = None

        self.n = 1
        # estimation
        self.rate_ = 0.  # spikes/s
        self.v_mean = -60.  # pop mean voltage
        self.is_adapting = False

    @property
    def rate_hz(self):
        return self.rate_

    @rate_hz.setter
    def rate_hz(self, value):
        self.rate_ = value

    @property
    def rate_ms(self):
        return self.rate_ / 1e3

    @rate_ms.setter
    def rate_ms(self, value):
        self.rate_ = value * 1e3

    @property
    def has_nmda(self):
        return any([s.is_nmda for s in self.sources])

    @property
    def J(self):
        """
        Linearization factor for NMDA
        Unitless
        """
        return 1 + self.params["gamma"] * np.exp(-self.params["beta"] * self.v_mean)

    @property
    def total_cond(self):
        """
        Gm * SE in [1]
        Units of S
        """
        return self.params["g_L"] + np.sum(s.conductance for s in self.sources)

    @property
    def tau_eff(self):
        """
        Seconds
        """
        return self.params["C_m"] / self.total_cond

    @property
    def mu(self):
        """
        Volt
        """
        return 1. / self.total_cond * np.sum(s.voltage_conductance for s in self.sources)

    @property
    def sigma_square(self):
        """
        Volt^2
        """
        if not self.noise:
            return 0.
        return (self.noise.g_base / self.params["C_m"] * (self.v_mean - self.noise.E_rev)) ** 2 * self.tau_eff * self.noise.g_dyn() * self.noise.noise_tau

    @property
    def rate_prediction(self):
        return 1e3 * self.phi_firing_func()

    @property
    def v_mean_prediction(self):
        """
        Volt
        """
        return self.params["E_L"] + self.mu - (self.params["V_th"] - self.params["V_reset"]) * self.rate_ms * self.tau_eff

    def phi_firing_func(self):

        sigma = np.sqrt(self.sigma_square)
        tau_eff = self.tau_eff

        beta = (self.params["V_reset"] - self.params["E_L"] - self.mu) / sigma
        alpha = -0.5 * self.noise.noise_tau / tau_eff \
                + 1.03 * np.sqrt(self.noise.noise_tau / tau_eff) \
                + (-self.mu - self.params["E_L"] + self.params["V_th"]) * (1. + (0.5 * self.noise.noise_tau / tau_eff)) / sigma

        def integrand(x):
            if x < -10.:
                return np.exp(10. ** 2) * (1. + erf(10.))
            if x > 10.:
                return 0.
            return np.exp(x ** 2) * (1. + erf(x))
        return 1. / (self.params["t_ref"] + tau_eff * np.sqrt(np.pi) * quad(integrand, beta, alpha)[0])

    def __repr__(self):
        return "MFpop [%s] <%s (%i sources, n: %i, rate: %.4f, v_mean: %.4f)>" % (id(self), self.name, len(self.sources), self.n, self.rate_hz, self.v_mean)

    def print_sys(self, mf=False):
        print("\t%s - tau_eff: %.1fms, mu: %.4f, sig^2: %.4f, rate_pred: %.4f, v_mean_pred: %.4f" % (self, self.tau_eff, self.mu, self.sigma_square, self.rate_prediction, self.v_mean_prediction))
        for s in self.sources:
            print("\t\t", s.print_sys())
