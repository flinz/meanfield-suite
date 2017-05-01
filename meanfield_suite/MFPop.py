from abc import abstractproperty, abstractmethod
from math import erf

import numpy as np
from brian2 import units, Equations, NeuronGroup, check_units
from scipy.integrate import quad

from MFParams import MFParams
from Utils import lazy
from params import NP, SP


class MFPop(object):

    def __init__(self, name, n, params=None):
        self.name = name
        self.n = n
        self.params = MFParams({}) if params is None else MFParams(params)

        self.sources = []
        self.noise = None
        # base estimation
        self._rate = 0. * units.Hz
        self._v_mean = -60. * units.mV

    def add_noise(self, noise):
        self.noise = noise

    def add_source(self, source):
        self.sources.append(source)

    @property
    @check_units(result=units.Hz)
    def rate(self):
        return self._rate

    @rate.setter
    @check_units(value=units.Hz)
    def rate(self, value):
        self._rate = value

    @property
    @check_units(result=units.volt)
    def v_mean(self):
        return self._v_mean

    @v_mean.setter
    @check_units(value=units.volt)
    def v_mean(self, value):
        self._v_mean = value

    #@property
    #@abstractmethod
    @abstractproperty
    def rate_prediction(self):
        pass

    #@property
    #@abstractmethod
    @abstractproperty
    def v_mean_prediction(self):
        pass


class MFLinearPop(MFPop):
    """pop: similar neurons"""

    def __init__(self, name, n, params):
        super().__init__(name, n, params)

        defaults = {}
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
        for i, s in enumerate(self.sources):
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

    @lazy
    def brian2(self, method='euler'):
        P = NeuronGroup(
            self.n,
            self.brian2_model(),
            method=method,
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
        return self.params[NP.GM] + np.sum(s.conductance for s in self.sources)

    @property
    @check_units(result=units.second)
    def tau_eff(self):
        """
        Seconds
        """
        print(self.params[NP.CM])
        print(self.total_cond)
        return self.params[NP.CM] / self.total_cond

    @property
    @check_units(result=units.volt)
    def mu(self):
        """
        Volt
        """
        return 1. / self.total_cond * np.sum(s.voltage_conductance for s in self.sources)

    @property
    @check_units(result=units.volt ** 2)
    def sigma_square(self):
        """
        Volt^2
        """
        if not self.noise:
            return 0. * units.volt
        return (self.noise.g_base / self.params[NP.CM] * (self.v_mean - self.noise.params[SP.VE])) ** 2 * self.tau_eff * self.noise.g_dyn() * self.noise.params[SP.TAU]

    def phi_firing_func(self):

        print(self.sigma_square)
        sigma = np.sqrt(self.sigma_square)
        tau_eff = self.tau_eff

        beta = (self.params[NP.VRES] - self.params[NP.VL] - self.mu) / sigma
        alpha = -0.5 * self.noise.params[SP.TAU] / tau_eff \
                + 1.03 * np.sqrt(self.noise.params[SP.TAU] / tau_eff) \
                + (-self.mu - self.params[NP.VL] + self.params[NP.VTHR]) * (
                1. + (0.5 * self.noise.params[SP.TAU] / tau_eff)) / sigma

        def integrand(x):
            if x < -10.:
                return np.exp(10. ** 2) * (1. + erf(10.))
            if x > 10.:
                return 0.
            return np.exp(x ** 2) * (1. + erf(x))

        return 1. / (self.params[NP.TAU_RP] + tau_eff * np.sqrt(np.pi) * quad(integrand, beta, alpha)[0])

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

    def __repr__(self):
        return "MFpop [{}] <{} ({} sources, n: {}, rate: {}, v_mean: {})>".format(id(self), self.name, len(self.sources), self.n, self.rate, self.v_mean)

    def print_sys(self):
        print("\t{} - tau_eff: {}, mu: {}, sig^2: {}, rate_pred: {}, v_mean_pred: {}".format(
            self, self.tau_eff, self.mu, self.sigma_square, self.rate_prediction, self.v_mean_prediction))
        for s in self.sources:
            print("\t\t", s.print_sys())


class MFNonLinearPop(MFLinearPop):
    """pop: similar neurons"""

    def __init__(self, name, n, params):
        super().__init__(name, n, params)

        defaults = {}
        expectations = {
            NP.GAMMA: 1,
            NP.BETA: 1,
        }

        self.params.fill(defaults)
        self.params.verify(expectations)

    @property
    def J(self):
        """Linearization factor for NMDA"""
        return 1 + self.params[NP.GAMMA] * np.exp(-self.params[NP.BETA] * self.v_mean)


