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
