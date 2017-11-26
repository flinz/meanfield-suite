from abc import abstractproperty, abstractmethod

from brian2 import units, check_units

from meanfield.parameters.MFParams import MFParams


class MFPop(object):

    def __init__(self, name, n, params=None):
        self.name = name
        self.n = n
        if params:
            if isinstance(params, MFParams):
                self.params = params
            else:
                self.params = MFParams(params)
        else:
            self.params = MFParams({})

        self.sources = []
        self.noise = []

        # base estimation
        self._rate = 0. * units.Hz
        self._v_mean = -60. * units.mV

    def add_noise(self, noise):
        if len(self.noise):
            raise NotImplementedError('multiple noise not supported yet')

        self.noise.append(noise)

    def add_source(self, source):
        # TODO control source and no duplicate
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

    @abstractmethod
    def rate_prediction(self):
        pass

    @abstractmethod
    def v_mean_prediction(self):
        pass

    @abstractmethod
    def brian2(self):
        pass
