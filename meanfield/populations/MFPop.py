from abc import abstractproperty, abstractmethod

from brian2 import units, check_units

from meanfield.parameters.MFParams import MFParams
from meanfield.utils import name2identifier


class MFPop(object):

    def __init__(self, name, n, params=None):
        self.name = name
        self.ref = name2identifier(name)

        self.n = n
        if params:
            if isinstance(params, MFParams):
                self.params = params
            else:
                self.params = MFParams(params)
        else:
            self.params = MFParams({})

        self.sources = []
        self.noises = []

        # base estimation
        self._rate = 200 * units.Hz
        self._v_mean = -55. * units.mV

    def add_noise(self, noise):
        if len(self.noises):
            raise NotImplementedError('multiple noise not supported yet')

        self.noises.append(noise)

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

    @abstractmethod
    def brian2_model(self):
        pass

    def introspect(self, indent=0) -> str:
        builder = []
        spaces = ' ' * indent
        builder.append(str(self))

        for noise in self.noises:
            builder.append('{}  - {}'.format(spaces, noise))

        for source in self.sources:
            builder.append('{}  - {}'.format(spaces, source))

        return '\n'.join(builder)


    def __repr__(self):
        return '{} [{}] ({} noises, {} sources, n: {}, rate: {}, v_mean: {})'.format(self.__class__.__name__, self.name, len(self.noises), len(self.sources), self.n, self.rate, self.v_mean)
