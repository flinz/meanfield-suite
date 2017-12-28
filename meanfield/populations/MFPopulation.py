from abc import abstractmethod
from types import MappingProxyType
from typing import Union, Optional

from brian2 import units, check_units, Equations, BrianObject

from meanfield.parameters.MFParams import MFParams
from meanfield.utils import create_identifier, create_name


class MFPopulation(object):

    arguments = MappingProxyType({})

    defaults = MappingProxyType({})

    @check_units(n=1)
    def __init__(self, n: int, parameters: Optional[Union[dict, MFParams]]=None, name: str=None):

        self.name = name if name else create_name(self)
        self.ref = create_identifier(self.name)
        self.n = n

        self.parameters = MFParams({}) if not parameters else MFParams(parameters)
        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

        self.inputs = []
        self.noises = []

        # base estimation
        self._rate = 200 * units.Hz
        self._v_mean = -55. * units.mV

    def __getitem__(self, key):
        return self.parameters[key]

    def add_noise(self, noise) -> 'MFPopulation':
        if len(self.noises):
            raise NotImplementedError('multiple noise not supported yet')

        self.noises.append(noise)
        return self

    def add_input(self, input) -> 'MFPopulation':
        # TODO control source and no duplicate
        self.inputs.append(input)
        return self

    @property
    @check_units(result=units.Hz)
    def rate(self) -> units.Hz:
        return self._rate

    @rate.setter
    @check_units(value=units.Hz)
    def rate(self, value) -> units.Hz:
        self._rate = value

    @property
    @check_units(result=units.volt)
    def v_mean(self) -> units.volt:
        return self._v_mean

    @v_mean.setter
    @check_units(value=units.volt)
    def v_mean(self, value) -> units.volt:
        self._v_mean = value

    def introspect(self, indent=0) -> str:
        builder = []
        spaces = ' ' * indent
        builder.append(str(self))

        for noise in self.noises:
            builder.append('{}  - {}'.format(spaces, noise))

        for inputs in self.inputs:
            builder.append('{}  - {}'.format(spaces, inputs))

        return '\n'.join(builder)

    def __repr__(self) -> str:
        return '{} [{}] ({} noises, {} inputs, n: {}, rate: {}, v_mean: {})'.format(self.__class__.__name__, self.name, len(self.noises), len(self.inputs), self.n, self.rate, self.v_mean)

    # Theory

    @abstractmethod
    @check_units(result=units.Hz)
    def rate_prediction(self) -> units.Hz:
        pass

    @abstractmethod
    @check_units(result=units.volt)
    def v_mean_prediction(self) -> units.volt:
        pass

    # Simulation

    @abstractmethod
    def brian2(self) -> BrianObject:
        pass

    @abstractmethod
    def brian2_model(self) -> Optional[Equations]:
        pass

