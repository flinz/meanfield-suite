from abc import abstractmethod
from types import MappingProxyType
from typing import Union, Optional

from brian2 import units, check_units, Equations, BrianObject

from meanfield.parameters.MFParameters import MFParameters
from meanfield.utils import create_identifier, create_name


class MFPopulation(object):

    arguments = MappingProxyType({})

    defaults = MappingProxyType({})

    @check_units(n=1)
    def __init__(self, n: int, parameters: Optional[Union[dict, MFParameters]]=None, name: str=None):

        self.name = name if name else create_name(self)
        self.ref = create_identifier(self.name)
        self.n = n

        self.parameters = MFParameters({}) if not parameters else MFParameters(parameters)
        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

        self.inputs = []
        self.noises = []

        # base estimation
        self.rate = 200 * units.Hz

    def __getitem__(self, key):
        return self.parameters[key]

    def __setitem__(self, key, value):
        self.parameters[key] = value

    def add_noise(self, noise) -> 'MFPopulation':
        if len(self.noises):
            raise NotImplementedError('multiple noise not supported yet')

        self.noises.append(noise)
        return self

    def add_input(self, input) -> 'MFPopulation':
        self.inputs.append(input)
        return self

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
        return '{} [{}] ({} noises, {} inputs, n: {}, rate: {})'.format(self.__class__.__name__, self.name, len(self.noises), len(self.inputs), self.n, self.rate)

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

