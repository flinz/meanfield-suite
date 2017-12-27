from abc import abstractmethod
from types import MappingProxyType
from typing import Union, Dict

from brian2 import units, check_units, Equations

from meanfield.parameters.MFParams import MFParams
from meanfield.utils import create_identifier, lazyproperty, create_name
from meanfield.parameters import IP, PP
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.parameters.Connection import ConnectionStrategy
from meanfield.parameters import Connection


class MFInput(object):
    """Source: a synapse coupled to pops"""

    arguments = MappingProxyType({
        IP.GM: units.siemens,
        IP.VREV: units.volt,
        IP.TAU: units.second,  # TODO : tau subclass?
    })

    defaults = MappingProxyType({})

    def __init__(self, target: MFPopulation, parameters: Union[Dict, MFParams]=None, name: str=None, connection: ConnectionStrategy=Connection.all_to_all(), add_as_input: bool=True):

        self.name = name if name else create_name(self)
        self.ref = create_identifier(self.name)

        self.parameters = MFParams({}) if not parameters else MFParams(parameters)
        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

        self.g_base = parameters[IP.GM]  # TODO : parametrize, solve for ? => general setter checking for units

        self.connection = connection

        # link to pop
        self.target = target
        if add_as_input:
            self.target.add_input(self)  # TODO consistent

    def __getitem__(self, key):
        return self.parameters[key]

    def __repr__(self):
        return "{} [{}] ({}, {})".format(self.__class__.__name__, self.name, self.parameters, self.connection)

    # Theory

    @property
    @check_units(result=units.siemens)
    def conductance(self):
        return self.g_dyn() * self.g_base

    @property
    @check_units(result=units.amp)
    def voltage_conductance(self):
        return self.conductance * (self[IP.VREV] - self.target[PP.VL])

    @abstractmethod
    def g_dyn(self):
        pass

    # Simulation

    @abstractmethod
    def brian2(self):
        """Builds lazily Brian2 synapse component once."""
        pass

    def brian2_model(self):
        """Returns Brian2 dynamic (Equations) affecting specified populations."""
        return Equations(
            '''
            I = g * (v - vrev) * s : amp
            ds / dt = - s / tau : 1
            ''',
            I=self.current_name,
            g=self[IP.GM],
            s=self.post_variable_name,
            vrev=self[IP.VREV],
            tau=self[IP.TAU],
        )

    @property
    def current_name(self):
        return 'I_' + self.ref

    @property
    def post_variable_name(self):
        return 's_' + self.ref

