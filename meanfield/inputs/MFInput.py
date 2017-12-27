from abc import abstractmethod
from typing import Union, Dict

from brian2 import units, check_units, Equations

from meanfield.parameters.MFParams import MFParams
from meanfield.utils import name2identifier, lazyproperty
from meanfield.parameters import IP, PP
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.parameters.Connection import ConnectionStrategy


class MFInput(object):
    """Source: a synapse coupled to pops"""

    def __init__(self, name: str, pop: MFPopulation, params: Union[Dict, MFParams], connection: ConnectionStrategy, add_as_input=True):

        self.name = name
        self.ref = name2identifier(name)
        self.params = MFParams(params)

        defaults = {}
        expectations = {
            IP.GM: units.siemens,
            IP.VREV: units.volt,
            IP.TAU: units.second,  # TODO : tau subclass?
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

        self.g_base = params[IP.GM]  # TODO : parametrize, solve for ? => general setter checking for units

        self.connection = connection

        # link to pop
        self.pop = pop
        if add_as_input:
            self.pop.add_input(self)  # TODO consistent

    def __repr__(self):
        return "{} [{}] ({}, {})".format(self.__class__.__name__, self.name, self.params, self.connection)

    # Theory

    @property
    @check_units(result=units.siemens)
    def conductance(self):
        return self.g_dyn() * self.g_base

    @property
    @check_units(result=units.amp)
    def voltage_conductance(self):
        return self.conductance * (self.params[IP.VREV] - self.pop.params[PP.VL])

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
            g=self.params[IP.GM],
            s=self.post_variable_name,
            vrev=self.params[IP.VREV],
            tau=self.params[IP.TAU],
        )

    @property
    def current_name(self):
        return 'I_' + self.ref

    @property
    def post_variable_name(self):
        return 's_' + self.ref

