from abc import abstractmethod
from typing import Union, Dict

from brian2 import units, check_units, Equations

from meanfield.parameters.MFParams import MFParams
from meanfield.utils import name2identifier, lazyproperty
from meanfield.parameters import SP, NP
from meanfield.populations.MFPop import MFPop


class MFSource(object):
    """Source: a synapse coupled to pops"""

    def __init__(self, name: str, pop: MFPop, params: Union[Dict, MFParams]):

        self.name = name
        self.ref = name2identifier(name)
        self.params = MFParams(params)

        defaults = {}
        expectations = {
            SP.GM: units.siemens,
            SP.VREV: units.volt,
            SP.TAU: units.second,  # TODO : tau subclass?
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

        self.g_base = params[SP.GM]  # TODO : parametrize, solve for ? => general setter checking for units

        # link to pop
        self.pop = pop
        self.pop.add_source(self)  # TODO consistent

    @property
    @check_units(result=units.siemens)
    def conductance(self):
        return self.g_dyn() * self.g_base

    @property
    @check_units(result=units.amp)
    def voltage_conductance(self):
        return self.conductance * (self.params[SP.VREV] - self.pop.params[NP.VL])

    @abstractmethod
    def g_dyn(self):
        pass

    @property
    def current_name(self):
        return 'I_' + self.ref

    @property
    def post_variable_name(self):
        return 's_' + self.ref

    @abstractmethod
    def b2_syn(self):
        """Builds lazily Brian2 synapse component once."""
        pass

    def b2_dyn(self):
        """Returns Brian2 dynamic (Equations) affecting specified populations."""
        return Equations(
            '''
            I = g * (v - vrev) * s : amp
            ds / dt = - s / tau : 1
            ''',
            I=self.current_name,
            g=self.params[SP.GM],
            s=self.post_variable_name,
            vrev=self.params[SP.VREV],
            tau=self.params[SP.TAU],
        )

    def __repr__(self):
        return "{} [{}] ({})".format(self.__class__.__name__, self.name, self.params)

