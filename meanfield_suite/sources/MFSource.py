from abc import abstractmethod

from brian2 import units, Equations, check_units

from MFParams import MFParams
from Utils import name2identifier
from params import SP, NP


class MFSource(object):
    """Source: a synapse coupled to pops"""
    # abstract

    def __init__(self, name, pop, params):

        self.name = name
        self.ref = name2identifier(name)
        self.g_base = 0. * units.siemens  # [nS]

        defaults = {}
        expectations = {
            SP.GM: units.siemens,
            SP.VE: units.volt,
            SP.TAU: units.second,
        }

        self.params = MFParams(params)
        self.params.fill(defaults)
        self.params.verify(expectations)

        self.g_base = params[SP.GM]  # TODO : parameterize, solve for ? => general setter checking for units

        # link to pop
        self.pop = pop
        self.pop.sources.append(self)  # TODO consistent

    @property
    @check_units(result=units.siemens)
    def conductance(self):
        return self.g_dyn() * self.g_base

    @property
    def voltage_conductance(self):
        cond = self.g_dyn() * self.g_base
        return cond * (self.params[SP.VE] - self.pop.params[NP.VL])


    def print_sys(self):
        print("{}: {}".format(self, self.conductance))

    def __repr__(self):
        return "MFSource [{}] <{}, nmda: {}, E_rev: {}>".format(id(self), self.name, self.is_nmda, self.params[SP.VE])

    @abstractmethod
    def g_dyn(self):
        pass

    @abstractmethod
    def brian2(self):
        pass

    @property
    def current_name(self):
        return 'I_' + self.ref

    @property
    def post_variable_name(self):
        return 's_' + self.ref

    def brian2_model(self):
        return Equations(
            '''
            I = g * (v - ve) * s : amp
            ds / dt = - s / tau : 1
            ''',
            s=self.post_variable_name,
            I=self.current_name,
            g=self.params[SP.GM],
            ve=self.params[SP.VE],
            tau=self.params[SP.TAU]
        )