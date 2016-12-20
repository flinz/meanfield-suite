from brian2 import units, Equations

from MFParams import MFParams
from params import SP, NP


class MFSource(object):
    """Source: a synapse coupled to pops"""

    def __init__(self, name, pop):

        self.name = name
        self.g_dyn = lambda: 0.  # this should be given as a function describing the synaptic conductance
        self.is_nmda = False
        self.g_base = 0.  # [nS]
        self.E_rev = 0.   # [mV] excitatory by default
        self.noise_tau = 0.   # [mV] excitatory by default

        self.params = MFParams({
            SP.GM: self.g_base * units.nS,
            SP.VE: self.E_rev * units.mvolt,
            SP.TAU_M: 10 * units.msecond
        })
        self.params.verify({
            SP.GM: units.siemens,
            SP.VE: units.volt,
            SP.TAU_M: units.second
        })

        # link to pop
        pop.sources.append(self)
        self.pop = pop

    def brian2_model(self, n, var='I'):
        return Equations(
            '''
            I = g * (v - ve) * s : amp
            ds / dt = - s / tau : 1
            ''',
            s='s'+n,
            I=var,
            g=self.params[SP.GM],
            ve=self.params[SP.VE],
            tau=self.params[SP.TAU_M]
        )

    @property
    def conductance(self):
        tmp = self.g_dyn() * self.g_base
        if self.is_nmda:
            J_ = 1./self.pop.J
            return tmp * J_ * (
                1. + (1.-J_) * self.pop.params[NP.BETA] * (self.pop.v_mean - self.E_rev)
            )
        return tmp

    @property
    def voltage_conductance(self):
        cond = self.g_dyn() * self.g_base
        tmp = cond * (self.E_rev - self.pop.params[NP.VL] / units.volt)
        if self.is_nmda:
            J_ = 1./self.pop.J
            return J_ * (
                tmp + (1.-J_) * cond * self.pop.params[NP.BETA] * (self.pop.v_mean - self.E_rev) * (self.pop.v_mean - self.pop.params["E_L"])
            )
        return tmp

    def print_sys(self):
        print("%s: %f" % (self, self.conductance))

    def __repr__(self):
        return "MFSource [%s] <%s, nmda: %s, E_rev: %.1f>" % (id(self),self.name, self.is_nmda, self.E_rev)
