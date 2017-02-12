from abc import abstractmethod

from brian2 import units, Equations, Synapses, PoissonInput

from MFParams import MFParams
from params import SP, NP

class MFSource(object):
    """Source: a synapse coupled to pops"""

    def __init__(self, name, pop, params):

        self.name = name
        self.is_nmda = False
        self.g_base = 0.  # [nS]
        self.E_rev = 0.   # [mV] excitatory by default
        self.noise_tau = 0.   # [mV] excitatory by default

        defaults = {}
        expectations = {
            SP.GM: units.siemens,
            SP.VE: units.volt,
            SP.TAU_M: units.second,
        }
        # TODO : which is static/dynamic

        self.params = MFParams(params)
        self.params.fill(defaults)
        self.params.verify(expectations)

        # TODO : linearize with unit or not ?
        self.g_base = params[SP.GM] / units.siemens
        self.E_rev = params[SP.VE] / units.volt
        self.noise_tau = params[SP.TAU_M] / units.second

        # link to pop
        self.pop = pop
        self.pop.sources.append(self) # TODO consistent

    @property
    def conductance(self):
        tmp = self.g_dyn() * self.g_base
        #if self.is_nmda:
        #    J_ = 1./self.pop.J
        #    return tmp * J_ * (
        #        1. + (1.-J_) * self.pop.params[NP.BETA] * (self.pop.v_mean - self.E_rev)
        #    )
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

    @abstractmethod
    def g_dyn(self):
        pass

    @abstractmethod
    def brian2(self):
        pass

    def brian2_model(self, n):
        current = 'I{}'.format(n)
        return current, Equations(
            '''
            I = g * (v - ve) * s : amp
            ds / dt = - s / tau : 1
            ''',
            s='s_{}'.format(self.name.replace(' ', '_')),
            I='I_{}'.format(self.name.replace(' ', '_')),
            g=self.params[SP.GM],
            ve=self.params[SP.VE],
            tau=self.params[SP.TAU_M]
        )

    # TODO store post-corresponding var s_{}


class MFStaticSource(MFSource):

    def __init__(self, name, pop, params, rate, n):
        super().__init__(name, pop, params)
        self.rate = rate
        self.n = n
        # TODO real param

    def g_dyn(self):
        # TODO : unit return ?
        return self.rate * self.n * self.params[SP.TAU_M] / units.second

    @property
    def brian2(self):
        return None


class MFDynamicSource(MFSource):

    def __init__(self, name, pop, params, from_pop):
        super().__init__(name, pop, params)
        self.from_pop = from_pop

        defaults = {
            SP.W: 1.,
            SP.FRAC: 1.
        }
        expectations = {
            SP.W: 1, # unitless
            SP.FRAC: 1
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    def g_dyn(self):
        # TODO : unit return ?
        return self.from_pop.n * self.from_pop.rate_ms * self.params[SP.W] * self.params[SP.FRAC] * self.params[SP.TAU_M] / units.second

    @property
    def brian2(self, mode='i != j'):
        model = Equations('w : 1')
        eqs_pre = '''
        s_AMPA += w
        s_NMDA += w
        '''
        C = Synapses(self.from_pop.brian2, self.pop.brian2, method='euler', model=model, on_pre=eqs_pre)
        C.connect(mode)
        C.w[:] = 1
        return C

