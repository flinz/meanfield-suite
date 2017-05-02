from brian2 import check_units, Equations, Synapses

from MFSource import MFSource
from Utils import lazy


class MFDynamicSource(MFSource):

    def __init__(self, name, pop, params, from_pop, synapse=None):
        super().__init__(name, pop, params)
        self.from_pop = from_pop
        self.synapse = synapse
        defaults = {
            SP.W: 1.,
            SP.FRAC: 1.
        }
        expectations = {
            SP.W: 1., # unitless
            SP.FRAC: 1.
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    @check_units(result=1)
    def g_dyn(self):
        activation = self.synapse(self.from_pop.rate) if self.synapse else self.from_pop.rate * self.params[SP.TAU]
        return self.from_pop.n * activation * self.params[SP.W]

    @lazy
    def brian2(self, mode='i != j'):
        model = Equations('w : 1')
        eqs_pre = '''
        {} += w
        '''.format(self.post_variable_name)
        C = Synapses(self.from_pop.brian2, self.pop.brian2, method='euler', model=model, on_pre=eqs_pre)
        C.connect(mode)
        C.w[:] = 1
        return C

    def __repr__(self):
        return "MFSource [{}] <{}, nmda: {}, E_rev: {}>".format(id(self), self.name, self.is_nmda, self.params[SP.VE])