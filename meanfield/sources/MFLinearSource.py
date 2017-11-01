from brian2 import check_units, Equations, Synapses

from meanfield.sources.MFSource import MFSource
from meanfield.utils import lazyproperty
from meanfield.parameters import SP


class MFLinearSource(MFSource):

    def __init__(self, name, pop, params, from_pop):
        super().__init__(name, pop, params)

        defaults = {
            SP.W: 1,
        }
        expectations = {
            SP.W: 1,
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

        self.from_pop = from_pop

    @check_units(result=1)
    def g_dyn(self):
        return self.from_pop.n * self.from_pop.rate * self.params[SP.TAU] * self.params[SP.W]

    @lazyproperty
    def b2_syn(self, method='euler', weight=1, **kv):
        model = Equations('w : 1')
        on_pre = '{} += w'.format(self.post_variable_name)
        print(self.from_pop)
        syn = Synapses(
            source=self.from_pop.brian2,
            target=self.pop.brian2,
            method=method,
            model=model,
            on_pre=on_pre
        )
        syn.connect(j='i')
        syn.w[:] = weight
        return syn
