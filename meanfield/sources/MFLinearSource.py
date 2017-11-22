from brian2 import check_units, Equations, Synapses

from meanfield.parameters.Connection import ConnectionStrategy
from meanfield.parameters import SP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPop import MFPop
from meanfield.sources.MFSource import MFSource
from meanfield.utils import lazyproperty
from meanfield.parameters import Connection


class MFLinearSource(MFSource):

    def __init__(self, name: str, pop: MFPop, params: MFParams, from_pop: MFPop, connection: ConnectionStrategy=Connection.all_to_all()):
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
        self.connection = connection

    @check_units(result=1)
    def g_dyn(self):
        return self.connection.theory(self.from_pop.n) * self.from_pop.rate * self.params[SP.TAU] * self.params[SP.W]

    @lazyproperty
    def b2_syn(self, method='euler', weight=1, **kv):
        model = Equations('w : 1')
        on_pre = '{} += w'.format(self.post_variable_name)
        syn = Synapses(
            source=self.from_pop.brian2,
            target=self.pop.brian2,
            method=method,
            model=model,
            on_pre=on_pre
        )
        self.connection.simulation(syn)
        syn.w[:] = self.params[SP.W]
        return syn
