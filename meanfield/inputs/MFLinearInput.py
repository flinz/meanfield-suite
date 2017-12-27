from typing import Union, Dict

from brian2 import check_units, Equations, Synapses

from meanfield.parameters import Connection
from meanfield.parameters import IP
from meanfield.parameters.Connection import ConnectionStrategy
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.inputs.MFInput import MFInput
from meanfield.utils import lazyproperty


class MFLinearInput(MFInput):

    def __init__(self, name: str, pop: MFPopulation, params: Union[Dict, MFParams], from_pop: MFPopulation, connection: ConnectionStrategy=Connection.all_to_all()):
        super().__init__(name, pop, params, connection)

        defaults = {
            IP.W: 1,
        }
        expectations = {
            IP.W: 1,
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

        self.from_pop = from_pop

    # Theory

    @check_units(result=1)
    def g_dyn(self):
        return self.connection.theory(self.from_pop.n) * self.from_pop.rate * self.params[IP.TAU] * self.params[IP.W]

    # Simulation

    @lazyproperty
    def brian2(self):
        model = Equations('w : 1')
        on_pre = '{} += w'.format(self.post_variable_name)
        syn = Synapses(
            source=self.from_pop.brian2,
            target=self.pop.brian2,
            method='euler',
            model=model,
            on_pre=on_pre,
            name=self.ref
        )
        self.connection.simulation(syn)
        syn.w[:] = self.params[IP.W]
        return syn

