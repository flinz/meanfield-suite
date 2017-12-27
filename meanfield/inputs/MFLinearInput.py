from types import MappingProxyType
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

    arguments = MappingProxyType({
        IP.W: 1,
    })

    defaults = MappingProxyType({
        IP.W: 1,
    })

    def __init__(self, name: str, pop: MFPopulation, parameters: Union[Dict, MFParams], from_pop: MFPopulation, connection: ConnectionStrategy=Connection.all_to_all()):
        super().__init__(name, pop, parameters, connection)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

        self.from_pop = from_pop

    # Theory

    @check_units(result=1)
    def g_dyn(self):
        return self.connection.theory(self.from_pop.n) * self.from_pop.rate * self[IP.TAU] * self[IP.W]

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
        syn.w[:] = self[IP.W]
        return syn

