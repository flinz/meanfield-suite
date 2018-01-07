from types import MappingProxyType
from typing import Union, Dict

from brian2 import check_units, Equations, Synapses, BrianObject, units

from meanfield.inputs.MFInput import MFInput
from meanfield.parameters import IP
from meanfield.parameters.MFParameters import MFParameters
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import lazyproperty


class MFLinearInput(MFInput):

    arguments = MappingProxyType({
        IP.GM: units.siemens,
        IP.VREV: units.volt,
        IP.TAU: units.second,
        IP.W: 1,
    })

    defaults = MappingProxyType({
        IP.W: 1,
    })

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParameters], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

        target.add_input(self)

    # Theory

    @check_units(result=1)
    def g_dyn(self):
        return self.connection.theory(self.origin.n) * self.origin.rate * self[IP.TAU] * self[IP.W]

    # Simulation

    @lazyproperty
    def brian2(self) -> BrianObject:
        model = Equations('w : 1')
        on_pre = f'{self.post_variable_name} += w'
        syn = Synapses(
            source=self.origin.brian2,
            target=self.target.brian2,
            method='euler',
            model=model,
            on_pre=on_pre,
            name=self.ref
        )
        self.connection.simulation(syn)
        syn.w[:] = self[IP.W]
        return syn

