from types import MappingProxyType
from typing import Union, Dict

from brian2 import check_units, Equations, Synapses, BrianObject, units

from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.inputs.MFLinearNMDAInput import MFLinearNMDAInput
from meanfield.inputs.modifiers.stp import STP
from meanfield.inputs.MFInput import MFInput
from meanfield.parameters import IP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import lazyproperty


class MFLinearSTPInput(MFLinearInput):
    arguments = MappingProxyType({})

    defaults = MappingProxyType({})

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParams], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    @check_units(result=1)
    def g_dyn(self):
        return self.connection.theory(self.origin.n) * self.synapse(self.origin.rate) * self[IP.W]

    # Simulation

    @lazyproperty
    def brian2(self) -> BrianObject:
        model = Equations('''
        x : 1
        u : 1
        w : 1
        ''')

        U = self.synapse.U
        tauf = self.synapse.tau_f
        taud = self.synapse.tau_x

        on_pre = f'''
        u = {U} + (u - {U}) * exp(- (t - lastupdate) / ({tauf / units.ms} * ms))
        x = 1 + (x - 1) * exp(- (t - lastupdate) / ({taud / units.ms} * ms))

        {self.post_variable_name} += w * u * x

        x *= (1 - u)
        u += {U} * (1 - u)
        '''
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
        syn.x[:] = 1
        syn.u[:] = 1
        return syn

