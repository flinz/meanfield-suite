from types import MappingProxyType
from typing import Dict, Union

from brian2 import Equations, Synapses, check_units, BrianObject

from meanfield.inputs.MFLinearNMDAInput import MFLinearNMDAInput
from meanfield.utils import lazyproperty
from meanfield.parameters import IP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.parameters import Connection
from meanfield.parameters.Connection import ConnectionStrategy


class MFNonLinearNMDAInput(MFLinearNMDAInput):

    arguments = MappingProxyType({
        IP.TAU_NMDA: 1.,
        IP.TAU_NMDA_RISE: 1.,
        IP.BETA: 1.,
    })

    defaults = MappingProxyType({})

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParams], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    # Simulation

    @property
    def post_nonlinear_name(self):
        return 'x_' + self.ref

    @lazyproperty
    def brian2(self) -> BrianObject:
        model = Equations(
            '''
            w : 1
            s_post = w * s : 1 (summed)
            ds / dt = - s / tau_decay + alpha * x * (1 - s) : 1 (clock-driven)
            dx / dt = - x / tau_rise : 1 (clock-driven)
            ''',
            s_post=self.post_variable_name + '_post',
            s=self.post_variable_name,
            x=self.post_nonlinear_name,
            tau_decay=self[IP.TAU_NMDA],
            tau_rise=self[IP.TAU_NMDA_RISE],
            alpha=self[IP.ALPHA], # TODO ALPHA ? 1 / ms
        )
        eqs_pre = '''
        {} += 1
        '''.format(self.post_nonlinear_name)
        C = Synapses(
            self.origin.brian2,
            self.target.brian2,
            method='euler',
            model=model,
            on_pre=eqs_pre,
            name=self.ref,
        )
        C.connect()
        C.w[:] = 1
        return C

    @property
    def brian2_model(self) -> Equations:
        return Equations(
            '''
            I = g * (v - ve) / (1 + gamma * exp(- beta * v) ) * s_post : amp
            s_post: 1
            ''',
            s_post=self.post_variable_name + '_post',
            I=self.current_name,
            g=self[IP.GM],
            ve=self[IP.VE],
            beta=self[IP.BETA]
        )
