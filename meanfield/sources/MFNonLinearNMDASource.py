from typing import Dict, Union

from brian2 import Equations, Synapses, check_units

from meanfield.sources.MFLinearNMDASource import MFLinearNMDASource
from meanfield.utils import lazyproperty
from meanfield.parameters import SP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPop import MFPop
from meanfield.parameters import Connection
from meanfield.parameters.Connection import ConnectionStrategy


class MFNonLinearNMDASource(MFLinearNMDASource):

    def __init__(self, name: str, pop: MFPop, params: Union[Dict, MFParams], from_pop: MFPop, connection: ConnectionStrategy=Connection.all_to_all()):
        super().__init__(name, pop, params, from_pop, connection)
        defaults = {

        }
        expectations = {
            SP.TAU_NMDA: 1.,
            SP.TAU_NMDA_RISE: 1.,
            SP.BETA: 1.,
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    def brian2_model(self):
        return Equations(
            '''
            I = g * (v - ve) / (1 + gamma * exp(- beta * v) ) * s_post : amp
            s_post: 1
            ''',
            s_post=self.post_variable_name + '_post',
            I=self.current_name,
            g=self.params[SP.GM],
            ve=self.params[SP.VE],
            beta=self.params[SP.BETA]
        )

    @property
    def post_nonlinear_name(self):
        return 'x_' + self.ref

    @lazyproperty
    def brian2(self):
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
            tau_decay=self.params[SP.TAU_NMDA],
            tau_rise=self.params[SP.TAU_NMDA_RISE],
            alpha=self.params[SP.ALPHA], # TODO ALPHA ? 1 / ms
        )
        eqs_pre = '''
        {} += 1
        '''.format(self.post_nonlinear_name)
        C = Synapses(self.from_pop.brian2, self.pop.brian2, method='euler', model=model, on_pre=eqs_pre)
        C.connect()
        C.w[:] = 1
        return C
