from typing import Union, Dict

from brian2 import Equations, Synapses, units

from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.utils import lazyproperty
from meanfield.parameters import IP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPop import MFPop
from meanfield.parameters import Connection
from meanfield.parameters.Connection import ConnectionStrategy


class MFLinear3TSInput(MFLinearInput):
    def __init__(self, name: str, pop: MFPop, params: Union[Dict, MFParams], from_pop: MFPop, connection: ConnectionStrategy=Connection.all_to_all()):
        super().__init__(name, pop, params, from_pop, connection)

        defaults = {}
        expectations = {
            IP.TAU_RISE: units.second,
            IP.TAU_D1: units.second,
            IP.TAU_D2: units.second,
            IP.ALPHA: 1
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    def g_dyn(self):
        taus = [
            self.params[IP.ALPHA] * self.params[IP.TAU_D1],
            (1. - self.params[IP.ALPHA]) * self.params[IP.TAU_D2],
            - self.params[IP.ALPHA] * self.params[IP.TAU_D1] * self.params[IP.TAU_RISE] / (self.params[IP.TAU_D1] + self.params[IP.TAU_RISE]),
            - (1. - self.params[IP.ALPHA]) * self.params[IP.TAU_D2] * self.params[IP.TAU_RISE] / (self.params[IP.TAU_D2] + self.params[IP.TAU_RISE])
        ]
        return self.from_pop.n * self.from_pop.rate * sum(taus)

    @lazyproperty
    def brian2(self):
        model = Equations('''
        w1 : 1
        w2 : 1
        w3 : 1
        w4 : 1
        ''')
        on_pre = '''
        {} += w1
        {} += w2
        {} += w3
        {} += w4
        '''.format(self.post_variable_name_1, self.post_variable_name_2, self.post_variable_name_3, self.post_variable_name_4)
        syn = Synapses(
            source=self.from_pop.brian2,
            target=self.pop.brian2,
            method='euler',
            model=model,
            on_pre=on_pre,
            name=self.ref,
        )
        syn.connect(j='i')
        syn.w1[:] = 1
        syn.w2[:] = 1
        syn.w3[:] = 1
        syn.w4[:] = 1
        return syn

    @property
    def post_variable_name_1(self):
        return self.post_variable_name + '_1'

    @property
    def post_variable_name_2(self):
        return self.post_variable_name + '_2'

    @property
    def post_variable_name_3(self):
        return self.post_variable_name + '_3'

    @property
    def post_variable_name_4(self):
        return self.post_variable_name + '_4'

    def brian2_model(self):

        tau_mix1 = (self.params[IP.TAU_RISE] * self.params[IP.TAU_D1]) / (self.params[IP.TAU_RISE] + self.params[IP.TAU_D1])
        tau_mix2 = (self.params[IP.TAU_RISE] * self.params[IP.TAU_D2]) / (self.params[IP.TAU_RISE] + self.params[IP.TAU_D2])

        return Equations(
            '''
            I = g * (v - vrev) * (a * s1 + (1 - a) * s2 - a * s3 - (1 - a) * s4) : amp
            ds1 / dt = - s1 / tau_d1 : 1
            ds2 / dt = - s2 / tau_d2 : 1
            ds3 / dt = - s3 / tau_mix1 : 1
            ds4 / dt = - s4 / tau_mix2 : 1
            ''',
            I=self.current_name,
            g=self.params[IP.GM],
            s1=self.post_variable_name_1,
            s2=self.post_variable_name_2,
            s3=self.post_variable_name_3,
            s4=self.post_variable_name_4,
            vrev=self.params[IP.VREV],
            tau_d1=self.params[IP.TAU_D1],
            tau_d2=self.params[IP.TAU_D2],
            tau_mix1=tau_mix1,
            tau_mix2=tau_mix2,
            a=self.params[IP.ALPHA]
        )

