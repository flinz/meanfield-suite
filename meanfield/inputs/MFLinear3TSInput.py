from types import MappingProxyType
from typing import Union, Dict

from brian2 import Equations, Synapses, units

from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.utils import lazyproperty
from meanfield.parameters import IP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.parameters import Connection
from meanfield.parameters.Connection import ConnectionStrategy


class MFLinear3TSInput(MFLinearInput):

    arguments = MappingProxyType({
        IP.TAU_RISE: units.second,
        IP.TAU_D1: units.second,
        IP.TAU_D2: units.second,
        IP.ALPHA: 1
    })

    defaults = MappingProxyType({})

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParams], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    def g_dyn(self):
        taus = [
            self[IP.ALPHA] * self[IP.TAU_D1],
            (1. - self[IP.ALPHA]) * self[IP.TAU_D2],
            - self[IP.ALPHA] * self[IP.TAU_D1] * self[IP.TAU_RISE] / (self[IP.TAU_D1] + self[IP.TAU_RISE]),
            - (1. - self[IP.ALPHA]) * self[IP.TAU_D2] * self[IP.TAU_RISE] / (self[IP.TAU_D2] + self[IP.TAU_RISE])
        ]
        return self.origin.n * self.origin.rate * sum(taus)

    # Simulation

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
            source=self.origin.brian2,
            target=self.target.brian2,
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

        tau_mix1 = (self[IP.TAU_RISE] * self[IP.TAU_D1]) / (self[IP.TAU_RISE] + self[IP.TAU_D1])
        tau_mix2 = (self[IP.TAU_RISE] * self[IP.TAU_D2]) / (self[IP.TAU_RISE] + self[IP.TAU_D2])

        return Equations(
            '''
            I = g * (v - vrev) * (a * s1 + (1 - a) * s2 - a * s3 - (1 - a) * s4) : amp
            ds1 / dt = - s1 / tau_d1 : 1
            ds2 / dt = - s2 / tau_d2 : 1
            ds3 / dt = - s3 / tau_mix1 : 1
            ds4 / dt = - s4 / tau_mix2 : 1
            ''',
            I=self.current_name,
            g=self[IP.GM],
            s1=self.post_variable_name_1,
            s2=self.post_variable_name_2,
            s3=self.post_variable_name_3,
            s4=self.post_variable_name_4,
            vrev=self[IP.VREV],
            tau_d1=self[IP.TAU_D1],
            tau_d2=self[IP.TAU_D2],
            tau_mix1=tau_mix1,
            tau_mix2=tau_mix2,
            a=self[IP.ALPHA]
        )

