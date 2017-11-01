from brian2 import Equations, Synapses, units

from MFLinearSource import MFLinearSource
from utils import lazy
from params import SP


class MFLinear3TSSource(MFLinearSource):
    def __init__(self, name, pop, params, from_pop):
        super().__init__(name, pop, params, from_pop)

        defaults = {}
        expectations = {
            SP.TAU_RISE: units.second,
            SP.TAU_D1: units.second,
            SP.TAU_D2: units.second,
            SP.ALPHA: 1
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    def g_dyn(self):
        taus = [
            self.params[SP.ALPHA] * self.params[SP.TAU_D1],
            (1. - self.params[SP.ALPHA]) * self.params[SP.TAU_D2],
            - self.params[SP.ALPHA] * self.params[SP.TAU_D1] * self.params[SP.TAU_RISE] / (self.params[SP.TAU_D1] + self.params[SP.TAU_RISE]),
            - (1. - self.params[SP.ALPHA]) * self.params[SP.TAU_D2] * self.params[SP.TAU_RISE] / (self.params[SP.TAU_D2] + self.params[SP.TAU_RISE])
        ]
        return self.from_pop.n * self.from_pop.rate * sum(taus)

    @lazy
    def b2_syn(self, mode='i != j', method='euler', weight1=1, weight2=1, weight3=3):
        model = Equations('''
        w1 : 1
        w2 : 1
        w3 : 1
        ''')
        on_pre = '''
        {} += w1
        {} += w2
        {} += w3
        '''.format(self.post_variable_name_1, self.post_variable_name_2, self.post_variable_name_3)
        print(self.from_pop)
        syn = Synapses(
            source=self.from_pop.brian2,
            target=self.pop.brian2,
            method=method,
            model=model,
            on_pre=on_pre
        )
        syn.connect(mode)
        syn.w1[:] = weight1
        syn.w2[:] = weight2
        syn.w3[:] = weight3
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

    def b2_dyn(self):
        return Equations(
            '''
            I = g * (v - vrev) * (1 - s1) * (a * s2 + (1 - a) * s3) : amp
            ds1 / dt = - s1 / tau_1 : 1
            ds2 / dt = - s2 / tau_2 : 1
            ds3 / dt = - s3 / tau_3 : 1
            ''',
            I=self.current_name,
            g=self.params[SP.GM],
            s1=self.post_variable_name_1,
            s2=self.post_variable_name_2,
            s3=self.post_variable_name_3,
            vrev=self.params[SP.VREV],
            tau_1=self.params[SP.TAU_RISE],
            tau_2=self.params[SP.TAU_D1],
            tau_3=self.params[SP.TAU_D2],
            a=self.params[SP.ALPHA]
        )

