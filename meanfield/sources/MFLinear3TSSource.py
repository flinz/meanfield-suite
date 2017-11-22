from brian2 import Equations, Synapses, units

from meanfield.sources.MFLinearSource import MFLinearSource
from meanfield.utils import lazyproperty
from meanfield.parameters import SP


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

    @lazyproperty
    def b2_syn(self, mode='i != j', method='euler', weight1=1, weight2=1, weight3=1):
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
            method=method,
            model=model,
            on_pre=on_pre
        )
        syn.connect(j='i')
        syn.w1[:] = weight1
        syn.w2[:] = weight2
        syn.w3[:] = weight3
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

    def b2_dyn(self):

        tau_mix1 = (self.params[SP.TAU_RISE] * self.params[SP.TAU_D1]) / (self.params[SP.TAU_RISE] + self.params[SP.TAU_D1])
        tau_mix2 = (self.params[SP.TAU_RISE] * self.params[SP.TAU_D2]) / (self.params[SP.TAU_RISE] + self.params[SP.TAU_D2])

        return Equations(
            '''
            I = g * (v - vrev) * (a * s1 + (1 - a) * s2 - a * s3 - (1 - a) * s4) : amp
            ds1 / dt = - s1 / tau_d1 : 1
            ds2 / dt = - s2 / tau_d2 : 1
            ds3 / dt = - s3 / tau_mix1 : 1
            ds4 / dt = - s4 / tau_mix2 : 1
            ''',
            I=self.current_name,
            g=self.params[SP.GM],
            s1=self.post_variable_name_1,
            s2=self.post_variable_name_2,
            s3=self.post_variable_name_3,
            s4=self.post_variable_name_4,
            vrev=self.params[SP.VREV],
            tau_d1=self.params[SP.TAU_D1],
            tau_d2=self.params[SP.TAU_D2],
            tau_mix1=tau_mix1,
            tau_mix2=tau_mix2,
            a=self.params[SP.ALPHA]
        )

