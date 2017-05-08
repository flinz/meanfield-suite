from brian2 import Equations

from MFLinearSource import MFLinearSource
from Utils import lazy
from params import SP


class MFTLinearTTSSource(MFLinearSource):
    def __init__(self, name, pop, params, from_pop):
        super().__init__(name, pop, params, from_pop)

        defaults = {}
        expectations = {
            SP.GM_1: 1,
            SP.GM_2: 1,
            SP.GM_3: 1,
            SP.TAU_1: 1,
            SP.TAU_2: 1,
            SP.TAU_3: 1,
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    @lazy
    def b2_syn(self):
        # TODO
        return None

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
            I = g1 * (v - vrev) * s1 + g2 * (v - vrev) * s2 + g3 * (v - vrev) * s3 : amp
            ds1 / dt = - s1 / tau1 : 1
            ds2 / dt = - s2 / tau2 : 1
            ds3 / dt = - s3 / tau3 : 1
            ''',
            I=self.current_name,
            g=self.params[SP.GM],
            s1=self.post_variable_name_1,
            s2=self.post_variable_name_2,
            s3=self.post_variable_name_3,
            vrev=self.params[SP.VREV],
            tau_1=self.params[SP.TAU_1],
            tau_2=self.params[SP.TAU_2],
            tau_3=self.params[SP.TAU_3],
        )

