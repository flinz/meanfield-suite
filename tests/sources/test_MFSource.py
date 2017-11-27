from brian2.units import *

from meanfield.populations.MFLinearPop import MFLinearPop
from meanfield.sources.MFSource import MFSource
from meanfield.parameters import NP
from meanfield.parameters import SP
from meanfield.parameters import Connection
from tests.utils import assert_equations, enable_cpp

params_pop = {
    NP.GAMMA: 0.280112,
    NP.BETA: 0.062,
    NP.GM: 25. * nS,
    NP.CM: 0.5 * nF,  # * 1e3,
    NP.VL: -70. * mV,
    NP.VTHR: -50. * mV,
    NP.VRES: -55. * mV,
    NP.TAU_RP: 2. * ms
}

params_source = {
    SP.GM: 0 * siemens,
    SP.VREV: 0 * volt,
    SP.TAU: 10 * ms,
}

enable_cpp()

class TestMFSource(object):

    def test_model_gen(self):
        pop = MFLinearPop("test", 1, params_pop)
        source = MFSource('test', pop, params_source, Connection.all_to_all())

        assert_equations(
            source.brian2_model(),
            '''
            I_test = (0. * siemens) * (v - (0. * volt)) * s_test : amp
            ds_test / dt = -s_test / (10. * msecond) : 1
            '''
        )

        # def testModelUpdate(self):
        #     pop = MFLinearPop("test", 1, parameters)
        #     source = MFSource('test', pop)
        #
        #     assert_equations(
        #         source.brian2(),
        #         '''
        #         I1 = (0. * siemens) * (v - (0. * volt)) * s1 : A
        #         ds1 / dt = -s1 / (10. * msecond) : 1
        #         ''',
        #         True
        #     )

