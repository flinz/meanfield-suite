from brian2.units import *

from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.inputs.MFInput import MFInput
from meanfield.parameters import PP
from meanfield.parameters import IP
from meanfield.parameters import Connection
from tests.utils import assert_equations, enable_cpp

params_pop = {
    PP.GAMMA: 0.280112,
    PP.BETA: 0.062,
    PP.GM: 25. * nS,
    PP.CM: 0.5 * nF,  # * 1e3,
    PP.VL: -70. * mV,
    PP.VTHR: -50. * mV,
    PP.VRES: -55. * mV,
    PP.TAU_RP: 2. * ms
}

params_source = {
    IP.GM: 0 * siemens,
    IP.VREV: 0 * volt,
    IP.TAU: 10 * ms,
}


class TestMFInput(object):

    def test_model_gen(self):
        enable_cpp()
        pop = MFLinearPopulation(1, params_pop)
        source = MFInput('test', pop, params_source, Connection.all_to_all())

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

