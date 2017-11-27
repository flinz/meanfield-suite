from brian2.units import *

from meanfield.populations.MFLinearPop import MFLinearPop
from meanfield.sources.MFSource import MFSource
from meanfield.parameters import NP
from meanfield.parameters import SP
from meanfield.parameters import Connection
from tests.utils import assert_equations

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

class TestMFPop(object):

    def test_model_gen_without_source(self):
        pop = MFLinearPop("test", 1, params_pop)

        assert_equations(
            pop.brian2_model(),
            '''
            I = 0 : amp
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad) : volt (unless refractory)
            '''
        )

    def test_model_gen_with_one_source(self):
        pop = MFLinearPop("test", 1, params_pop)
        _ = MFSource('test', pop, params_source, Connection.all_to_all())

        assert_equations(
            pop.brian2_model(),
            '''
            I_test = (0. * siemens) * (v - (0. * volt)) * s_test : amp
            I = I_test : amp
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad) : volt (unless refractory)
            ds_test/dt = - s_test / (10. * msecond) : 1
            '''
        )

    def test_model_gen_with_two_sources(self):
        pop = MFLinearPop("test", 1, params_pop)
        MFSource('test1', pop, params_source, Connection.all_to_all())
        MFSource('test2', pop, params_source, Connection.all_to_all())

        assert_equations(
            pop.brian2_model(),
            '''
            I_test1 = (0. * siemens) * (v - (0. * volt)) * s_test1 : amp
            I_test2 = (0. * siemens) * (v - (0. * volt)) * s_test2 : amp
            I= I_test1 + I_test2 : amp
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad) : volt (unless refractory)
            ds_test1/dt = - s_test1 / (10. * msecond)  : 1
            ds_test2/dt = - s_test2 / (10. * msecond)  : 1
            '''
        )

