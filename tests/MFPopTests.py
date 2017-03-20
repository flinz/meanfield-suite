import unittest

from brian2 import *

from MFPop import MFLinearPop
from MFSource import MFSource
from params import NP
from params import SP
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
    SP.VE: 0 * volt,
    SP.TAU: 10 * ms,
}

class MFPopTests(unittest.TestCase):

    def testModelGenWithoutSource(self):
        pop = MFLinearPop("test", 1, params_pop)

        assert_equations(
            pop.brian2_model(),
            '''
            I = 0 : A
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad)  : V (unless refractory)
            '''
        )

    def testModelGenWithOneSource(self):
        pop = MFLinearPop("test", 1, params_pop)
        _ = MFSource('test', pop, params_source)

        assert_equations(
            pop.brian2_model(),
            '''
            I_test = (0. * siemens) * (v - (0. * volt)) * s_test : A
            I = I_test : A
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad) : V (unless refractory)
            ds_test/dt = - s_test / (10. * msecond) : 1
            '''
        )

    def testModelGenWithTwoSources(self):
        pop = MFLinearPop("test", 1, params_pop)
        _ = MFSource('test1', pop, params_source)
        _ = MFSource('test2', pop, params_source)

        assert_equations(
            pop.brian2_model(),
            '''
            I_test1 = (0. * siemens) * (v - (0. * volt)) * s_test1 : A
            I_test2 = (0. * siemens) * (v - (0. * volt)) * s_test2 : A
            I= I_test1 + I_test2 : A
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad) : V (unless refractory)
            ds_test1/dt = - s_test1 / (10. * msecond)  : 1
            ds_test2/dt = - s_test2 / (10. * msecond)  : 1
            '''
        )

