import unittest

from brian2 import units

from populations.MFLinearPop import MFLinearPop
from sources.MFSource import MFSource
from parameters import NP
from parameters import SP
from tests.utils import assert_equations

params_pop = {
    NP.GAMMA: 0.280112,
    NP.BETA: 0.062,
    NP.GM: 25. * units.nS,
    NP.CM: 0.5 * units.nF,  # * 1e3,
    NP.VL: -70. * units.mV,
    NP.VTHR: -50. * units.mV,
    NP.VRES: -55. * units.mV,
    NP.TAU_RP: 2. * units.ms
}

params_source = {
    SP.GM: 0 * units.siemens,
    SP.VREV: 0 * units.volt,
    SP.TAU: 10 * units.ms,
}

class MFPopTests(unittest.TestCase):

    def testModelGenWithoutSource(self):
        pop = MFLinearPop("test", 1, params_pop)

        assert_equations(
            pop.brian2_model(),
            '''
            I = 0 : amp
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad) : volt (unless refractory)
            '''
        )

    def testModelGenWithOneSource(self):
        pop = MFLinearPop("test", 1, params_pop)
        _ = MFSource('test', pop, params_source)

        assert_equations(
            pop.brian2_model(),
            '''
            I_test = (0. * siemens) * (v - (0. * volt)) * s_test : amp
            I = I_test : amp
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad) : volt (unless refractory)
            ds_test/dt = - s_test / (10. * msecond) : 1
            '''
        )

    def testModelGenWithTwoSources(self):
        pop = MFLinearPop("test", 1, params_pop)
        MFSource('test1', pop, params_source)
        MFSource('test2', pop, params_source)

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

