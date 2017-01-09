import unittest

from brian2 import *

from MFParams import MFParams
from MFPop import MFLinearPop
from MFSource import MFSource
from params import NP
from tests.utils import assert_equations

params = MFParams({
    NP.GAMMA: 0.280112,
    NP.BETA: 0.062,
    NP.GM: 25. * nS,
    NP.CM: 0.5 * nF,  # * 1e3,
    NP.VL: -70. * mV,
    NP.VTHR: -50. * mV,
    NP.VRES: -55. * mV,
    NP.TAU_RP: 2. * ms
})

class MFPopTests(unittest.TestCase):

    def testModelGenWithoutSource(self):
        pop = MFLinearPop("test", 1, params)

        assert_equations(
            pop.brian_v(),
            '''
            I = 0 : A
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad)  : V (unless refractory)
            '''
        )

    def testModelGenWithOneSource(self):
        pop = MFLinearPop("test", 1, params)
        source = MFSource('test', pop)

        assert_equations(
            pop.brian_v(),
            '''
            I0 = (0. * siemens) * (v - (0. * volt)) * s0 : A
            I = I0 : A
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad) : V (unless refractory)
            ds0/dt = - s0 / (10. * msecond) : 1
            '''
        )

    def testModelGenWithTwoSources(self):
        pop = MFLinearPop("test", 1, params)
        source1 = MFSource('test1', pop)
        source2 = MFSource('test2', pop)

        assert_equations(
            pop.brian_v(),
            '''
            I0 = (0. * siemens) * (v - (0. * volt)) * s0  : A
            I1 = (0. * siemens) * (v - (0. * volt)) * s1  : A
            I= I0 + I1 : A
            dv/dt = (- (25. * nsiemens) * (v - (-70. * mvolt)) - I) / (0.5 * nfarad)  : V (unless refractory)
            ds0/dt = - s0 / (10. * msecond)  : 1
            ds1/dt = - s1 / (10. * msecond)  : 1
            '''
        )

