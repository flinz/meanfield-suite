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

class MFSourceTests(unittest.TestCase):

    def testModelGen(self):
        pop = MFLinearPop("test", 1, params)
        source = MFSource('test', pop)

        current, eqs = source.brian2_model(1)

        assert_equations(
            eqs,
            '''
            I1 = (0. * siemens) * (v - (0. * volt)) * s1 : A
            ds1 / dt = -s1 / (10. * msecond) : 1
            '''
        )

    def testModelUpdate(self):
        pop = MFLinearPop("test", 1, params)
        source = MFSource('test', pop)

        assert_equations(
            source.brian_link(),
            '''
            I1 = (0. * siemens) * (v - (0. * volt)) * s1 : A
            ds1 / dt = -s1 / (10. * msecond) : 1
            ''',
            True
        )

