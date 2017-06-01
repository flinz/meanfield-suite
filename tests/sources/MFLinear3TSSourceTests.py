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

class MFStaticSourceTests(unittest.TestCase):

    def testModelGen(self):
        pass


        # simulate poisson spike train (100 neurons) connect all with same source, record s (should be linear function of input), 10 rates * 10 reps
        # s = tau * nu

        ## s ~= g_dyn (theory)
        # TTS different slope but still linear
