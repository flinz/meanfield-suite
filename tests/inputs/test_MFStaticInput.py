from brian2.units import *

from meanfield.parameters import IP
from meanfield.parameters import PP

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
    IP.VE: 0 * volt,
    IP.TAU: 10 * ms,
}

class TestMFStaticInput(object):

    def test_model_gen(self):
        pass

