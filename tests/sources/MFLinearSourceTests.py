import unittest

from brian2 import *
from brian2 import units

from MFLinearPop import MFLinearPop
from MFLinearSource import MFLinearSource
from MFSolver import MFSolverRatesVoltages
from MFSource import MFSource
from MFSystem import MFSystem
from Utils import lazy
from params import NP
from params import SP

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


    def testSimulation(self):

        class MFPoissonSource:
            @check_units(rate=units.Hz)
            def __init__(self, name, n, rate):
                self.n = n
                self.rate = rate
            @lazy
            def brian2(self):
                return PoissonGroup(self.n, self.rate)


        i = MFPoissonSource('poisson', 10, 10 * Hz)

        p_params = {
            NP.GM: 10 * nsiemens,
            NP.VL: 0 * mV,
            NP.CM: 5 * nfarad,
            NP.VTHR: 0 * mV,
            NP.VRES: 0 * mV,
            NP.TAU_RP: 15 * ms
        }
        p = MFLinearPop('p', 10, p_params)

        s_params = {
            SP.GM: 10 * nsiemens,
            SP.VREV: 0 * volt,
            SP.TAU: 20 * ms,
        }
        s = MFLinearSource('s', p, s_params, i)

        system = MFSystem('test')
        system.pops += [p]
        solver = MFSolverRatesVoltages(system, solver='mse')
        solver.run()

        net = Network()
        net.add(i.brian2)
        net.add(p.brian2)
        net.add(s.b2_syn)

        m = StateMonitor(s.b2_syn, [s.post_variable_name], record=True)

        net.add(m)

        net.run(3000 * ms)

        print(m.t)
        print(m.__getattr__('s_s')[0])


