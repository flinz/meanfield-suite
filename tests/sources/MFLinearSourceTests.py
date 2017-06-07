import unittest

from brian2 import *
from brian2 import units

from MFLinearPop import MFLinearPop
from MFLinearSource import MFLinearSource
from MFPoissonSource import MFPoissonSource
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

        #noise = MFStaticSource("E_noise1", pop_e1, {
        #    SP.GM: g_AMPA_ext_E,
        #    SP.VE: 0. * mV,
        #    SP.TAU: tau_AMPA
        #}, rate, C_ext)

        t = 3000 * ms
        dt = 0.01 * ms
        n = 100

        i = MFPoissonSource('poisson', n, n * 10 * Hz)

        p_params = {
            NP.GM: 10 * nsiemens,
            NP.VL: 0 * mV,
            NP.CM: 5 * nfarad,
            NP.VTHR: 0 * mV,
            NP.VRES: 0 * mV,
            NP.TAU_RP: 15 * ms
        }
        p = MFLinearPop('p', n, p_params)

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
        theory = s.g_dyn() / s.from_pop.n

        m = StateMonitor(s.b2_syn, s.post_variable_name, record=True)

        defaultclock.dt = dt
        net = Network()
        net.add(i.brian2)
        net.add(p.brian2)
        net.add(s.b2_syn)
        net.add(m)
        net.run(t)

        stable_t = int(t / dt * 0.1)
        simulation = m.__getattr__(s.post_variable_name)[:, stable_t:]
        simulation_mean = np.mean(simulation)

        assert np.isclose(theory, simulation_mean, rtol=0.5, atol=0.5)
        print(simulation_mean)
        print(20 * ms * 10 * Hz)
        # TODO : post_variable = tau * nu


