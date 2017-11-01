import unittest

from brian2 import *

from sources.MFLinearNMDASource import MFLinearNMDASource
from populations.MFLinearPop import MFLinearPop
from populations.MFPoissonSource import MFPoissonSource
from solvers.MFSolver import MFSolverRatesVoltages
from MFSystem import MFSystem
from parameters import NP
from parameters import SP

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

class MFLinearNMDASourceTests(unittest.TestCase):

    def testSimulationVsTheory(self):
        set_device('cpp_standalone')
        #prefs.codegen.target = 'cython'

        t = 3000 * ms
        dt = 0.01 * ms
        n = 100

        poisson = MFPoissonSource('poisson', n, n * 10 * Hz)
        pop = MFLinearPop('pop', n, {
            NP.GM: 10 * nsiemens,
            NP.VL: 0 * mV,
            NP.CM: 5 * nfarad,
            NP.VTHR: 0 * mV,
            NP.VRES: 0 * mV,
            NP.TAU_RP: 15 * ms
        })
        syn = MFLinearNMDASource('syn', pop, {
            SP.GM: 10 * nsiemens,
            SP.VREV: 0 * volt,
            SP.TAU: 20 * ms,
            SP.TAU_NMDA: 30 * ms,
            SP.ALPHA: 1, # TODO git aalpha ?
            SP.BETA: 1,
            SP.GAMMA: 1,
        }, poisson)

        system = MFSystem('test')
        system.pops += [pop]
        solver = MFSolverRatesVoltages(system, solver='mse')
        solver.run()
        theory = syn.g_dyn() / syn.from_pop.n

        m = StateMonitor(syn.b2_syn, syn.post_variable_name, record=range(100))
        defaultclock.dt = dt
        net = Network()
        net.add(poisson.brian2)
        net.add(pop.brian2)
        net.add(syn.b2_syn)
        net.add(m)
        net.run(t)

        stable_t = int(t / dt * 0.1)
        simulation = m.__getattr__(syn.post_variable_name)[:, stable_t:]
        simulation_mean = np.mean(simulation)
        print(simulation)

        assert np.isclose(theory, simulation_mean, rtol=0.5, atol=0.5)
        print(simulation_mean)
        print(20 * ms * 10 * Hz)
        # TODO : post_variable = tau * nu

