from brian2 import StateMonitor, defaultclock, Network
from brian2.units import *
import numpy as np

from meanfield.sources.MFLinear3TSSource import MFLinear3TSSource
from meanfield.populations.MFLinearPop import MFLinearPop
from meanfield.populations.MFPoissonSource import MFPoissonSource
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.MFSystem import MFSystem
from meanfield.parameters import NP
from meanfield.parameters import SP

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

class TestMFLinear3TSSource(object):

    def test_simulation_theory(self):

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
        syn = MFLinear3TSSource('syn', pop, {
            SP.GM: 10 * nsiemens,
            SP.VREV: 0 * volt,
            SP.TAU: 0 * ms, # unused
            SP.TAU_RISE: 10 * ms,
            SP.TAU_D1: 20 * ms,
            SP.TAU_D2: 30 * ms,
            SP.ALPHA: 1
        }, poisson)

        system = MFSystem('test')
        system.pops += [pop]
        solver = MFSolverRatesVoltages(system, solver='mse')
        solver.run()
        theory = syn.g_dyn() / syn.from_pop.n

        print([syn.post_variable_name_1, syn.post_variable_name_2, syn.post_variable_name_3])
        print(syn.b2_syn)

        m = StateMonitor(syn.b2_syn, [syn.post_variable_name_1, syn.post_variable_name_2, syn.post_variable_name_3], record=True)
        defaultclock.dt = dt
        net = Network()
        net.add(poisson.brian2)
        net.add(pop.brian2)
        net.add(syn.b2_syn)
        net.add(m)
        net.run(t)

        stable_t = int(t / dt * 0.1)
        simulation_1 = m.__getattr__(syn.post_variable_name_1)[:, stable_t:]
        simulation_mean_1 = np.mean(simulation_1)
        simulation_2 = m.__getattr__(syn.post_variable_name_2)[:, stable_t:]
        simulation_mean_2 = np.mean(simulation_2)
        simulation_3 = m.__getattr__(syn.post_variable_name_3)[:, stable_t:]
        simulation_mean_3 = np.mean(simulation_3)

        #assert np.isclose(theory, simulation_mean, rtol=0.5, atol=0.5)
        print(simulation_mean_1)
        print(simulation_mean_2)
        print(simulation_mean_3)
        print(20 * ms * 10 * Hz)
        # TODO : post_variable = tau * nu


    def test_model_gen(self):
        pass


        # simulate poisson spike train (100 neurons) connect all with same source, record s (should be linear function of input), 10 rates * 10 reps
        # s = tau * nu

        ## s ~= g_dyn (theory)
        # TTS different slope but still linear
