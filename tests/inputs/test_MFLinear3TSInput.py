from brian2 import StateMonitor, defaultclock, Network
from brian2.units import *
import numpy as np

from meanfield.inputs.MFLinear3TSInput import MFLinear3TSInput
from meanfield.populations.MFLinearPop import MFLinearPop
from meanfield.populations.MFPoissonPop import MFPoissonPop
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.MFSystem import MFSystem
from meanfield.parameters import PP
from meanfield.parameters import IP
from tests.utils import enable_cpp

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

class TestMFLinear3TSInput(object):

    def test_simulation_theory(self):

        enable_cpp()

        t = 3000 * ms
        dt = 0.01 * ms
        n = 100

        alpha = 0.5

        poisson = MFPoissonPop('poisson', n, n * 10 * Hz, {
            PP.GM: 1 * nsiemens,
            PP.VRES: 0 * mV,
            PP.TAU_RP: 2 * ms
        })
        pop = MFLinearPop('pop', n, {
            PP.GM: 10 * nsiemens,
            PP.VL: 0 * mV,
            PP.CM: 5 * nfarad,
            PP.VTHR: 0 * mV,
            PP.VRES: 0 * mV,
            PP.TAU_RP: 15 * ms
        })
        syn = MFLinear3TSInput('syn', pop, {
            IP.GM: 10 * nsiemens,
            IP.VREV: 0 * mvolt,
            IP.TAU: 0 * ms, # unused
            IP.TAU_RISE: 2 * ms,
            IP.TAU_D1: 20 * ms,
            IP.TAU_D2: 30 * ms,
            IP.ALPHA: alpha
        }, poisson)

        system = MFSystem('test')
        system.pops += [pop]
        solver = MFSolverRatesVoltages(system, solver='mse')
        solver.run()
        theory = syn.g_dyn() / syn.from_pop.n

        print([syn.post_variable_name_1, syn.post_variable_name_2, syn.post_variable_name_3, syn.post_variable_name_4])
        print(syn.brian2)

        m1 = StateMonitor(syn.brian2, syn.post_variable_name_1, record=range(100))
        m2 = StateMonitor(syn.brian2, syn.post_variable_name_2, record=range(100))
        m3 = StateMonitor(syn.brian2, syn.post_variable_name_3, record=range(100))
        m4 = StateMonitor(syn.brian2, syn.post_variable_name_4, record=range(100))
        defaultclock.dt = dt
        net = Network()
        net.add(poisson.brian2)
        net.add(pop.brian2)
        net.add(syn.brian2)
        net.add(m1)
        net.add(m2)
        net.add(m3)
        net.add(m4)

        # FIXME
        print(Network.__instances__())
        print(net.objects)
        net.run(t)

        stable_t = int(t / dt * 0.1)
        simulation_1 = m1.__getattr__(syn.post_variable_name_1)[:, stable_t:]
        simulation_mean_1 = np.mean(simulation_1)
        simulation_2 = m2.__getattr__(syn.post_variable_name_2)[:, stable_t:]
        simulation_mean_2 = np.mean(simulation_2)
        simulation_3 = m3.__getattr__(syn.post_variable_name_3)[:, stable_t:]
        simulation_mean_3 = np.mean(simulation_3)
        simulation_4 = m4.__getattr__(syn.post_variable_name_4)[:, stable_t:]
        simulation_mean_4 = np.mean(simulation_4)

        simulation_mean = alpha * simulation_mean_1 + (1 - alpha) * simulation_mean_2 + - alpha * simulation_mean_3 - (1 - alpha) * simulation_mean_4

        # TODO : post_variable = tau * nu
        assert np.isclose(theory, simulation_mean, rtol=0.5, atol=0.5)



    def test_model_gen(self):
        pass


        # simulate poisson spike train (100 neurons) connect all with same source, record s (should be linear function of input), 10 rates * 10 reps
        # s = tau * nu

        ## s ~= g_dyn (theory)
        # TTS different slope but still linear
