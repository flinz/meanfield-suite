from brian2 import StateMonitor, Network, defaultclock
from brian2.units import *
import numpy as np

from meanfield.populations.MFLinearPop import MFLinearPop
from meanfield.sources.MFLinearSource import MFLinearSource
from meanfield.populations.MFPoissonPop import MFPoissonPop
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.MFSystem import MFSystem
from meanfield.parameters import NP
from meanfield.parameters import SP
from meanfield.parameters import Connection
from tests.utils import enable_cpp

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


class TestMFLinearSource(object):

    def test_simulation_theory(self):
        enable_cpp()

        t = 3000 * ms
        dt = 0.01 * ms
        n = 100

        poisson = MFPoissonPop('poisson', n, 10 * Hz, {
            NP.GM: 0 * nsiemens,
            NP.VRES: 0 * mV,
            NP.TAU_RP: 0 * ms
        })
        pop = MFLinearPop('pop', n, {
            NP.GM: 10 * nsiemens,
            NP.VL: 0 * mV,
            NP.CM: 5 * nfarad,
            NP.VTHR: 0 * mV,
            NP.VRES: 0 * mV,
            NP.TAU_RP: 15 * ms
        })
        syn = MFLinearSource('syn', pop, {
            SP.GM: 10 * nsiemens,
            SP.VREV: 0 * mV,
            SP.TAU: 20 * ms,
        }, poisson, Connection.one_to_one())

        system = MFSystem('test')
        system.pops += [pop]
        solver = MFSolverRatesVoltages(system, solver='mse')
        solver.run()
        theory = syn.g_dyn() / syn.from_pop.n

        m = StateMonitor(syn.brian2, syn.post_variable_name, record=range(100))
        defaultclock.dt = dt
        net = Network()
        net.add(poisson.brian2)
        net.add(pop.brian2)
        net.add(syn.brian2)
        net.add(m)
        net.run(t)

        stable_t = int(t / dt * 0.1)
        simulation = m.__getattr__(syn.post_variable_name)[:, stable_t:]
        simulation_mean = np.mean(simulation)

        assert np.isclose(theory, simulation_mean, rtol=0.5, atol=0.5)
        print(simulation_mean)
        print(20 * ms * 10 * Hz)
        # TODO : post_variable = tau * nu

        # FIXME test for all


