from brian2 import StateMonitor, Network, defaultclock
from brian2.units import *
import numpy as np

from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.populations.MFPoissonPopulation import MFPoissonPopulation
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.MFSystem import MFSystem
from meanfield.parameters import PP
from meanfield.parameters import IP
from meanfield.parameters import Connection
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


class TestMFLinearInput(object):

    def test_simulation_theory(self):
        enable_cpp()

        t = 3000 * ms
        dt = 0.01 * ms
        n = 100

        poisson = MFPoissonPopulation(n, 10 * Hz, {
            PP.GM: 0 * nsiemens,
            PP.VRES: 0 * mV,
            PP.TAU_RP: 0 * ms
        })
        pop = MFLinearPopulation(n, {
            PP.GM: 10 * nsiemens,
            PP.VL: 0 * mV,
            PP.CM: 5 * nfarad,
            PP.VTHR: 0 * mV,
            PP.VRES: 0 * mV,
            PP.TAU_RP: 15 * ms
        })
        syn = MFLinearInput(poisson, pop, {
            IP.GM: 10 * nsiemens,
            IP.VREV: 0 * mV,
            IP.TAU: 20 * ms,
        }, connection=Connection.one_to_one())

        system = MFSystem('test')
        system.pops += [pop]
        solver = MFSolverRatesVoltages(system, solver='mse')
        solver.run()
        theory = syn.g_dyn() / syn.origin.n

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


