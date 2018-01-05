import pytest
from brian2 import PopulationRateMonitor, Network, defaultclock, NeuronGroup, Synapses, PoissonInput
from brian2.units import *
import matplotlib.pyplot as plt
import numpy as np

from meanfield.MFSystem import MFSystem
from meanfield.parameters import PP
from meanfield.parameters import IP
from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.inputs.MFStaticInput import MFStaticInput
from meanfield.inputs.MFLinearInput import MFLinearInput
from tests.utils import enable_cpp
from meanfield.utils import brian2_introspect, reset_brian2


class TestRecurrent(object):

    def test_simulation_theory(self):
        enable_cpp()

        t = 1000 * ms
        dt = 0.01 * ms
        defaultclock.dt = dt

        pop = MFLinearPopulation(100, {
            PP.GM: 25. * nS,
            PP.CM: 0.5 * nF,
            PP.VL: -70. * mV,
            PP.VTHR: -50. * mV,
            PP.VRES: -55. * mV,
            PP.TAU_RP: 2. * ms
        })
        pop.rate = 1 * Hz

        noise = MFStaticInput(1000, 1 * Hz, pop, {
            IP.GM: 2 * nS,
            IP.VREV: 0 * volt,
            IP.TAU: 2. * ms,
        })

        rec = MFLinearInput(pop, pop, {
            IP.GM: 0.973 / 100 * nS,
            IP.VREV: -70 * volt,
            IP.TAU: 10. * ms,
        })

        system = MFSystem(pop)

        solver = MFSolverRatesVoltages(system, solver='mse', maxiter=1)
        sol = solver.run()
        theory = sol.state[0]

        rate = PopulationRateMonitor(pop.brian2)

        net = system.collect_brian2_network(rate)
        net.run(t)

        #brian2_introspect(net, globals())

        stable_t = int(t / dt * 0.1)
        isolated = np.array(rate.rate)[stable_t:-stable_t]
        print(isolated.mean())

        if False:
            plt.plot(rate.t / ms, rate.smooth_rate(width=25 * ms) / Hz)
            plt.plot(np.ones(10000) * isolated.mean(), label='simulation mean')
            plt.plot(np.ones(10000) * theory, label='theory mean')
            plt.ylabel('Population rate (Hz) per 100')
            plt.xlabel('Simulation time (ms)')
            plt.title('Poisson noise 5 Hz per 1000')
            plt.legend()
            plt.show()

    @pytest.mark.skip(reason="fitting")
    def test_theory(self):

        n_pop = 25
        n_noise = 100

        def for_rate(rate):
            reset_brian2()

            pop = MFLinearPopulation(n_pop, {
                PP.GM: 25. * nS,
                PP.CM: 0.5 * nF,
                PP.VL: -70. * mV,
                PP.VTHR: -50. * mV,
                PP.VRES: -55. * mV,
                PP.TAU_RP: 2. * ms
            })
            pop.rate = 1 * Hz

            MFStaticInput(n_noise, rate * Hz, pop, {
                IP.GM: 2 * nS,
                IP.VREV: 0 * volt,
                IP.TAU: 2. * ms,
            })

            t = 500 * ms
            dt = 0.1 * ms

            system = MFSystem(pop)
            rate = PopulationRateMonitor(pop.brian2)
            net = system.collect_brian2_network(rate)
            net.run(t)
            stable_t = int(t / dt * 0.1)
            isolated = np.array(rate.rate)[stable_t:-stable_t]
            sim = np.mean(isolated)

            solver = MFSolverRatesVoltages(system, solver='simplex', maxiter=1)
            sol = solver.run()
            return [sol.state[0], sim]


        rates = np.linspace(1, 100, 20)
        dom = np.array([for_rate(r) for r in rates])
        print(dom)
        plt.plot(rates, dom[:, 0], label='theory')
        plt.plot(rates, dom[:, 1], label='simulation')
        plt.xlabel(f'Noise rate (Hz) on {n_pop}')
        plt.ylabel(f'Population rate (Hz) from {n_noise}')
        plt.legend()
        plt.show()

        print(dom[:, 0] - dom[:, 1])
        plt.plot(rates, dom[:, 0] - dom[:, 1])
        plt.show()
