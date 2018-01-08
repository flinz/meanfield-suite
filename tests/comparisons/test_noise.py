import pytest
from brian2 import PopulationRateMonitor, defaultclock
from brian2.units import *
import matplotlib.pyplot as plt
import numpy as np

from core.MFSystem import MFSystem
from meanfield.parameters import PP
from meanfield.parameters import IP
from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.solvers.MFSolver import MFSolver
from meanfield.inputs.MFStaticInput import MFStaticInput
from tests.utils import enable_cpp


class TestNoise(object):

    def test_simulation_theory(self):
        enable_cpp()

        t = 10000 * ms
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
        pop.rate = 10 * Hz

        noise = MFStaticInput(1000, 3 * Hz, pop, {
            IP.GM: 2 * nS,
            IP.VREV: 0 * volt,
            IP.TAU: 2. * ms,
        })

        system = MFSystem(pop)

        solver = MFSolver.rates_voltages(system, solver='mse', maxiter=1)
        sol = solver.run()

        if False:
            theory = sol.state[0]
            rate = PopulationRateMonitor(pop.brian2)

            net = system.collect_brian2_network(rate)
            net.run(t)

            stable_t = int(t / dt * 0.1)
            isolated = np.array(rate.rate)[stable_t:-stable_t]
            print(isolated.mean())

            plt.plot(rate.t / ms, rate.smooth_rate(width=50 * ms) / Hz)
            plt.plot(np.ones(10000) * isolated.mean(), label='simulation mean')
            plt.plot(np.ones(10000) * theory, label='theory mean')
            plt.ylabel('Population rate (Hz) per 100')
            plt.xlabel('Simulation time (ms)')
            plt.title('Poisson noise 3 Hz per 1000')
            plt.legend()
            plt.show()

    @pytest.mark.skip(reason="plotting")
    def test_theory(self):
        enable_cpp()

        def for_rate(rate):
            pop = MFLinearPopulation(100, {
                PP.GM: 25. * nS,
                PP.CM: 0.5 * nF,
                PP.VL: -70. * mV,
                PP.VTHR: -50. * mV,
                PP.VRES: -55. * mV,
                PP.TAU_RP: 2. * ms
            })
            pop.rate = 10 * Hz

            noise = MFStaticInput(1000, rate * Hz, pop, {
                IP.GM: 2 * nS,
                IP.VREV: 0 * volt,
                IP.TAU: 2. * ms,
            })
            pop.add_noise(noise)

            system = MFSystem(pop)

            solver = MFSolver.rates_voltages(system, solver='mse', maxiter=1)
            print(solver.state)
            sol = solver.run()
            return sol.state[0]

        rates = np.linspace(1, 10, 20)
        dom = np.array([for_rate(r) for r in rates])
        plt.plot(rates, dom)
        plt.xlabel('Noise rate (Hz) per 1000')
        plt.ylabel('Population rate (Hz) per 100')
        plt.show()


