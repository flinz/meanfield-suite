import pytest
from brian2 import PopulationRateMonitor, Network, defaultclock
from brian2.units import *
import matplotlib.pyplot as plt
import numpy as np

from meanfield.MFSystem import MFSystem
from meanfield.parameters import NP
from meanfield.parameters import SP
from meanfield.populations.MFLinearPop import MFLinearPop
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.sources.MFStaticSource import MFStaticSource
from tests.utils import enable_cpp


class TestNoise(object):

    def test_simulation_theory(self):
        enable_cpp()

        t = 10000 * ms
        dt = 0.01 * ms
        defaultclock.dt = dt

        pop = MFLinearPop("pop", 100, {
            NP.GM: 25. * nS,
            NP.CM: 0.5 * nF,
            NP.VL: -70. * mV,
            NP.VTHR: -50. * mV,
            NP.VRES: -55. * mV,
            NP.TAU_RP: 2. * ms
        })
        pop.rate = 10 * Hz

        noise = MFStaticSource("noise", pop, 1000, 3 * Hz, {
            SP.GM: 2 * nS,
            SP.VREV: 0 * volt,
            SP.TAU: 2. * ms,
        })
        pop.add_noise(noise)

        system = MFSystem("pop noise")
        system.pops += [pop]

        solver = MFSolverRatesVoltages(system, solver='mse', maxiter=1)
        sol = solver.run()

        if False:
            theory = sol.state[0]
            rate = PopulationRateMonitor(pop.brian2)

            net = Network()
            net.add(pop.brian2)
            net.add(noise.brian2)
            net.add(rate)
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
            pop = MFLinearPop("pop", 100, {
                NP.GM: 25. * nS,
                NP.CM: 0.5 * nF,
                NP.VL: -70. * mV,
                NP.VTHR: -50. * mV,
                NP.VRES: -55. * mV,
                NP.TAU_RP: 2. * ms
            })
            pop.rate = 10 * Hz

            noise = MFStaticSource("noise", pop, 1000, rate * Hz, {
                SP.GM: 2 * nS,
                SP.VREV: 0 * volt,
                SP.TAU: 2. * ms,
            })
            pop.add_noise(noise)

            system = MFSystem("pop noise")
            system.pops += [pop]

            solver = MFSolverRatesVoltages(system, solver='mse', maxiter=1)
            print(solver.state)
            sol = solver.run()
            return sol.state[0]

        rates = np.linspace(1, 10, 20)
        dom = np.array([for_rate(r) for r in rates])
        plt.plot(rates, dom)
        plt.xlabel('Noise rate (Hz) per 1000')
        plt.ylabel('Population rate (Hz) per 100')
        plt.show()


