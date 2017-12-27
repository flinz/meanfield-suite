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
from meanfield.utils import brian2_introspect


class TestRecurrent(object):

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

        noise = MFStaticInput(1000, 5 * Hz, pop, {
            IP.GM: 2 * nS,
            IP.VREV: 0 * volt,
            IP.TAU: 2. * ms,
        })
        pop.add_noise(noise)

        rec = MFLinearInput(pop, pop, {
            IP.GM: 0.973 / 100 * nS,
            IP.VREV: -70 * volt,
            IP.TAU: 10. * ms,
        })

        system = MFSystem("pop noise rec")
        system.pops += [pop]

        solver = MFSolverRatesVoltages(system, solver='mse', maxiter=1)
        sol = solver.run()
        theory = sol.state[0]

        rate = PopulationRateMonitor(pop.brian2)

        net = Network()
        net.add(pop.brian2)
        net.add(noise.brian2)
        net.add(rec.brian2)
        net.add(rate)



        #brian2_introspect(net, globals())


        if False:
            net.run(t)
            stable_t = int(t / dt * 0.1)
            isolated = np.array(rate.rate)[stable_t:-stable_t]
            print(isolated.mean())

            plt.plot(rate.t / ms, rate.smooth_rate(width=25 * ms) / Hz)
            plt.plot(np.ones(10000) * isolated.mean(), label='simulation mean')
            plt.plot(np.ones(10000) * theory, label='theory mean')
            plt.ylabel('Population rate (Hz) per 100')
            plt.xlabel('Simulation time (ms)')
            plt.title('Poisson noise 5 Hz per 1000')
            plt.legend()
            plt.show()

    @pytest.mark.skip(reason="plotting")
    def test_theory(self):
        enable_cpp()

        def for_rate(rate):
            pop = MFLinearPopulation("pop", 100, {
                PP.GM: 25. * nS,
                PP.CM: 0.5 * nF,
                PP.VL: -70. * mV,
                PP.VTHR: -50. * mV,
                PP.VRES: -55. * mV,
                PP.TAU_RP: 2. * ms
            })
            pop.rate = 10 * Hz

            noise = MFStaticInput("noise", pop, 1000, rate * Hz, {
                IP.GM: 2 * nS,
                IP.VREV: 0 * volt,
                IP.TAU: 2. * ms,
            })
            pop.add_noise(noise)

            rec = MFLinearInput("rec", pop, {
                IP.GM: 0.973 / 100 * nS,
                IP.VREV: -70 * volt,
                IP.TAU: 10. * ms,
            }, pop)

            system = MFSystem("pop noise rec")
            system.pops += [pop]

            solver = MFSolverRatesVoltages(system, solver='mse', maxiter=1)
            sol = solver.run()
            return sol.state[0]

        rates = np.linspace(1, 25, 25)
        dom = np.array([for_rate(r) for r in rates])
        plt.plot(rates, dom)
        plt.xlabel('Noise rate (Hz) per 1000')
        plt.ylabel('Population rate (Hz) per 100')
        plt.show()
