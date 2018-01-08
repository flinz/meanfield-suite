import matplotlib.pyplot as plt
import numpy as np
from brian2 import PopulationRateMonitor, defaultclock
from brian2.units import *

from meanfield.core.MFSystem import MFSystem
from meanfield.inputs.MFStaticInput import MFStaticInput
from meanfield.parameters import IP
from meanfield.parameters import PP
from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.solvers.MFSolver import MFSolver
from meanfield.utils import reset_brian2

n_pop = 25
n_noise = 1000
t = 10 * second
samples = 20

def with_rate(rate):
    reset_brian2()

    pop = MFLinearPopulation(n_pop, {
        PP.GM: 25. * nS,
        PP.CM: 0.5 * nF,
        PP.VL: -70. * mV,
        PP.VTHR: -50. * mV,
        PP.VRES: -55. * mV,
        PP.TAU_RP: 2. * ms
    })
    pop.rate = 10 * Hz

    MFStaticInput(n_noise, rate * Hz, pop, {
        IP.GM: 2.08 * nS,
        IP.VREV: 0 * volt,
        IP.TAU: 1.5 * ms,
    })

    system = MFSystem(pop)

    rate = PopulationRateMonitor(pop.brian2)
    net = system.collect_brian2_network(rate)
    net.run(t)

    margin = round(0.5 * second / defaultclock.dt)
    stable = rate.smooth_rate(width=20 * ms)[margin:-margin]

    mean_simu = np.mean(stable)
    std_simu = np.reshape(stable, (samples, -1)).mean(axis=1).std() / np.sqrt(samples) * 1.96
    print(std_simu)

    solver = MFSolver.rates_voltages(system, solver='simplex', maxiter=1)
    sol = solver.run()
    mean_theo = sol.state[0]

    return [mean_theo, mean_simu, std_simu]


rates = np.linspace(2, 10, 20)
theory, simulation, error = zip(*[with_rate(r) for r in rates])

plt.errorbar(rates, simulation, yerr=error, fmt='.', label='simulation (95%)')
plt.plot(rates, theory, label='theory')

plt.xlabel(f'Noise rate (Hz) from {n_noise}')
plt.ylabel(f'Population rate (Hz) on {n_pop}')
plt.yscale('log')
plt.legend()
plt.show()

