from meanfield import MFSolver, PP
from meanfield.zoo import brunel_wang_2001_stp
from brian2 import *


def with_g(g_e, g_i):
    # 25, 20 default

    system = brunel_wang_2001_stp.no_subpopulation()
    pop_e, pop_i = system.populations

    pop_e[PP.GM] = g_e * nsiemens
    pop_i[PP.GM] = g_i * nsiemens

    solver = MFSolver.rates_voltages(system, solver='simplex', maxiter=1)
    sol = solver.run()
    print(sol)
    r_e, _, r_i, _ = sol.state
    return [max(0, r_e), max(0, r_i)]


vals = 20
xx = np.linspace(15, 35, vals)
yy = np.linspace(10, 40, vals)


def sim():
    system = brunel_wang_2001_stp.no_subpopulation()
    pop_e, pop_i = system.populations

    r_E = PopulationRateMonitor(pop_e.brian2)
    r_I = PopulationRateMonitor(pop_i.brian2)

    net = system.collect_brian2_network(r_E, r_I)
    net.run(2000 * ms, report='stdout')


    xlabel('ms')
    ylabel('Hz')

    plot(r_E.t / ms, r_E.smooth_rate(width=25 * ms) / Hz, label='E')
    plot(r_I.t / ms, r_I.smooth_rate(width=25 * ms) / Hz, label='I')

    legend()
    show()