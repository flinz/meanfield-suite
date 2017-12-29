from brian2 import *

from meanfield.MFSystem import MFSystem
from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.inputs.MFStaticInput import MFStaticInput
from meanfield.parameters import Connection
from meanfield.parameters import IP
from meanfield.parameters import PP
from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.solvers.MFSolver import MFSolverRatesVoltages

BrianLogger.log_level_debug()
set_device('cpp_standalone')

# neurons
N = 1000
N_E = int(N * 0.8)  # pyramidal neurons
N_I = int(N * 0.2)  # interneurons

# voltage
V_L = -70. * mV  # resting
V_thr = -50. * mV
V_reset = -55. * mV
V_E = 0. * mV
V_I = -70. * mV

# membrane capacitance
C_m_E = 0.5 * nF
C_m_I = 0.2 * nF

# membrane leak
g_m_E = 25. * nS
g_m_I = 20. * nS

# refactorty period
tau_rp_E = 2. * ms
tau_rp_I = 1. * ms

# external stimuli
rate = 3 * Hz
C_ext = 800

# synapses
C_E = N_E
C_I = N_I

# AMPA (excitatory)
g_AMPA_ext_E = 2.08 * nS
g_AMPA_rec_E = 0.104 * nS * 800. / N_E
g_AMPA_ext_I = 1.62 * nS
g_AMPA_rec_I = 0.081 * nS * 800. / N_E
tau_AMPA = 2. * ms

# NMDA (excitatory)
g_NMDA_E = 0.01 * 0.327 * nS * 800. / N_E
g_NMDA_I = 0.01 * 0.258 * nS * 800. / N_E
tau_NMDA_rise = 2. * ms
tau_NMDA_decay = 100. * ms
alpha = 0. / ms # 0.5 / ms
Mg2 = 0. # 1.

# GABAergic (inhibitory)
g_GABA_E = 1.25 * nS * 200. / N_I
g_GABA_I = 0.973 * nS * 200. / N_I
tau_GABA = 10. * ms

# subpopulations
f = 0.1
p = 1
N_sub = int(N_E * f)
N_non = int(N_E * (1. - f * p))
w_plus = 1#2.1
w_minus = 1. - f * (w_plus - 1.) / (1. - f)

E_params = {
    PP.GM: g_m_E,
    PP.CM: C_m_E,
    PP.VL: V_L,
    PP.VTHR: V_thr,
    PP.VRES: V_reset,
    PP.TAU_RP: tau_rp_E
}

I_params = {
    PP.GM: g_m_I,
    PP.CM: C_m_I,
    PP.VL: V_L,
    PP.VTHR: V_thr,
    PP.VRES: V_reset,
    PP.TAU_RP: tau_rp_I,
}

# NMDA
pop_e1 = MFLinearPopulation(N_non, E_params, name="E")

# AMPA
pop_e2 = MFLinearPopulation(N_sub, E_params, name="Edown")

# GABA
pop_i = MFLinearPopulation(N_I, I_params, name="I")


# noise pops
MFStaticInput(C_ext, rate, pop_e1, {
    IP.GM: g_AMPA_ext_E,
    IP.VREV: 0 * mV,
    IP.TAU: tau_AMPA,
}, name="E_noise1")

MFStaticInput(C_ext, rate, pop_e2, {
    IP.GM: g_AMPA_ext_E,
    IP.VREV: 0 * mV,
    IP.TAU: tau_AMPA,
}, name="E_noise2")

MFStaticInput(C_ext, rate, pop_i, {
    IP.GM: g_AMPA_ext_I,
    IP.VREV: 0 * mV,
    IP.TAU: tau_AMPA,
}, name="I_noise")

# E->E NMDA
MFLinearInput(pop_e1, pop_e1, {
    IP.GM: g_NMDA_E,
    IP.VREV: 0 * mV,
    IP.TAU: tau_NMDA_decay,
}, name='EE NMDA 1', connection=Connection.all_to_others())

MFLinearInput(pop_e1, pop_e2, {
    IP.GM: g_NMDA_E,
    IP.VREV: 0 * mV,
    IP.TAU: tau_NMDA_decay,
}, name='EE NMDA 2')

# E->E AMPA
MFLinearInput(pop_e2, pop_e1, {
    IP.GM: g_AMPA_rec_E,
    IP.VREV: 0 * mvolt,
    IP.TAU: tau_AMPA,
}, name='EE AMPA 1')

MFLinearInput(pop_e2, pop_e2, {
    IP.GM: g_AMPA_rec_E,
    IP.VREV: 0 * mvolt,
    IP.TAU: tau_AMPA,
}, name='EE AMPA 2', connection=Connection.all_to_others())

# E->I NMDA
MFLinearInput(pop_e1, pop_i, {
    IP.GM: g_NMDA_I,
    IP.VREV: 0 * mvolt,
    IP.TAU: tau_NMDA_decay,
}, name='EI NMDA')

# E->I AMPA
MFLinearInput(pop_e2, pop_i, {
    IP.GM: g_AMPA_rec_E,
    IP.VREV: 0 * mvolt,
    IP.TAU: tau_AMPA,
}, name='EI AMPA')

# I->I GABA
MFLinearInput(pop_i, pop_i, {
    IP.GM: g_GABA_I,
    IP.VREV: -70 * mvolt,
    IP.TAU: tau_GABA,
}, name='II GABA', connection=Connection.all_to_others())

# I->E GABA
MFLinearInput(pop_i, pop_e1, {
    IP.GM: g_GABA_E,
    IP.VREV: -70 * mvolt,
    IP.TAU: tau_GABA,
}, name='IE GABA 1')

MFLinearInput(pop_i, pop_e2, {
    IP.GM: g_GABA_E,
    IP.VREV: -70 * mvolt,
    IP.TAU: tau_GABA,
}, name='IE GABA 2')


system = MFSystem(pop_e1, pop_e2, pop_i, name="Brunel Wang simplified")

solver = MFSolverRatesVoltages(system, solver='mse')
sol = solver.run()

sp1 = SpikeMonitor(pop_e1.brian2[:10])
sp2 = SpikeMonitor(pop_e2.brian2[:10])
sp3 = SpikeMonitor(pop_i.brian2[:10])
rate1 = PopulationRateMonitor(pop_e1.brian2)
rate2 = PopulationRateMonitor(pop_e2.brian2)
rate3 = PopulationRateMonitor(pop_i.brian2)
s = StateMonitor(pop_i.brian2, ['v'], record=[100])


net = system.collect_brian2_network(sp1, sp2, sp3, rate1, rate2, rate3, s)
net.run(2000 * ms)

#brian2_introspect(net, globals())
#print()
#system.print_introspect()


if True:

    #suptitle('new')

    #subplot(211)
    plot(rate1.t / ms, rate1.smooth_rate(width=50 * ms) / Hz, label='simulation pyramidal')
    plot(rate2.t / ms, rate2.smooth_rate(width=50 * ms) / Hz, label='simulation pyramidal')
    plot(rate3.t / ms, rate3.smooth_rate(width=50 * ms) / Hz, label='simulation interneuron')
    plot(ones(2000) * sol.state[0], label='theory pyramidal')
    plot(ones(2000) * sol.state[1], label='theory pyramidal')
    plot(ones(2000) * sol.state[2], label='theory interneuron')
    plt.ylabel('Population rate (Hz)')
    plt.xlabel('Simulation time (ms)')
    legend()

    #subplot(212)
    #title('spikes')
    #plot(sp1.t / ms, sp1.i, '.', markersize=5, label='pyramidal')
    #plot(sp2.t / ms, sp2.i, '.', markersize=5, label='pyramidal')
    #plot(sp3.t / ms, sp3.i, '.', markersize=5, label='interneuron')
    #legend()

    show()
