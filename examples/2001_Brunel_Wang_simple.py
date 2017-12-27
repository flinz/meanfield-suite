from brian2 import *
from meanfield.parameters import PP
from meanfield.parameters import IP

from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.MFSystem import MFSystem
from meanfield.parameters import Connection
from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.inputs.MFStaticInput import MFStaticInput
from meanfield.utils import brian2_introspect

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
source_e_noise1 = MFStaticInput("E_noise1", pop_e1, C_ext, rate, {
    IP.GM: g_AMPA_ext_E,
    IP.VREV: 0 * mV,
    IP.TAU: tau_AMPA,
})
pop_e1.add_noise(source_e_noise1)

source_e_noise2 = MFStaticInput("E_noise2", pop_e2, C_ext, rate, {
    IP.GM: g_AMPA_ext_E,
    IP.VREV: 0 * mV,
    IP.TAU: tau_AMPA,
})
pop_e2.add_noise(source_e_noise2)

source_i_noise = MFStaticInput("I_noise", pop_i, C_ext, rate, {
    IP.GM: g_AMPA_ext_I,
    IP.VREV: 0 * mV,
    IP.TAU: tau_AMPA,
})
pop_i.add_noise(source_i_noise)

# E->E NMDA
source_ee_nmda1 = MFLinearInput('EE NMDA 1', pop_e1, {
    IP.GM: g_NMDA_E,
    IP.VREV: 0 * mV,
    IP.TAU: tau_NMDA_decay,
}, pop_e1, Connection.all_to_others())

source_ee_nmda2 = MFLinearInput('EE NMDA 2', pop_e2, {
    IP.GM: g_NMDA_E,
    IP.VREV: 0 * mV,
    IP.TAU: tau_NMDA_decay,
}, pop_e1)

# E->E AMPA
source_ee_ampa1 = MFLinearInput('EE AMPA 1', pop_e1, {
    IP.GM: g_AMPA_rec_E,
    IP.VREV: 0 * mvolt,
    IP.TAU: tau_AMPA,
}, pop_e2)

source_ee_ampa2 = MFLinearInput('EE AMPA 2', pop_e2, {
    IP.GM: g_AMPA_rec_E,
    IP.VREV: 0 * mvolt,
    IP.TAU: tau_AMPA,
}, pop_e2, Connection.all_to_others())

# E->I NMDA
source_ie_nmda = MFLinearInput('EI NMDA', pop_i, {
    IP.GM: g_NMDA_I,
    IP.VREV: 0 * mvolt,
    IP.TAU: tau_NMDA_decay,
}, pop_e1)

# E->I AMPA
source_ie_ampa = MFLinearInput('EI AMPA', pop_i, {
    IP.GM: g_AMPA_rec_E,
    IP.VREV: 0 * mvolt,
    IP.TAU: tau_AMPA,
}, pop_e2)

# I->I GABA
source_ii_gaba = MFLinearInput('II GABA', pop_i, {
    IP.GM: g_GABA_I,
    IP.VREV: -70 * mvolt,
    IP.TAU: tau_GABA,
}, pop_i, Connection.all_to_others())

# I->E GABA
source_ie_gaba1 = MFLinearInput('IE GABA 1', pop_e1, {
    IP.GM: g_GABA_E,
    IP.VREV: -70 * mvolt,
    IP.TAU: tau_GABA,
}, pop_i)

source_ie_gaba2 = MFLinearInput('IE GABA 2', pop_e2, {
    IP.GM: g_GABA_E,
    IP.VREV: -70 * mvolt,
    IP.TAU: tau_GABA,
}, pop_i)


system = MFSystem("Brunel Wang simplified")
system.pops += [pop_e1, pop_e2, pop_i]

solver = MFSolverRatesVoltages(system, solver='mse')
sol = solver.run()

sp1 = SpikeMonitor(pop_e1.brian2[:10])
sp2 = SpikeMonitor(pop_e2.brian2[:10])
sp3 = SpikeMonitor(pop_i.brian2[:10])
rate1 = PopulationRateMonitor(pop_e1.brian2)
rate2 = PopulationRateMonitor(pop_e2.brian2)
rate3 = PopulationRateMonitor(pop_i.brian2)
s = StateMonitor(pop_i.brian2, ['v'], record=[100])

net = Network()
net.add(pop_e1.brian2)
net.add(pop_e2.brian2)
net.add(pop_i.brian2)
net.add(source_e_noise1.brian2)
net.add(source_e_noise2.brian2)
net.add(source_i_noise.brian2)
net.add(source_ee_nmda1.brian2)
net.add(source_ee_nmda2.brian2)
net.add(source_ee_ampa1.brian2)
net.add(source_ee_ampa2.brian2)
net.add(source_ie_nmda.brian2)
net.add(source_ie_ampa.brian2)
net.add(source_ii_gaba.brian2)
net.add(source_ie_gaba1.brian2)
net.add(source_ie_gaba2.brian2)
net.add(sp1)
net.add(sp2)
net.add(sp3)
net.add(rate1)
net.add(rate2)
net.add(rate3)
net.add(s)

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
