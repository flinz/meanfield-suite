from time import sleep

from brian2 import *
from meanfield.parameters import NP
from meanfield.parameters import SP

from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFLinearPop import MFLinearPop
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.MFSystem import MFSystem
from meanfield.parameters import Connection
from meanfield.sources.MFLinearSource import MFLinearSource
from meanfield.sources.MFStaticSource import MFStaticSource

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
g_NMDA_E = 0.327 * nS * 800. / N_E
g_NMDA_I = 0.258 * nS * 800. / N_E
tau_NMDA_rise = 2. * ms
tau_NMDA_decay = 100. * ms
alpha = 0. / ms # 0.5 / ms
Mg2 = 0. # 1.

# GABAergic (inhibitory)
g_GABA_E = 0.1 * 1.25 * nS * 200. / N_I # FIXME
g_GABA_I = 0.973 * nS * 200. / N_I
tau_GABA = 10. * ms

# subpopulations
f = 0.1
p = 1
N_sub = int(N_E * f)
N_non = int(N_E * (1. - f * p))
w_plus = 2.1
w_minus = 1. - f * (w_plus - 1.) / (1. - f)


E_params = MFParams({
    #NP.GAMMA: 0.280112,
    #NP.BETA: 0.062,
    NP.GM: g_m_E,
    NP.CM: C_m_E * 1e3,
    NP.VL: V_L,
    NP.VTHR: V_thr,
    NP.VRES: V_reset,
    NP.TAU_RP: tau_rp_E
})

I_params = MFParams({
    #gamma: 0.280112,
    #beta: 0.062,
    NP.GM: g_m_I,
    NP.CM: C_m_I * 1e3,
    NP.VL: V_L,
    NP.VTHR: V_thr,
    NP.VRES: V_reset,
    NP.TAU_RP: tau_rp_I
})

nu_e = 0.01
nu_i = 0.01


# NMDA
pop_e1 = MFLinearPop("E", N_non, E_params)
pop_e1.rate_ms = nu_e * Hz
pop_e1.v_mean = -52. * mV

# AMPA
pop_e2 = MFLinearPop("Edown", N_sub, E_params)
pop_e2.rate_ms = nu_e * Hz
pop_e2.v_mean = -52. * mV

# GABA
pop_i = MFLinearPop("I", N_I, I_params)
pop_i.rate_ms = nu_i * Hz
pop_i.v_mean = -52. * mV


# noise pops
source_e_noise1 = MFStaticSource("E_noise1", pop_e1, C_ext, rate, {
    SP.GM: g_AMPA_ext_E,
    SP.VREV: 0 * volt,
    SP.TAU: tau_AMPA,
})
pop_e1.add_noise(source_e_noise1)

source_e_noise2 = MFStaticSource("E_noise2", pop_e2, C_ext, rate, {
    SP.GM: g_AMPA_ext_E,
    SP.VREV: 0 * volt,
    SP.TAU: tau_AMPA,
})
pop_e2.add_noise(source_e_noise2)

source_i_noise = MFStaticSource("I_noise", pop_i, C_ext, rate, {
    SP.GM: g_AMPA_ext_I,
    SP.VREV: 0 * volt,
    SP.TAU: tau_AMPA,
})
pop_i.add_noise(source_i_noise)

# E->E NMDA
source_ee_nmda1 = MFLinearSource('EE NMDA 1', pop_e1, {
    SP.GM: g_NMDA_E,
    SP.VREV: 0 * volt,
    SP.TAU: tau_NMDA_decay,
}, pop_e1, Connection.all_to_others())

source_ee_nmda2 = MFLinearSource('EE NMDA 2', pop_e2, {
    SP.GM: g_NMDA_E,
    SP.VREV: 0 * volt,
    SP.TAU: tau_NMDA_decay,
}, pop_e1)

# E->E AMPA
source_ee_ampa1 = MFLinearSource('EE AMPA 1', pop_e1, {
    SP.GM: g_AMPA_rec_E,
    SP.VREV: 0 * volt,
    SP.TAU: tau_AMPA,
}, pop_e2)

source_ee_ampa2 = MFLinearSource('EE AMPA 2', pop_e2, {
    SP.GM: g_AMPA_rec_E,
    SP.VREV: 0 * volt,
    SP.TAU: tau_AMPA,
}, pop_e2, Connection.all_to_others())

# E->I NMDA
source_ie_nmda = MFLinearSource('EI NMDA', pop_i, {
    SP.GM: g_NMDA_I,
    SP.VREV: 0 * volt,
    SP.TAU: tau_NMDA_decay,
}, pop_e1)

# E->I AMPA
source_ie_ampa = MFLinearSource('EI AMPA', pop_i, {
    SP.GM: g_AMPA_rec_E,
    SP.VREV: 0 * volt,
    SP.TAU: tau_AMPA,
}, pop_e2)

# I->I GABA
source_ii_gaba = MFLinearSource('II GABA', pop_i, {
    SP.GM: g_GABA_I,
    SP.VREV: -70 * volt,
    SP.TAU: tau_GABA,
}, pop_i, Connection.all_to_others())

# I->E GABA
source_ie_gaba1 = MFLinearSource('IE GABA 1', pop_e1, {
    SP.GM: g_GABA_E,
    SP.VREV: -70 * volt,
    SP.TAU: tau_GABA,
}, pop_i)

source_ie_gaba2 = MFLinearSource('IE GABA 2', pop_e2, {
    SP.GM: g_GABA_E,
    SP.VREV: -70 * volt,
    SP.TAU: tau_GABA,
}, pop_i)


system = MFSystem("Brunel Wang simplified")
system.pops += [pop_e1, pop_e2, pop_i]

solver = MFSolverRatesVoltages(system, solver='gradient')
solver.run()


sp1 = SpikeMonitor(pop_e1.brian2[:40])
sp2 = SpikeMonitor(pop_e2.brian2[:40])
sp3 = SpikeMonitor(pop_i.brian2[:40])
rate1 = PopulationRateMonitor(pop_e1.brian2)
rate2 = PopulationRateMonitor(pop_e2.brian2)
rate3 = PopulationRateMonitor(pop_i.brian2)
s = StateMonitor(pop_i.brian2, ['v'], record=[100])

system.print_introspect()

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



subplot(311)
plot(rate1.t / ms, rate1.smooth_rate(width=25 * ms) / Hz, label='pyramidal')
plot(rate2.t / ms, rate2.smooth_rate(width=25 * ms) / Hz, label='pyramidal')
plot(rate3.t / ms, rate3.smooth_rate(width=25 * ms) / Hz, label='interneuron')
legend()

subplot(312)
plot(sp1.t / ms, sp1.i, '.', markersize=5, label='pyramidal')
plot(sp2.t / ms, sp2.i, '.', markersize=5, label='pyramidal')
plot(sp3.t / ms, sp3.i, '.', markersize=5, label='interneuron')
legend()

subplot(313)
#plot(sp3.t / ms, sp3.i, '.', markersize=5)
plot(s.t / ms, s.v[0] / mV)

show()



# eqs_E = '''
# dv / dt = (- g_m_E * nS * (v - V_L * mV) - I_syn) / (C_m_E * nF) : volt (unless refractory)
#
# I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp
#
# I_AMPA_ext = g_AMPA_ext_E * nS * (v - V_E * mV) * s_AMPA_ext : amp
# I_AMPA_rec = g_AMPA_rec_E * nS * (v - V_E * mV) * s_AMPA : amp
# ds_AMPA_ext / dt = - s_AMPA_ext / (tau_AMPA * ms) : 1
# ds_AMPA / dt = - s_AMPA / (tau_AMPA * ms) : 1
#
# I_NMDA_rec = g_NMDA_E * nS * (v - V_E * mV) * s_NMDA : amp
# ds_NMDA / dt = - s_NMDA / (tau_NMDA_decay * ms) : 1
#
# I_GABA_rec = g_GABA_E * nS * (v - V_I * mV) * s_GABA : amp
# ds_GABA / dt = - s_GABA / (tau_GABA * ms) : 1
# '''
#
# eqs_I = '''
# dv / dt = (- g_m_I * nS * (v - V_L * mV) - I_syn) / (C_m_I * nF) : volt (unless refractory)
#
# I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp
#
# I_AMPA_ext = g_AMPA_ext_I * nS * (v - V_E * mV) * s_AMPA_ext : amp
# I_AMPA_rec = g_AMPA_rec_I * nS * (v - V_E * mV) * s_AMPA : amp
# ds_AMPA_ext / dt = - s_AMPA_ext / (tau_AMPA * ms) : 1
# ds_AMPA / dt = - s_AMPA / (tau_AMPA * ms) : 1
#
# I_NMDA_rec = g_NMDA_I * nS * (v - V_E * mV) * s_NMDA : amp
# ds_NMDA / dt = - s_NMDA / (tau_NMDA_decay * ms) : 1
#
# I_GABA_rec = g_GABA_I * nS * (v - V_I * mV) * s_GABA : amp
# ds_GABA / dt = - s_GABA / (tau_GABA * ms) : 1
# '''
#
# P_E = NeuronGroup(N_E, eqs_E, method='euler', threshold='v > V_thr * mV', reset='v = V_reset * mV', refractory=tau_rp_E * ms)
# P_E.v = V_L * mV
# P_I = NeuronGroup(N_I, eqs_I, method='euler', threshold='v > V_thr * mV', reset='v = V_reset * mV', refractory=tau_rp_I * ms)
# P_I.v = V_L * mV
#
# eqs_glut = '''
# w : 1
# '''
#
# eqs_pre_glut = '''
# s_AMPA += w
# s_NMDA += w
# '''
#
# eqs_pre_gaba = '''
# s_GABA += 1
# '''
#
# eqs_pre_ext = '''
# s_AMPA_ext += 1
# '''
#
# # recurrent E to E
# C_E_E = Synapses(P_E, P_E, method='euler', model=eqs_glut, on_pre=eqs_pre_glut)
# C_E_E.connect('i != j')
# C_E_E.w[:] = 1
#
# for pi in range(N_non, N_non + p * N_sub, N_sub):
#
#     # internal other subpopulation to current nonselective
#     C_E_E.w[C_E_E.indices[:, pi:pi + N_sub]] = w_minus
#
#     # internal current subpopulation to current subpopulation
#     C_E_E.w[C_E_E.indices[pi:pi + N_sub, pi:pi + N_sub]] = w_plus
#
# # E to I
# C_E_I = Synapses(P_E, P_I, method='euler', model=eqs_glut, on_pre=eqs_pre_glut)
# C_E_I.connect()
# C_E_I.w[:] = 1
#
# # recurrent I to I
# C_I_I = Synapses(P_I, P_I, on_pre=eqs_pre_gaba)
# C_I_I.connect('i != j')
#
# # I to E
# C_I_E = Synapses(P_I, P_E, on_pre=eqs_pre_gaba)
# C_I_E.connect()
#
# # external
# C_P_E = PoissonInput(P_E, 's_AMPA_ext', C_ext, rate * Hz, 500)
# C_P_I = PoissonInput(P_I, 's_AMPA_ext', C_ext, rate * Hz, 500)
#
# # monitors
# sp_E = SpikeMonitor(P_E[:40])
# sp_I = SpikeMonitor(P_I[:40])
# r_E = PopulationRateMonitor(P_E)
# r_I = PopulationRateMonitor(P_I)
#
# run(3000 * ms, report='text', report_period=0.5 * second)
#
# subplot(311)
# plot(r_E.t / ms, r_E.smooth_rate(width=25 * ms) / Hz, label='pyramidal neuron')
# plot(r_I.t / ms, r_I.smooth_rate(width=25 * ms) / Hz, label='interneuron')
# legend()
#
# subplot(312)
# plot(sp_E.t / ms, sp_E.i, '.', markersize=5, label='nonselective')
#
# subplot(313)
# plot(sp_I.t / ms, sp_I.i, '.', markersize=5)
#
# show()


