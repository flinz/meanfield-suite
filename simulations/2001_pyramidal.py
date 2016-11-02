from brian2 import *

BrianLogger.log_level_debug()

N = 1000 / 5
N_E = int(N * 0.8) # pyramidals
N_I = int(N * 0.2) # interneurons

C_E = N_E #800
C_I = N_I #200
C_ext = 800

V_L = -70. * mV # resting
V_thr = -50. * mV
V_reset = -55. * mV

C_m_E = 0.5 * nF
C_m_I = 0.2 * nF

g_m_E = 25. * nS
g_m_I = 20. * nS

tau_rp_E = 2. * ms
tau_rp_I = 1. * ms

#tau_m_E = C_m_E / g_m_E
#tau_m_I = C_m_I / g_m_I

g_AMPA_ext_E = 2.08 * nS
g_AMPA_rec_E = 0.104 * nS  * 800. / N_E
g_NMDA_E = 0.327 * nS  * 800. / N_E
g_GABA_E = 1.25 * nS  * 200. / N_I

g_AMPA_ext_I = 1.62 * nS
g_AMPA_rec_I = 0.081 * nS * 800. / N_E
g_NMDA_I = 0.258 * nS * 800. / N_E
g_GABA_I = 0.973 * nS * 200. / N_I

V_E = 0. * mV # TODO expected?
V_I = -70. * mV
Mg2 = 1.

tau_AMPA = 2. * ms
tau_NMDA_decay = 100. * ms
tau_NMDA_rise = 2. * ms
tau_GABA = 10. * ms
alpha = 0.5 / ms

f = 0.1
p = 5
N_p = int(N_E * f)
w_plus = 2.1
w_minus = 1 - f * (w_plus - 1) / (1 - f)
rate = 3 * Hz

eqs_E = '''
dv / dt = (- g_m_E * (v - V_L) - I_syn) / C_m_E : volt (unless refractory)

I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp

I_AMPA_ext = g_AMPA_ext_E * (v - V_E) * s_AMPA_ext : amp
I_AMPA_rec = g_AMPA_rec_E * (v - V_E) * 1 * s_AMPA : amp
ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : 1
ds_AMPA / dt = - s_AMPA / tau_AMPA : 1

I_NMDA_rec = g_NMDA_E * (v - V_E) / (1 + Mg2 * exp(-0.062 * v / mV) / 3.57) * s_NMDA_tot : amp
s_NMDA_tot : 1

I_GABA_rec = g_GABA_E * (v - V_I) * s_GABA : amp
ds_GABA / dt = - s_GABA / tau_GABA : 1
'''

eqs_I = '''
dv / dt = (- g_m_I * (v - V_L) - I_syn) / C_m_I : volt (unless refractory)

I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp

I_AMPA_ext = g_AMPA_ext_I * (v - V_E) * s_AMPA_ext : amp
I_AMPA_rec = g_AMPA_rec_I * (v - V_E) * 1 * s_AMPA : amp
ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : 1
ds_AMPA / dt = - s_AMPA / tau_AMPA : 1

I_NMDA_rec = g_NMDA_I * (v - V_E) / (1 + Mg2 * exp(-0.062 * v / mV) / 3.57) * s_NMDA_tot : amp
s_NMDA_tot : 1

I_GABA_rec = g_GABA_I * (v - V_I) * s_GABA : amp
ds_GABA / dt = - s_GABA / tau_GABA : 1
'''

eqs_glut = '''
s_NMDA_tot_post = w * s_NMDA : 1 (summed)
ds_NMDA / dt = - s_NMDA / tau_NMDA_decay + alpha * x * (1 - s_NMDA) : 1 (clock-driven)
dx / dt = - x / tau_NMDA_rise : 1 (clock-driven)
w : 1
'''

eqs_pre_glut = '''
s_AMPA += w
x += 1
'''

eqs_pre_gaba = '''
s_GABA += 1
'''

eqs_pre_ext = '''
s_AMPA_ext += 1
'''

P_E = NeuronGroup(N_E, eqs_E, method='rk4', threshold='v > V_thr', reset='v = V_reset;', refractory=tau_rp_E)
P_I = NeuronGroup(N_I, eqs_I, method='rk4', threshold='v > V_thr', reset='v = V_reset;', refractory=tau_rp_I)

P_E.v = V_L
P_I.v = V_L

# irresponsive = P_E[:N_p]
# C_E_rec = Synapses(
#     irresponsive,
#     irresponsive,
#     on_pre=eqs_pre
# )
# C_E_rec.connect(condition='i != j')
#
# C_E_rec2 = Synapses(
#     P_E[N_p:],
#     irresponsive,
#     on_pre=eqs_pre
# )
# C_E_rec2.connect()
#
# for pi in range(p):
#     pi = (pi + 1) * N_p
#     group = P_E[pi:pi + N_p]
#     next = P_E[pi + N_p:]
#
#     C_E_rec3 = Synapses(
#         group,
#         group,
#         on_pre='''
#         w = w_plus
#         ''' + eqs_pre
#     )
#     C_E_rec3.connect(condition='i != j')
#
#     C_E_rec4 = Synapses(
#         next,
#         group,
#         on_pre='''
#         w = w_minus
#         ''' + eqs_pre
#     )
#    C_E_rec4.connect()

# # NMDA + AMPA
C_E_E = Synapses(P_E, P_E, model=eqs_glut, on_pre=eqs_pre_glut)
C_E_E.connect(condition='i != j')
C_E_E.w[:] = 1
#
# NMDA + AMPA
C_E_I = Synapses(P_E, P_I, model=eqs_glut, on_pre=eqs_pre_glut)
C_E_I.connect()
C_E_I.w[:] = 1
#
# # GABA
C_I_E = Synapses(P_I, P_E, on_pre=eqs_pre_gaba)
C_I_E.connect()

# # GABA
C_I_I = Synapses(P_I, P_I, on_pre=eqs_pre_gaba)
C_I_I.connect(condition='i != j')

# AMPA + NMDA
C_P_E = PoissonInput(P_E, 's_AMPA_ext', C_ext, rate, 1)
C_P_I = PoissonInput(P_I, 's_AMPA_ext', C_ext, rate, 1)
# P_PE = PoissonGroup(C_ext, rate)
# C_P_E = Synapses(P_PE, P_E, on_pre=eqs_pre_ext)
# C_P_E.connect()
# C_P_E.w[:] = w_minus
#
# P_PI = PoissonGroup(C_ext, rate)
# C_P_I = Synapses(P_PI, P_I, on_pre=eqs_pre_ext)
# C_P_I.connect()
# C_P_I.w[:] = w_minus


S_E = StateMonitor(P_E[:10], ['s_AMPA', 's_GABA', 's_AMPA_ext', 's_NMDA_tot'], record=True)
S_I = StateMonitor(P_I[:10], ['s_AMPA', 's_GABA', 's_AMPA_ext', 's_NMDA_tot'], record=True)

M_E = SpikeMonitor(P_E[:10])
M_I = SpikeMonitor(P_I[:10])

R_E = PopulationRateMonitor(P_E)
R_I = PopulationRateMonitor(P_I)


run(3000 * ms)

subplot(1,2,1)
plot(R_E.t / ms, R_E.smooth_rate(width=10 * ms) / Hz)
subplot(1,2,2)
plot(R_I.t / ms, R_I.smooth_rate(width=10 * ms) / Hz)

show()

subplot(221)
title('pyramidal voltage (10)')
for i in range(1):
    plot(S_E.t / ms, S_E.s_GABA[i], label="gaba")
    plot(S_E.t / ms, S_E.s_AMPA_ext[i], label="ext")
    plot(S_E.t / ms, S_E.s_AMPA[i], label="ampa")
    plot(S_E.t / ms, S_E.s_NMDA_tot[i], label="nmda")
legend()

subplot(222)
title('interneuron voltage (10)')
for i in range(1):
    plot(S_I.t / ms, S_I.s_GABA[i], label="gaba")
    plot(S_I.t / ms, S_I.s_AMPA_ext[i], label="ext")
    plot(S_I.t / ms, S_I.s_AMPA[i], label="ampa")
    plot(S_I.t / ms, S_I.s_NMDA_tot[i], label="nmda")
legend()

subplot(223)
title('pyramidal spikes')
plot(M_E.t / ms, M_E.i, '.', markersize=5)

subplot(224)
title('interneuron spikes')
plot(M_I.t / ms, M_I.i, '.', markersize=5)

show()

# TODO
# no linear decomposition
