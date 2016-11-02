from brian2 import *

BrianLogger.log_level_debug()

N = 1000/10
N_E = int(N * 0.8) # pyramidals
N_I = int(N * 0.2) # interneurons

C_E = N_E #800
C_ext = 800

V_L = -70. * mV # resting
V_thr = -50. * mV
V_reset = -55. * mV

C_m_E = 0.5 * nF
g_m_E = 25. * nS
tau_rp_E = 2. * ms

g_AMPA_ext_E = 2.08 * nS
g_AMPA_rec_E = 0.104 * nS
g_NMDA_E = 0.327 * nS
g_GABA_E = 1.25 * nS

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

eqs = '''
dv / dt = (- g_m_E * (v - V_L) - I_syn) / C_m_E : volt (unless refractory)

I_syn = I_AMPA_ext : amp # + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp

I_AMPA_ext = g_AMPA_ext_E * (v - V_E) * 1 * s_AMPA : amp
I_AMPA_rec = g_AMPA_rec_E * (v - V_E) * s_AMPA : amp
ds_AMPA / dt = - s_AMPA / tau_AMPA : 1

I_NMDA_rec = g_NMDA_E * (v - V_E) / (1 + Mg2 * exp(-0.062 * v / volt) / 3.57) * s_NMDA_tot : amp
s_NMDA_tot : 1

I_GABA_rec = g_GABA_E * (v - V_I) * s_GABA : amp
ds_GABA / dt = - s_GABA / tau_GABA : 1
'''

eqs_glut = '''
s_NMDA_tot_post = 1 * s_NMDA : 1 (summed)
ds_NMDA / dt = - s_NMDA / tau_NMDA_decay + alpha * x * (1 - s_NMDA) : 1 (clock-driven)
dx / dt = - x / tau_NMDA_rise : 1 (clock-driven)
'''

eqs_pre_ext = '''s_AMPA += 1'''
eqs_pre_glut = '''s_AMPA += 1'''
eqs_pre_gaba = '''s_GABA += 1'''

P_E = NeuronGroup(N_E, eqs, method='rk4', threshold='v > V_thr', reset='v = V_reset;', refractory=tau_rp_E)

P_E.v = V_L

# NMDA + AMPA
C_E_E = Synapses(P_E, P_E, model=eqs_glut,on_pre=eqs_pre_glut)
C_E_E.connect(condition='i != j')
#C_E_E.w[:] = 1

# AMPA + NMDA
C_P_E = PoissonInput(P_E, 's_AMPA', C_ext, rate, 1)

S_E = StateMonitor(P_E, ['s_AMPA', 's_GABA', 's_NMDA_tot'], record=True)

M_E = SpikeMonitor(P_E)

R_E = PopulationRateMonitor(P_E)

run(5000 * ms)

plot(R_E.t / ms, R_E.rate / Hz / C_E)

show()

subplot(221)
title('pyramidal voltage (10)')
for i in range(10):
    #plot(S_E.t / ms, S_E.s_GABA[i])
    #plot(S_E.t / ms, S_E.s_AMPA[i])
    plot(S_E.t / ms, S_E.s_AMPA[i])

subplot(222)
title('interneuron voltage (10)')
for i in range(10):
    #plot(S_I.t / ms, S_I.s_GABA[i])
    #plot(S_I.t / ms, S_I.s_AMPA[i])
    plot(S_E.t / ms, S_E.s_NMDA_tot[i])

subplot(223)
title('pyramidal spikes')
plot(M_E.t / ms, M_E.i, '.')

subplot(224)
title('interneuron spikes')

show()

