"""
2001 Brunel Wang simplified
Effects of Neuromodulation in a Cortical Network Model of Object Working Memory Dominated by Recurrent Inhibition.
Journal of Computational Neuroscience 11, 63-85, 2001.
"""

import pandas as pd
from brian2 import *
BrianLogger.log_level_info()

# neurons
N = 1000 / 4
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
g_AMPA_rec_E = 0.104 * nS  * 800. / N_E
g_AMPA_ext_I = 1.62 * nS
g_AMPA_rec_I = 0.081 * nS * 800. / N_E
tau_AMPA = 2. * ms

# NMDA (excitatory)
g_NMDA_E = 0.327 * nS * 800. / N_E
g_NMDA_I = 0.258 * nS * 800. / N_E
tau_NMDA_rise = 2. * ms
tau_NMDA_decay = 100. * ms

# GABAergic (inhibitory)
g_GABA_E = 1.25 * nS * 200. / N_I
g_GABA_I = 0.973 * nS * 200. / N_I
tau_GABA = 10. * ms

# subpopulations
f = 0.1
p = 1
N_sub = int(N_E * f)
N_non = int(N_E * (1. - f * p))
w_plus = 2.1
w_minus = 1. - f * (w_plus - 1.) / (1. - f)

eqs_E = '''
dv / dt = (- g_m_E * (v - V_L) - I_syn) / C_m_E : volt (unless refractory)

I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp

I_AMPA_ext = g_AMPA_ext_E * (v - V_E) * s_AMPA_ext : amp
I_AMPA_rec = g_AMPA_rec_E * (v - V_E) * s_AMPA : amp
ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : 1
ds_AMPA / dt = - s_AMPA / tau_AMPA : 1

I_NMDA_rec = g_NMDA_E * (v - V_E) * s_NMDA : amp
ds_NMDA / dt = - s_NMDA / tau_NMDA_decay : 1

I_GABA_rec = g_GABA_E * (v - V_I) * s_GABA : amp
ds_GABA / dt = - s_GABA / tau_GABA : 1
'''

eqs_I = '''
dv / dt = (- g_m_I * (v - V_L) - I_syn) / C_m_I : volt (unless refractory)

I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp

I_AMPA_ext = g_AMPA_ext_I * (v - V_E) * s_AMPA_ext : amp
I_AMPA_rec = g_AMPA_rec_I * (v - V_E) * s_AMPA : amp
ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : 1
ds_AMPA / dt = - s_AMPA / tau_AMPA : 1

I_NMDA_rec = g_NMDA_I * (v - V_E) * s_NMDA : amp
ds_NMDA / dt = - s_NMDA / tau_NMDA_decay : 1

I_GABA_rec = g_GABA_I * (v - V_I) * s_GABA : amp
ds_GABA / dt = - s_GABA / tau_GABA : 1
'''

P_E = NeuronGroup(N_E, eqs_E, method='euler', threshold='v > V_thr', reset='v = V_reset', refractory=tau_rp_E)
P_E.v = V_L
P_I = NeuronGroup(N_I, eqs_I, method='euler', threshold='v > V_thr', reset='v = V_reset', refractory=tau_rp_I)
P_I.v = V_L

eqs_glut = '''
w : 1
'''

eqs_pre_glut = '''
s_AMPA += w
s_NMDA += w
'''

eqs_pre_gaba = '''
s_GABA += 1
'''

eqs_pre_ext = '''
s_AMPA_ext += 1
'''

# recurrent E to E
C_E_E = Synapses(P_E, P_E, method='euler', model=eqs_glut, on_pre=eqs_pre_glut)
C_E_E.connect('i != j')
C_E_E.w[:] = 1

for pi in range(N_non, N_non + p * N_sub, N_sub):

    # internal other subpopulation to current nonselective
    C_E_E.w[C_E_E.indices[:, pi:pi + N_sub]] = w_minus

    # internal current subpopulation to current subpopulation
    C_E_E.w[C_E_E.indices[pi:pi + N_sub, pi:pi + N_sub]] = w_plus

# E to I
C_E_I = Synapses(P_E, P_I, method='euler', model=eqs_glut, on_pre=eqs_pre_glut)
C_E_I.connect()
C_E_I.w[:] = 1

# recurrent I to I
C_I_I = Synapses(P_I, P_I, on_pre=eqs_pre_gaba)
C_I_I.connect('i != j')

# I to E
C_I_E = Synapses(P_I, P_E, on_pre=eqs_pre_gaba)
C_I_E.connect()

# external
C_P_E = PoissonInput(P_E, 's_AMPA_ext', C_ext, rate, 1)
C_P_I = PoissonInput(P_I, 's_AMPA_ext', C_ext, rate, 1)

#stimuli = TimedArray(np.r_[np.zeros(10), np.ones(2) * 50, np.zeros(50)] * Hz, dt=100 * ms)
#C_PG_sti = PoissonGroup(100, 'stimuli(t)')
#C_PG_E = Synapses(C_PG_sti, P_E[N_non: N_non + N_sub], on_pre='s_AMPA += 1')
#C_PG_E.connect()

# monitors

r_E = SpikeMonitor(P_E, record=False)
r_I = SpikeMonitor(P_I, record=False)

store()

resolution = 25
ges = np.linspace(0.1, 0.8, resolution)
gis = np.linspace(0.1, 0.6, resolution)
rates = []

for ge in ges:
    for gi in gis:
        print('Run %f/%f' % (ge, gi))
        restore()

        g_NMDA_E = ge * nS * 800. / N_E
        g_GABA_E = 4. * g_NMDA_E
        g_NMDA_I = gi * nS * 800. / N_E
        g_GABA_I = 4. * g_NMDA_I

        run(500 * ms, report='text', report_period=0.5 * second)

        #plot(r_E.t / ms, r_E.smooth_rate(width=10 * ms) / Hz)
        #plot(r_I.t / ms, r_I.smooth_rate(width=10 * ms) / Hz)
        #show()

        rates.append({
            'g_E': ge,
            'g_I': gi,
            'rate_E': r_E.num_spikes / N_E / 2.,
            'rate_I': r_I.num_spikes / N_I / 2.
        })

rates = pd.DataFrame(rates)
pd.to_pickle(rates, '2001_Brunel_Wang_simplified_mean_%d.p' % resolution)

rates_E = rates.pivot_table(index='g_I', columns='g_E', values='rate_E')
rates_I = rates.pivot_table(index='g_I', columns='g_E', values='rate_I')

subplot(121)
title('excitatory')
imshow(rates_E, aspect='auto', interpolation='nearest', extent=[rates_E.columns.values[0], rates_E.columns.values[-1], rates_E.index.values[0], rates_E.index.values[-1]], origin="lower")
colorbar()
print(rates_E)

subplot(122)
title('inhibitory')
imshow(rates_I, aspect='auto', interpolation='nearest', extent=[rates_I.columns.values[0], rates_I.columns.values[-1], rates_I.index.values[0], rates_I.index.values[-1]], origin="lower")
colorbar()
print(rates_I)

show()
