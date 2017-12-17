"""
2001 Brunel Wang simplified
Effects of Neuromodulation in a Cortical Network Model of Object Working Memory Dominated by Recurrent Inhibition.
Journal of Computational Neuroscience 11, 63-85, 2001.
"""

from brian2 import *

from utils import brian2_introspect

BrianLogger.log_level_debug()

# neurons
N = 1000
N_I = int(N * 0.2)  # interneurons

# voltage
V_L = -70. * mV  # resting
V_thr = -50. * mV
V_reset = -55. * mV
V_E = 0 * mV
V_I = -70. * mV

# membrane capacitance
C_m_I = 0.2 * nF

# membrane leak
g_m_I = 20. * nS

# refactorty period
tau_rp_E = 2. * ms
tau_rp_I = 1. * ms

# external stimuli
rate = 3 * Hz
C_ext = 800

# synapses
C_I = N_I

# AMPA (excitatory)
g_AMPA_ext_E = 2.08 * nS
g_AMPA_ext_I = 1.62 * nS
tau_AMPA = 2. * ms

# GABAergic (inhibitory)
g_GABA_E = 0*1.25 * nS * 200. / N_I
g_GABA_I = 0.01*0.973 * nS * 200. / N_I
tau_GABA = 10. * ms



eqs_I = '''
dv / dt = (- g_m_I * (v - V_L) - I_syn) / C_m_I : volt (unless refractory)

I_syn = I_AMPA_ext + I_GABA_rec : amp

I_AMPA_ext = g_AMPA_ext_I * (v - V_E) * s_AMPA_ext : amp
ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : 1


I_GABA_rec = g_GABA_I * (v - V_I) * s_GABA : amp
ds_GABA / dt = - s_GABA / tau_GABA : 1
'''

P_I = NeuronGroup(N_I, eqs_I, method='euler', threshold='v > V_thr', reset='v = V_reset;', refractory=tau_rp_I)
P_I.v = V_L

eqs_pre_gaba = '''
s_GABA += 1
'''


# recurrent I to I
C_I_I = Synapses(P_I, P_I, on_pre=eqs_pre_gaba)
C_I_I.connect('i != j')

# external
C_P_I = PoissonInput(P_I, 's_AMPA_ext', C_ext, rate, 1)

sp_I = SpikeMonitor(P_I[:10])
r_I = PopulationRateMonitor(P_I)

run(2000 * ms)

magic_network._update_magic_objects(0)
brian2_introspect(magic_network, globals())

suptitle('old')

subplot(211)
title('rates')
plot(r_I.t / ms, r_I.smooth_rate(width=10 * ms) / Hz, label='interneuron')
legend()

subplot(212)
title('spikes')
plot(sp_I.t / ms, sp_I.i, '.', markersize=5, label='interneuron')
legend()

show()
