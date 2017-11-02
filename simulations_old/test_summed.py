from brian2 import *
BrianLogger.log_level_debug()

V_L = -70. * mV
V_thr = -50. * mV
V_reset = -55. * mV
V_E = 0. * mV

C_m_E = 0.5 * nF
g_m_E = 25. * nS
tau_rp_E = 2. * ms
rate = 3 * Hz
C_ext = 800

g_AMPA_ext_E = 2.08 * nS
g_AMPA_rec_E = 0.104 * nS  * 800. / 2.
g_AMPA_ext_I = 1.62 * nS
g_AMPA_rec_I = 0.081 * nS * 800. / 2.
tau_AMPA = 2. * ms

g_NMDA_E = 0.327 * nS  * 800. / 2.
g_NMDA_I = 0.258 * nS * 800. / 2.
tau_NMDA_rise = 2. * ms
tau_NMDA_decay = 100. * ms
alpha = 0.5 / ms
Mg2 = 1.

eqs_E = '''
dv / dt = (- g_m_E * (v - V_L) - I_syn) / C_m_E : volt (unless refractory)

I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec : amp

I_AMPA_ext = g_AMPA_ext_E * (v - V_E) * s_AMPA_ext : amp
I_AMPA_rec = g_AMPA_rec_E * (v - V_E) * 1 * s_AMPA : amp
ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : 1
ds_AMPA / dt = - s_AMPA / tau_AMPA : 1

I_NMDA_rec = g_NMDA_E * (v - V_E) / (1 + Mg2 * exp(-0.062 * v / mV) / 3.57) * s_NMDA_tot : amp
s_NMDA_tot : 1
'''

n = NeuronGroup(1, eqs_E, method='rk4', threshold='v > V_thr', reset='v = V_reset;', refractory=tau_rp_E)
n.v = V_L

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

indices = array([0, 1, 0, 1])
times = array([10, 500, 1000, 1000]) * ms
sg = SpikeGeneratorGroup(2, indices, times)

C = Synapses(sg, n, method='rk4', model=eqs_glut, on_pre=eqs_pre_glut)
C.connect()
C.w[0] = 2.
C.w[1] = 1.

sc1 = StateMonitor(C, ['s_NMDA', 'x', 'w'], record=True)
s = StateMonitor(n, ['s_AMPA', 's_AMPA_ext', 's_NMDA_tot'], record=True)
sp = SpikeMonitor(n)
sp_sg = SpikeMonitor(sg)

run(2000 * ms)

subplot(221)
plot(sc1.t / ms, sc1.s_NMDA[0], label='nmda')
plot(sc1.t / ms, sc1.x[0], label='x')
plot(sc1.t / ms, sc1.w[0], label='w')
legend()

subplot(223)
plot(s.t / ms, s.s_AMPA_ext[0], label="ext")
plot(s.t / ms, s.s_AMPA[0], label="ampa")
plot(s.t / ms, s.s_NMDA_tot[0], label="nmda")
legend()

subplot(224)
plot(sp.t / ms, sp.i + 0.25, '.', markersize=5, label='neuron')
plot(sp_sg.t / ms, sp_sg.i / 4. + 0.5, 'r.', markersize=5, label='generator')
ylim(0, 1)
legend()

show()
