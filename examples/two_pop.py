from brian2 import *
from meanfield.parameters import PP
from meanfield.parameters import IP

from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.MFSystem import MFSystem
from meanfield.inputs.MFStaticInput import MFStaticInput

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
g_AMPA_rec_E = 0*0.104 * nS * 800. / N_E
g_AMPA_ext_I = 1.62 * nS
g_AMPA_rec_I = 0*0.081 * nS * 800. / N_E
tau_AMPA = 2. * ms

# NMDA (excitatory)
g_NMDA_E = 0 * 0.327 * nS * 800. / N_E
g_NMDA_I = 0 * 0.258 * nS * 800. / N_E
tau_NMDA_rise = 2. * ms
tau_NMDA_decay = 100. * ms
alpha = 0. / ms # 0.5 / ms
Mg2 = 0. # 1.

# GABAergic (inhibitory)
g_GABA_E = 0*1.25 * nS * 200. / N_I
g_GABA_I = 0*0.973 * nS * 200. / N_I
tau_GABA = 10. * ms

# subpopulations
f = 0.1
p = 1
N_sub = int(N_E * f)
N_non = int(N_E * (1. - f * p))
w_plus = 2.1
w_minus = 1. - f * (w_plus - 1.) / (1. - f)

nu_e = 0.01
nu_i = 0.01


exc = MFLinearPopulation("exc", N_non, {
    PP.GM: g_m_E,
    PP.CM: C_m_E,
    PP.VL: V_L,
    PP.VTHR: V_thr,
    PP.VRES: V_reset,
    PP.TAU_RP: tau_rp_E
})
exc.rate_ms = nu_e * Hz
exc.v_mean = -52. * mV

ini = MFLinearPopulation("ini", N_I, {
    PP.GM: g_m_I,
    PP.CM: C_m_I,
    PP.VL: V_L,
    PP.VTHR: V_thr,
    PP.VRES: V_reset,
    PP.TAU_RP: tau_rp_I
})
ini.rate_ms = nu_i * Hz
ini.v_mean = -52. * mV


# noise

source_e_noise = MFStaticInput("E noise", exc, C_ext, rate, {
    IP.GM: g_AMPA_ext_E,
    IP.VREV: 0 * volt,
    IP.TAU: tau_AMPA,
})
exc.add_noise(source_e_noise)

source_i_noise = MFStaticInput("I noise", ini, C_ext, rate, {
    IP.GM: g_AMPA_ext_I,
    IP.VREV: 0 * volt,
    IP.TAU: tau_AMPA,
})
ini.add_noise(source_i_noise)


# I->I

#source_ii_gaba = MFLinearSource('II GABA', ini, {
#    SP.GM: g_GABA_I,
#    SP.VREV: -70 * volt,
#    SP.TAU: tau_GABA,
#}, ini, Connection.all_to_others())

# I->E

#source_ie_gaba1 = MFLinearSource('IE GABA', exc, {
#    SP.GM: g_GABA_E,
#    SP.VREV: -70 * volt,
#    SP.TAU: tau_GABA,
#}, ini)


system = MFSystem("Two pop")
system.pops += [exc]

solver = MFSolverRatesVoltages(system, solver='mse')
print(solver.state)
solver.run()

print(exc.brian2_model())


sp1 = SpikeMonitor(exc.brian2[:40])
sp2 = SpikeMonitor(ini.brian2[:40])
rate1 = PopulationRateMonitor(exc.brian2)
rate2 = PopulationRateMonitor(ini.brian2)
s1 = StateMonitor(ini.brian2, ['v'], record=[100])
s2 = StateMonitor(ini.brian2, ['v'], record=[100])

system.print_introspect()

net = Network()
net.add(exc.brian2)
net.add(ini.brian2)
net.add(source_e_noise.brian2)
net.add(source_i_noise.brian2)
#net.add(source_ii_gaba.brian2)
#net.add(source_ie_gaba1.brian2)
net.add(sp1)
net.add(sp2)
net.add(rate1)
net.add(rate2)
net.add(s1)
net.add(s2)
net.run(2000 * ms)



subplot(311)
plot(rate1.t / ms, rate1.smooth_rate(width=25 * ms) / Hz, label='pyramidal')
plot(rate2.t / ms, rate2.smooth_rate(width=25 * ms) / Hz, label='interneuron')
legend()

subplot(312)
plot(sp1.t / ms, sp1.i, '.', markersize=5, label='pyramidal')
plot(sp2.t / ms, sp2.i, '.', markersize=5, label='interneuron')
legend()

subplot(313)
plot(s1.t / ms, s1.v[0] / mV, label='pyramidal')
plot(s2.t / ms, s2.v[0] / mV, label='interneuron')
legend()

show()
