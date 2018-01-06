'''
Sample-specific persistent activity
-----------------------------------

Five subpopulations with three selective and one reset stimuli example.
Analog to figure 6b in the paper.

BRUNEL, Nicolas et WANG, Xiao-Jing. Effects of neuromodulation in a cortical network model of object working memory
dominated by recurrent inhibition. Journal of computational neuroscience, 2001, vol. 11, no 1, p. 63-85.
'''


from brian2 import *

from meanfield.inputs.MFNonLinearNMDAInput import MFNonLinearNMDAInput
from meanfield.MFSystem import MFSystem
from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.inputs.MFStaticInput import MFStaticInput
from meanfield.parameters import Connection
from meanfield.parameters import IP
from meanfield.parameters import PP
from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.utils import reset_brian2
from meanfield.parameters.MFParams import MFParams

BrianLogger.log_level_debug()
reset_brian2()

# populations
N = 250
N_E = int(N * 0.8)  # pyramidal neurons
N_I = int(N * 0.2)  # interneurons

# voltage
V_L = -70. * mV
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

# refractory period
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
alpha = 0.5 / ms
beta = 0.062
gamma = 1. / 3.57
Mg2 = 1.

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

# modeling
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
    IP.VREV: V_E,
    IP.TAU: tau_AMPA,
}, name="E_noise1")

MFStaticInput(C_ext, rate, pop_e2, {
    IP.GM: g_AMPA_ext_E,
    IP.VREV: V_E,
    IP.TAU: tau_AMPA,
}, name="E_noise2")

MFStaticInput(C_ext, rate, pop_i, {
    IP.GM: g_AMPA_ext_I,
    IP.VREV: V_E,
    IP.TAU: tau_AMPA,
}, name="I_noise")

# E->E NMDA
ee_nmda = MFParams({
    IP.GM: g_NMDA_E,
    IP.VREV: V_E,
    IP.TAU: 0 * ms,
    IP.TAU_NMDA: tau_NMDA_decay,
    IP.TAU_NMDA_RISE: tau_NMDA_rise,
    IP.ALPHA: alpha,
    IP.BETA: beta,
    IP.GAMMA: gamma,
    IP.MG: Mg2,
})

MFNonLinearNMDAInput(pop_e1, pop_e1, ee_nmda + {IP.W: 1}, name='EE NMDA 1', connection=Connection.all_to_others())

MFNonLinearNMDAInput(pop_e1, pop_e2, ee_nmda + {IP.W: w_minus}, name='EE NMDA 2')

MFNonLinearNMDAInput(pop_e2, pop_e2, ee_nmda + {IP.W: w_plus}, name='EE NMDA 3', connection=Connection.all_to_others())

MFNonLinearNMDAInput(pop_e2, pop_e1, ee_nmda + {IP.W: 1}, name='EE NMDA 4')


# E->E AMPA
ee_ampa = MFParams({
    IP.GM: g_AMPA_rec_E,
    IP.VREV: V_E,
    IP.TAU: tau_AMPA,
})

MFLinearInput(pop_e1, pop_e1, ee_ampa + {IP.W: 1}, name='EE AMPA 1')

MFLinearInput(pop_e1, pop_e2, ee_ampa + {IP.W: w_minus}, name='EE AMPA 2', connection=Connection.all_to_others())

MFLinearInput(pop_e2, pop_e1, ee_ampa + {IP.W: w_plus}, name='EE AMPA 3')

MFLinearInput(pop_e2, pop_e2, ee_ampa + {IP.W: 1}, name='EE AMPA 4', connection=Connection.all_to_others())

# E->I NMDA
ei_nmda = {
    IP.GM: g_NMDA_I,
    IP.VREV: V_E,
    IP.TAU: 0 * ms,
    IP.W: 1,
    IP.TAU_NMDA: tau_NMDA_decay,
    IP.TAU_NMDA_RISE: tau_NMDA_rise,
    IP.ALPHA: alpha,
    IP.BETA: beta,
    IP.GAMMA: gamma,
    IP.MG: Mg2,
}

MFNonLinearNMDAInput(pop_e1, pop_i, ei_nmda, name='EI NMDA')

MFNonLinearNMDAInput(pop_e2, pop_i, ei_nmda, name='EI NMDA 2')

# E->I AMPA
ei_ampa = {
    IP.GM: g_AMPA_rec_E,
    IP.VREV: V_E,
    IP.TAU: tau_AMPA,
}

MFLinearInput(pop_e1, pop_i, ei_ampa, name='EI AMPA')

MFLinearInput(pop_e2, pop_i, ei_ampa, name='EI AMPA 2')

# I->I GABA
MFLinearInput(pop_i, pop_i, {
    IP.GM: g_GABA_I,
    IP.VREV: V_I,
    IP.TAU: tau_GABA,
}, name='II GABA', connection=Connection.all_to_others())

# I->E GABA
ie_gaba = {
    IP.GM: g_GABA_E,
    IP.VREV: V_I,
    IP.TAU: tau_GABA,
}

MFLinearInput(pop_i, pop_e1, ie_gaba, name='IE GABA 1')

MFLinearInput(pop_i, pop_e2, ie_gaba, name='IE GABA 2')


# monitors
N_activity_plot = 15

sp_E_sel = SpikeMonitor(pop_e2.brian2[:N_activity_plot])
sp_E = SpikeMonitor(pop_e1.brian2[:N_activity_plot])
sp_I = SpikeMonitor(pop_i.brian2[:N_activity_plot])

r_E_sel = PopulationRateMonitor(pop_e2.brian2)
r_E = PopulationRateMonitor(pop_e1.brian2)
r_I = PopulationRateMonitor(pop_i.brian2)

# simulate, can be long >120s
system = MFSystem(pop_e1, pop_e2, pop_i, name="Brunel Wang 2001")

pop_e1.rate = 1.3 * Hz
pop_e2.rate = 1.3 * Hz
pop_i.rate = 6 * Hz

solver = MFSolverRatesVoltages(system, solver='simplex', maxiter=2)
sol = solver.run()

#system.graph().view(cleanup=True)


#print(pop_e1.brian2_model)
#print()
#print(pop_e2.brian2_model)
#print()
#print(pop_i.brian2_model)


net = system.collect_brian2_network(sp_E, sp_I, sp_E_sel, r_E_sel, r_E, r_I)
net.run(2000 * ms, report='stdout')


# plotting
suptitle('Meanfield')
title('Population rates')
xlabel('ms')
ylabel('Hz')

plot(r_E.t / ms, r_E.smooth_rate(width=25 * ms) / Hz, label='nonselective')
plot(r_I.t / ms, r_I.smooth_rate(width=25 * ms) / Hz, label='inhibitory')

plot(r_E_sel.t / ms, r_E_sel.smooth_rate(width=25 * ms) / Hz, label='selective 1')

legend()
show()

suptitle('Meanfield')
title('Population activities ({} neurons/pop)'.format(N_activity_plot))
xlabel('ms')
yticks([])

plot(sp_E.t / ms, sp_E.i + (p + 1) * N_activity_plot, '.', markersize=2, label='nonselective')
plot(sp_I.t / ms, sp_I.i + p * N_activity_plot, '.', markersize=2, label='inhibitory')

plot(sp_E_sel.t / ms, sp_E_sel.i, '.', markersize=2, label='selective 1')

legend()
show()