'''
Sample-specific persistent activity
-----------------------------------

Five subpopulations with three selective and one reset stimuli example.
Analog to figure 6b in the paper.

BRUNEL, Nicolas et WANG, Xiao-Jing. Effects of neuromodulation in a cortical network model of object working memory
dominated by recurrent inhibition. Journal of computational neuroscience, 2001, vol. 11, no 1, p. 63-85.
'''

import matplotlib.pyplot as plt
import numpy as np
from brian2 import BrianLogger, TimedArray, PoissonInput
from brian2.units import *

from meanfield import plots, modelling
from core.MFSystem import MFSystem
from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.inputs.MFNonLinearNMDAInput import MFNonLinearNMDAInput
from meanfield.inputs.MFStaticInput import MFStaticInput
from meanfield.parameters import Connection
from meanfield.parameters import IP
from meanfield.parameters import PP
from meanfield.parameters.MFParameters import MFParameters
from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.solvers.MFSolver import MFSolver
from meanfield.utils import reset_brian2

BrianLogger.log_level_debug()
reset_brian2()

# populations
N = 1000
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
p = 5
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

pop_e = MFLinearPopulation(N_non, E_params, name="E")
pop_e_sel = modelling.multiple_populations(p, MFLinearPopulation, N_sub, E_params, name="E sel")

#pop_e_el = MFLinearPopulation(N_sub, E_params, name="E sel")

I_params = {
    PP.GM: g_m_I,
    PP.CM: C_m_I,
    PP.VL: V_L,
    PP.VTHR: V_thr,
    PP.VRES: V_reset,
    PP.TAU_RP: tau_rp_I,
}

pop_i = MFLinearPopulation(N_I, I_params, name="I")

# noise pops
e_ampa = {
    IP.GM: g_AMPA_ext_E,
    IP.VREV: V_E,
    IP.TAU: tau_AMPA,
}

modelling.noise_connected(MFStaticInput, C_ext, rate, [pop_e] + pop_e_sel, e_ampa, name='noise')

i_ampa = {
    IP.GM: g_AMPA_ext_I,
    IP.VREV: V_E,
    IP.TAU: tau_AMPA,
}

modelling.noise_connected(MFStaticInput, C_ext, rate, pop_i, i_ampa, name='noise')

# E->E NMDA
ee_nmda = MFParameters({
    IP.GM: g_NMDA_E,
    IP.VREV: V_E,
    IP.TAU: tau_NMDA_decay,
    IP.TAU_NMDA_RISE: tau_NMDA_rise,
    IP.ALPHA: alpha,
    IP.BETA: beta,
    IP.GAMMA: gamma,
    IP.MG: Mg2,
    IP.W: 1,
})

modelling.fully_connected(MFNonLinearNMDAInput, pop_e, pop_e, ee_nmda, name='NMDA')
modelling.fully_connected(MFNonLinearNMDAInput, pop_e, pop_e_sel, ee_nmda, name='NMDA')
modelling.fully_connected(MFNonLinearNMDAInput, pop_e_sel, pop_e, ee_nmda, name='NMDA')
modelling.fully_connected(MFNonLinearNMDAInput, pop_e_sel, pop_e_sel,
                          parameters=ee_nmda + {IP.W: w_minus},
                          self_loop_parameters=ee_nmda + {IP.W: w_plus},
                          name='NMDA')

# E->E AMPA
ee_ampa = MFParameters({
    IP.GM: g_AMPA_rec_E,
    IP.VREV: V_E,
    IP.TAU: tau_AMPA,
})

modelling.fully_connected(MFLinearInput, pop_e, pop_e, ee_ampa, name='AMPA')
modelling.fully_connected(MFLinearInput, pop_e, pop_e_sel, ee_ampa, name='AMPA')
modelling.fully_connected(MFLinearInput, pop_e_sel, pop_e, ee_ampa, name='AMPA')
modelling.fully_connected(MFLinearInput, pop_e_sel, pop_e_sel,
                          parameters=ee_ampa + {IP.W: w_minus},
                          self_loop_parameters=ee_ampa + {IP.W: w_plus},
                          name='AMPA')

# E->I NMDA
ei_nmda = {
    IP.GM: g_NMDA_I,
    IP.VREV: V_E,
    IP.W: 1,
    IP.TAU: tau_NMDA_decay,
    IP.TAU_NMDA_RISE: tau_NMDA_rise,
    IP.ALPHA: alpha,
    IP.BETA: beta,
    IP.GAMMA: gamma,
    IP.MG: Mg2,
}

modelling.fully_connected(MFNonLinearNMDAInput, [pop_e] + pop_e_sel, pop_i, ei_nmda, name='NMDA')

# E->I AMPA
ei_ampa = {
    IP.GM: g_AMPA_rec_E,
    IP.VREV: V_E,
    IP.TAU: tau_AMPA,
}

modelling.fully_connected(MFLinearInput, [pop_e] + pop_e_sel, pop_i, ei_ampa, name='AMPA')

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

modelling.fully_connected(MFLinearInput, pop_i, [pop_e] + pop_e_sel, ie_gaba, name='GABA')

# monitors
N_activity_plot = 15

sp_E_sel = modelling.brian2_spike_monitors(pop_e_sel, n=N_activity_plot)
sp_E = modelling.brian2_spike_monitors(pop_e, n=N_activity_plot)
sp_I = modelling.brian2_spike_monitors(pop_i, n=N_activity_plot)

all_sp = sp_E + sp_I + sp_E_sel

rm_E_sel = modelling.brian2_rate_monitors(pop_e_sel)
rm_E = modelling.brian2_rate_monitors(pop_e)
rm_I = modelling.brian2_rate_monitors(pop_i)

all_rm = rm_E + rm_I + rm_E_sel

# simulate, can be long >120s
system = MFSystem(pop_e, *pop_e_sel, pop_i, name="Brunel Wang 2001")

pop_e.rate = 1 * Hz
pop_i.rate = 6 * Hz
for p in pop_e_sel:
    p.rate = 1 * Hz

solver = MFSolver.rates_voltages(system, solver='simplex', maxiter=1)
#sol = solver.run()
#print(sol)

system.graph().view(cleanup=True)
sdfdsf
# at 1s, select population 1
C_selection = int(f * C_ext)
rate_selection = 50 * Hz
stimuli1 = TimedArray(np.r_[np.zeros(40), np.ones(2), np.zeros(1000)], dt=25 * ms)
input1 = PoissonInput(pop_e_sel[0].brian2, 's_noise_E_sel_0', C_selection, rate_selection, 'stimuli1(t)')

# at 2s, select population 2
stimuli2 = TimedArray(np.r_[np.zeros(80), np.ones(2), np.zeros(100)], dt=25 * ms)
input2 = PoissonInput(pop_e_sel[1].brian2, 's_noise_E_sel_1', C_selection, rate_selection, 'stimuli2(t)')

net = system.collect_brian2_network(*all_sp, *all_rm, input1, input2)
net.run(4 * second, report='stdout')

# plotting

plots.rates(all_rm, 25 * ms)
#plt.xlim([0, 2000])
#plt.ylim([0, 9])
plt.show()

plots.activities(all_sp, N_activity_plot)
plt.show()

