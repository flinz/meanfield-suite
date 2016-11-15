from MFPop import MFpop
from MFSolver import MFSolver, MFSolverRatesVoltages
from MFSource import MFSource
from MFSystem import MFSystem
from Synapse import Synapse
from meanfield_classes import params_standard

# neurons
N = 1000 / 4
N_E = int(N * 0.8)  # pyramidal neurons
N_I = int(N * 0.2)  # interneurons

# voltage
V_L = -70.  # resting
V_thr = -50.
V_reset = -55.
V_E = 0.
V_I = -70.

# membrane capacitance
C_m_E = 0.5
C_m_I = 0.2 # TODO : nano

# membrane leak
g_m_E = 25.
g_m_I = 20.

# refactorty period
tau_rp_E = 2.
tau_rp_I = 1.

# external stimuli
rate = 3
C_ext = 800

# synapses
C_E = N_E
C_I = N_I

# AMPA (excitatory)
g_AMPA_ext_E = 2.08
g_AMPA_rec_E = 0.104 * 800. / N_E
g_AMPA_ext_I = 1.62
g_AMPA_rec_I = 0.081 * 800. / N_E
tau_AMPA = 2.

# NMDA (excitatory)
g_NMDA_E = 0.327 * 800. / N_E
g_NMDA_I = 0.258 * 800. / N_E
tau_NMDA_rise = 2.
tau_NMDA_decay = 100.
alpha = 0. # 0.5 / ms
Mg2 = 0. # 1.

# GABAergic (inhibitory)
g_GABA_E = 1.25 * 200. / N_I
g_GABA_I = 0.973 * 200. / N_I
tau_GABA = 10.

# subpopulations
f = 0.1
p = 1
N_sub = int(N_E * f)
N_non = int(N_E * (1. - f * p))
w_plus = 2.1
w_minus = 1. - f * (w_plus - 1.) / (1. - f)

params_standard = {
    "NMDA": {
        "gamma": 0.280112,
        "beta": 0.062,
    },
    "E": {
        "gamma": 0.280112,
        "beta": 0.062,
        "VE": V_E,
        "E_L": V_L,
        "VI": V_I,
        "V_th": V_thr,
        "V_reset": V_reset,
        "tau_AMPA": tau_AMPA,
        "t_ref": tau_rp_E,
        "C_m": C_m_E * 1e3,
        "g_L": g_m_E,
        "Cext": C_ext,
        "nu_ext": 0.0024,
        "gAMPA": g_AMPA_ext_E,
        "gNMDA": g_NMDA_E,
        "gGABA": g_GABA_E
    },
    "I": {
        "gamma": 0.280112,
        "beta": 0.062,
        "VE": V_E,
        "E_L": V_L,
        "VI": V_I,
        "V_th": V_thr,
        "V_reset": V_reset,
        "tau_AMPA": tau_AMPA,
        "t_ref": tau_rp_I,
        "C_m": C_m_I * 1e3,
        "g_L": g_m_I,
        "Cext": C_ext,
        "nu_ext": 0.0024,
        "gAMPA": g_AMPA_ext_I,
        "gNMDA": g_NMDA_I,
        "gGABA": g_GABA_I
    }
}

initials = {'nu_plus': .025, 'nu_min': .0015, 'nu_i': .009}
# TODO : start values

system = MFSystem("Brunel Wang simplified")

pop_e1 = MFpop("E", params_standard["E"])
pop_e1.n = N_non
pop_e1.rate_ms = initials["nu_plus"]
pop_e1.v_mean = -51. # TODO ?

pop_e2 = MFpop("Edown", params_standard["E"])
pop_e2.n = N_sub
pop_e2.rate_ms = initials["nu_min"]
pop_e2.v_mean = -55.

pop_i = MFpop("I", params_standard["I"])
pop_i.n = N_I
pop_i.rate_ms = initials["nu_i"]
pop_i.v_mean = -55.

system.pops += [pop_e1, pop_e2, pop_i]

# noise pops
source_e_noise1 = MFSource("E_noise1", pop_e1)
source_e_noise1.g_base = params_standard["E"]["gAMPA"]
source_e_noise1.g_dyn = lambda: params_standard["E"]["nu_ext"] * params_standard["E"]["Cext"] * params_standard["E"]["tau_AMPA"]
source_e_noise1.noise_tau = params_standard["E"]["tau_AMPA"]
pop_e1.noise = source_e_noise1

source_e_noise2 = MFSource("E_noise2", pop_e2)
source_e_noise2.g_base = params_standard["E"]["gAMPA"]
source_e_noise2.g_dyn = lambda: params_standard["E"]["nu_ext"] * params_standard["E"]["Cext"] * params_standard["E"]["tau_AMPA"]
source_e_noise2.noise_tau = params_standard["E"]["tau_AMPA"]
pop_e2.noise = source_e_noise2

source_i_noise = MFSource("I_noise", pop_i)
source_i_noise.g_base = params_standard["I"]["gAMPA"]
source_i_noise.g_dyn = lambda: params_standard["I"]["nu_ext"] * params_standard["I"]["Cext"] * params_standard["I"]["tau_AMPA"]
source_i_noise.noise_tau = params_standard["I"]["tau_AMPA"]
pop_i.noise = source_i_noise

# E->E NMDA
syn_spec_nmda = {
    "tau_syn_rise": tau_NMDA_rise,
    "tau_syn_d1": tau_NMDA_decay,
    "tau_syn_d2": tau_NMDA_decay,
    "balance": alpha,
    "tau_x": tau_NMDA_rise,  # depressing
}
syn_ee_nmda = Synapse(**syn_spec_nmda)

source_ee_nmda1 = MFSource('EE Nmda1', pop_e1)
source_ee_nmda1.g_base = params_standard["E"]["gNMDA"]
source_ee_nmda1.g_dyn = lambda: (
        w_plus * syn_ee_nmda(pop_e1.rate_ms) * f + (1. - f) * w_minus * syn_ee_nmda(pop_e2.rate_ms)
    ) * pop_e1.n
#source_ee_nmda1.is_nmda = True
# TODO how to add ampa

source_ee_nmda2 = MFSource('EE Nmda2', pop_e2)
source_ee_nmda2.g_base = params_standard["E"]["gNMDA"]
source_ee_nmda2.g_dyn = lambda: (
    # TODO (f * w_plus + (1. - 2. * f) * w_minus)
        w_minus * syn_ee_nmda(pop_e1.rate_ms) * f + (1. - f) * w_minus * syn_ee_nmda(pop_e2.rate_ms)
    ) * pop_e2.n
source_ee_nmda2.is_nmda = True

# E->I NMDA
syn_ie_nmda = Synapse(**syn_spec_nmda)
source_ie_nmda = MFSource('IE Nmda', pop_i)
source_ie_nmda.g_base = params_standard["I"]["gNMDA"]
source_ie_nmda.g_dyn = lambda: (
        syn_ie_nmda(pop_e1.rate_ms) * f + (1. - f) * syn_ie_nmda(pop_e2.rate_ms)
    ) * pop_e1.n
#source_ie_nmda.is_nmda = True

# I->I GABA
syn_spec_gaba = {
    "tau_syn_rise": 0.,
    "tau_syn_d1": tau_GABA,
    "tau_syn_d2": tau_GABA,
    "balance": 0.,
}
syn_ii_gaba = Synapse(**syn_spec_gaba)
source_ii_gaba = MFSource('II Gaba', pop_i)
source_ii_gaba.g_base = params_standard["I"]["gGABA"]
source_ii_gaba.g_dyn = lambda: syn_ii_gaba(pop_i.rate_ms) * pop_i.n
source_ii_gaba.E_rev = -70.

# I->E GABA
syn_ei_gaba = Synapse(**syn_spec_gaba)
source_ei_gaba1 = MFSource('EI Gaba', pop_e1)
source_ei_gaba1.g_base = params_standard["E"]["gGABA"]
source_ei_gaba1.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
source_ei_gaba1.E_rev = -70. # TODO : adjust

source_ei_gaba2 = MFSource('EI Gaba', pop_e2)
source_ei_gaba2.g_base = params_standard["E"]["gGABA"]
source_ei_gaba2.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
source_ei_gaba2.E_rev = -70.

# TODO : error when using no subpop

solver = MFSolverRatesVoltages(system)
res = solver.run()
