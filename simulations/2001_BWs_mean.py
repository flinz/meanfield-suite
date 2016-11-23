from MFPop import MFpop
from MFSolver import MFSolver, MFSolverRatesVoltages
from MFSource import MFSource
from MFSystem import MFSystem

# neurons
N = 1000
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
rate = 0.003
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
w_plus = 1
w_minus = 1. - f * (w_plus - 1.) / (1. - f)

print(N_non + N_sub)

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
        "nu_ext": rate,
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
        "nu_ext": rate,
        "gAMPA": g_AMPA_ext_I,
        "gNMDA": g_NMDA_I,
        "gGABA": g_GABA_I
    }
}

initials = {'nu_e': 0.01, 'nu_i': 0.01}

system = MFSystem("Brunel Wang simplified")

pop_e1 = MFpop("E", params_standard["E"])
pop_e1.n = N_non
pop_e1.rate_ms = initials["nu_e"]
pop_e1.v_mean = -55.

pop_e2 = MFpop("Edown", params_standard["E"])
pop_e2.n = N_sub
pop_e2.rate_ms = initials["nu_e"]
pop_e2.v_mean = -55.

pop_i = MFpop("I", params_standard["I"])
pop_i.n = N_I
pop_i.rate_ms = initials["nu_i"]
pop_i.v_mean = -55.

system.pops += [pop_e1, pop_e2, pop_i]

# noise pops
source_e_noise1 = MFSource("E_noise1", pop_e1)
source_e_noise1.noise_tau = params_standard["E"]["tau_AMPA"]
source_e_noise1.g_base = params_standard["E"]["gAMPA"]
source_e_noise1.g_dyn = lambda: params_standard["E"]["nu_ext"] * params_standard["E"]["Cext"] * params_standard["E"]["tau_AMPA"]
pop_e1.noise = source_e_noise1

source_e_noise2 = MFSource("E_noise2", pop_e2)
source_e_noise2.noise_tau = params_standard["E"]["tau_AMPA"]
source_e_noise2.g_base = params_standard["E"]["gAMPA"]
source_e_noise2.g_dyn = lambda: params_standard["E"]["nu_ext"] * params_standard["E"]["Cext"] * params_standard["E"]["tau_AMPA"]
pop_e2.noise = source_e_noise2

source_i_noise = MFSource("I_noise", pop_i)
source_i_noise.noise_tau = params_standard["I"]["tau_AMPA"]
source_i_noise.g_base = params_standard["I"]["gAMPA"]
source_i_noise.g_dyn = lambda: params_standard["I"]["nu_ext"] * params_standard["I"]["Cext"] * params_standard["I"]["tau_AMPA"]
pop_i.noise = source_i_noise

# E->E NMDA
source_ee_nmda1 = MFSource('EE NMDA 1', pop_e1)
source_ee_nmda1.is_nmda = False
source_ee_nmda1.g_base = params_standard["E"]["gNMDA"]
source_ee_nmda1.g_dyn = lambda: (
    pop_e1.n * f * w_plus * pop_e1.rate_ms * tau_NMDA_decay +
    pop_e2.n * (1. - f) * w_minus * pop_e2.rate_ms * tau_NMDA_decay
)

source_ee_nmda2 = MFSource('EE NMDA 2', pop_e2)
source_ee_nmda2.is_nmda = False
source_ee_nmda2.g_base = params_standard["E"]["gNMDA"]
source_ee_nmda2.g_dyn = lambda: (
    pop_e1.n * f * w_minus * pop_e1.rate_ms * tau_NMDA_decay +
    pop_e2.n * (1. - f) * w_plus * pop_e2.rate_ms * tau_NMDA_decay
)

# E->E AMPA
source_ee_ampa1 = MFSource('EE AMPA 1', pop_e1)
source_ee_ampa1.is_nmda = False
source_ee_ampa1.g_base = params_standard["E"]["gAMPA"]
source_ee_ampa1.g_dyn = lambda: (
    pop_e1.n * f * w_plus * pop_e1.rate_ms * tau_AMPA +
    pop_e2.n * (1. - f) * w_minus * pop_e2.rate_ms * tau_AMPA
)

source_ee_ampa2 = MFSource('EE AMPA 2', pop_e2)
source_ee_ampa2.is_nmda = False
source_ee_ampa2.g_base = params_standard["E"]["gAMPA"]
source_ee_ampa2.g_dyn = lambda: (
    pop_e1.n * f * w_minus * pop_e1.rate_ms * tau_AMPA +
    pop_e2.n * (1. - f) * w_plus * pop_e2.rate_ms * tau_AMPA
)

# E->I NMDA
source_ie_nmda = MFSource('IE NMDA', pop_i)
source_ie_nmda.is_nmda = False
source_ie_nmda.g_base = params_standard["I"]["gNMDA"]
source_ie_nmda.g_dyn = lambda: (
    pop_e1.n * f * pop_e1.rate_ms * tau_NMDA_decay +
    pop_e2.n * (1. - f) * pop_e2.rate_ms * tau_NMDA_decay
)

# E->I AMPA
source_ie_ampa = MFSource('IE AMPA', pop_i)
source_ie_ampa.is_nmda = False
source_ie_ampa.g_base = params_standard["I"]["gAMPA"]
source_ie_ampa.g_dyn = lambda: (
    pop_e1.n * f * pop_e1.rate_ms * tau_AMPA +
    pop_e2.n * (1. - f) * pop_e2.rate_ms * tau_AMPA
)

# I->I GABA
source_ii_gaba = MFSource('II GABA', pop_i)
source_ii_gaba.E_rev = -70.
source_ii_gaba.g_base = params_standard["I"]["gGABA"]
source_ii_gaba.g_dyn = lambda: pop_i.n * pop_i.rate_ms * tau_GABA

# I->E GABA
source_ei_gaba1 = MFSource('EI Gaba', pop_e1)
source_ei_gaba1.E_rev = -70. # TODO : adjust
source_ei_gaba1.g_base = params_standard["E"]["gGABA"]
source_ei_gaba1.g_dyn = lambda: pop_i.n * pop_i.rate_ms * tau_GABA

source_ei_gaba2 = MFSource('EI Gaba', pop_e2)
source_ei_gaba2.E_rev = -70.
source_ei_gaba2.g_base = params_standard["E"]["gGABA"]
source_ei_gaba2.g_dyn = lambda: pop_i.n * pop_i.rate_ms * tau_GABA


# TODO : error when using no subpop

solver = MFSolverRatesVoltages(system, maxiter=100)#, solver="gradient")
print(solver.mfstate.state)

res = solver.run()
