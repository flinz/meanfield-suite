from brian2 import *

from meanfield.MFSystem import MFSystem
from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.inputs.MFNonLinearNMDAInput import MFNonLinearNMDAInput
from meanfield.inputs.MFStaticInput import MFStaticInput
from meanfield.parameters import Connection
from meanfield.parameters import IP
from meanfield.parameters import PP
from meanfield.parameters.MFParameters import MFParameters
from meanfield.populations.MFLinearPopulation import MFLinearPopulation


# BRUNEL, Nicolas et WANG, Xiao-Jing. Effects of neuromodulation in a cortical network model of object working memory dominated by recurrent inhibition. Journal of computational neuroscience, 2001, vol. 11, no 1, p. 63-85.

def one_subpopulation():

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

    pop_e1 = MFLinearPopulation(N_non, E_params, name="E")
    pop_e2 = MFLinearPopulation(N_sub, E_params, name="Edown")
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
    ee_nmda = MFParameters({
        IP.GM: g_NMDA_E,
        IP.VREV: V_E,
        IP.TAU: tau_NMDA_decay,
        IP.TAU_NMDA_RISE: tau_NMDA_rise,
        IP.ALPHA: alpha,
        IP.BETA: beta,
        IP.GAMMA: gamma,
        IP.MG: Mg2,
    })

    MFNonLinearNMDAInput(pop_e1, pop_e1, ee_nmda + {IP.W: 1}, name='EE NMDA 1', connection=Connection.all_to_others())

    MFNonLinearNMDAInput(pop_e1, pop_e2, ee_nmda + {IP.W: w_minus}, name='EE NMDA 2')

    MFNonLinearNMDAInput(pop_e2, pop_e2, ee_nmda + {IP.W: w_plus}, name='EE NMDA 3',
                         connection=Connection.all_to_others())

    MFNonLinearNMDAInput(pop_e2, pop_e1, ee_nmda + {IP.W: 1}, name='EE NMDA 4')

    # E->E AMPA
    ee_ampa = MFParameters({
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
        IP.W: 1,
        IP.TAU: tau_NMDA_decay,
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

    return MFSystem(pop_e1, pop_e2, pop_i, name="Brunel Wang 2001")
