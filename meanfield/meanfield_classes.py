from brian2.units import *

from meanfield.populations.MFLinearPop import MFLinearPop
from meanfield.sources.MFStaticSource import MFStaticSource
from meanfield.sources.MFNonLinearNMDASource import MFNonLinearNMDASource
from meanfield.MFSystem import MFSystem
from meanfield.parameters import NP, SP

params_standard = {
    'NMDA': {
        'gamma': 0.280112,
        'beta': 0.062,
    },
    'E': {
        'gamma': 0.280112,
        'beta': 0.062,
        'VE': 0.,
        'V_L': -70. * mV,
        'V_th': -50. * mV,
        'V_reset': -60. * mV,
        'tau_AMPA': 2. * ms,
        't_ref': 2. * ms,
        'C_m': 500. * pF,
        'g_L': 25. * nsiemens,
        'Cext': 1000,
        'nu_ext': 2.4 * Hz,
        'gAMPA': 2.08 * nsiemens,
        'gNMDA': 0.327 * nsiemens,
        'gGABA': 1.25 * nsiemens
    },
    'I': {
        'gamma': 0.280112,
        'beta': 0.062,
        'VE': 0.,
        'V_L': -70. * mV,
        'V_th': -50. * mV,
        'V_reset': -60. * mV,
        'tau_AMPA': 2. * ms,
        't_ref': 1. * ms,
        'C_m': 200. * pF,
        'g_L': 20. * nsiemens,
        'Cext': 1000,
        'nu_ext': 2.4 * Hz,
        'gAMPA': 1.62 * nsiemens,
        'gNMDA': 0.258 * nsiemens,
        'gGABA': 0.973 * nsiemens
    }
}


def setup_brunel99_nmda():

    # brunel 1999 system, one up state pop
    ff = .1
    w_plus = 2.5
    w_min = 1. - ff * (w_plus - 1.) / (1. - ff)
    initials = {'nu_plus': 25 * Hz, 'nu_min': 1.5 * Hz, 'nu_i': 9 * Hz}

    system = MFSystem('Brunel')
    pop_e1 = MFNonLinearNMDASource('Eup', 800 * ff , {
        NP.GM: params_standard['E']['g_L'],
        NP.VL: params_standard['E']['V_L'],
        NP.CM: params_standard['E']['C_m'],
        NP.VTHR: params_standard['E']['V_th'],
        NP.VRES: params_standard['E']['V_reset'],
        NP.TAU_RP: params_standard['E']['tau_AMPA']
    })
    pop_e2 = MFNonLinearPop('Edown', 800 * (1. - ff), {
        NP.GM: params_standard['E']['g_L'],
        NP.VL: params_standard['E']['V_L'],
        NP.CM: params_standard['E']['C_m'],
        NP.VTHR: params_standard['E']['V_th'],
        NP.VRES: params_standard['E']['V_reset'],
        NP.TAU_RP: params_standard['E']['tau_AMPA']
    })
    pop_i = MFNonLinearPop('I', 200, {
        NP.GM: params_standard['I']['g_L'],
        NP.VL: params_standard['I']['V_L'],
        NP.CM: params_standard['I']['C_m'],
        NP.VTHR: params_standard['I']['V_th'],
        NP.VRES: params_standard['I']['V_reset'],
        NP.TAU_RP: params_standard['I']['tau_AMPA']
    })
    pop_e1.rate = initials['nu_plus']
    pop_e1.v_mean = -51. * mV
    pop_e2.rate = initials['nu_min']
    pop_e2.v_mean = -55. * mV
    pop_i.rate = initials['nu_i']
    pop_i.v_mean = -55. * mV

    system.pops += [pop_e1, pop_e2, pop_i]

    source_e_noise1 = MFStaticSource('E_noise1', pop_e1, {
        SP.GM: params_standard['E']['gAMPA'],
        SP.VE: 0 * mV,
        SP.TAU: params_standard['E']['tau_AMPA'],
    }, params_standard['E']['nu_ext'], params_standard['E']['Cext'])
    #source_e_noise1.g_base = params_standard['E']['gAMPA']
    #source_e_noise1.g_dyn = lambda: params_standard['E']['nu_ext'] * params_standard['E']['Cext'] * params_standard['E']['tau_AMPA']
    #source_e_noise1.noise_tau = params_standard['E']['tau_AMPA']
    pop_e1.noise = source_e_noise1

    source_e_noise2 = MFStaticSource('E_noise2', pop_e2, {
        SP.GM: params_standard['E']['gAMPA'],
        SP.VE: 0 * mV,
        SP.TAU: params_standard['E']['tau_AMPA'],
    }, params_standard['E']['nu_ext'], params_standard['E']['Cext'])
    #source_e_noise2.g_base = params_standard['E']['gAMPA']
    #source_e_noise2.g_dyn = lambda: params_standard['E']['nu_ext'] * params_standard['E']['Cext'] * params_standard['E']['tau_AMPA']
    #source_e_noise2.noise_tau = params_standard['E']['tau_AMPA']
    pop_e2.noise = source_e_noise2

    source_i_noise = MFStaticSource('I_noise', pop_i, {
        SP.GM: params_standard['I']['gAMPA'],
        SP.VE: 0 * mV,
        SP.TAU: params_standard['I']['tau_AMPA'],
    }, params_standard['I']['nu_ext'],  params_standard['I']['Cext'])
    #source_i_noise.g_base = params_standard['I']['gAMPA']
    #source_i_noise.g_dyn = lambda: params_standard['I']['nu_ext'] * params_standard['I']['Cext'] * params_standard['I']['tau_AMPA']
    #source_i_noise.noise_tau = params_standard['I']['tau_AMPA']
    pop_i.noise = source_i_noise

    # E->E NMDA
    syn_spec_nmda = {
        'tau_syn_rise': 1. * ms,
        'tau_syn_d1': 100. * ms,
        'tau_syn_d2': 100. * ms,
        'balance': .5,
        'tau_x': 150. * ms,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda1 = MFDynamicSource('EE Nmda1', pop_e1, {
        SP.GM: params_standard['E']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: w_plus
    }, from_pop=pop_e1, synapse=syn_ee_nmda)
    source_ee_nmda12 = MFDynamicSource('EE Nmda12', pop_e1, {
        SP.GM: params_standard['E']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: w_min
    }, from_pop=pop_e2, synapse=syn_ee_nmda)
    #source_ee_nmda1.g_base = params_standard['E']['gNMDA']
    #source_ee_nmda1.g_dyn = lambda: (
    #        ff * w_plus * syn_ee_nmda(pop_e1.rate_ms) + (1. - ff) * w_min * syn_ee_nmda(pop_e2.rate_ms)
    #    ) * pop_e1.n
    #source_ee_nmda1.is_nmda = True

    source_ee_nmda2 = MFDynamicSource('EE Nmda2', pop_e2, {
        SP.GM: params_standard['E']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: w_min
    }, from_pop=pop_e1, synapse=syn_ee_nmda)
    source_ee_nmda22 = MFDynamicSource('EE Nmda22', pop_e2, {
        SP.GM: params_standard['E']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: (w_plus * ff + (1. - 2. * ff) * w_min) / (1 - ff)
    }, from_pop=pop_e2, synapse=syn_ee_nmda)
    #source_ee_nmda2.g_base = params_standard['E']['gNMDA']
    #source_ee_nmda2.g_dyn = lambda: (
    #        ff * w_min * syn_ee_nmda(pop_e1.rate_ms) + (ff * w_plus + (1. - 2.*ff) * w_min) * syn_ee_nmda(pop_e2.rate_ms)
    #    ) * pop_e2.n
    #source_ee_nmda2.is_nmda = True

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFDynamicSource('IE Nmda', pop_i, {
        SP.GM: params_standard['I']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: 1
    }, from_pop=pop_e1, synapse=syn_ie_nmda)
    source_ie_nmda2 = MFDynamicSource('IE Nmda2', pop_i, {
        SP.GM: params_standard['I']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: 1
    }, from_pop=pop_e2, synapse=syn_ie_nmda)
    #source_ie_nmda.g_base = params_standard['I']['gNMDA']
    #source_ie_nmda.g_dyn = lambda: (
    #        ff * syn_ie_nmda(pop_e1.rate_ms) + (1. - ff) * syn_ie_nmda(pop_e2.rate_ms)
    #    ) * pop_e1.n
    #source_ie_nmda.is_nmda = True

    # I->I GABA
    syn_spec_gaba = {
        'tau_syn_rise': 0. * ms,
        'tau_syn_d1': 10. * ms,
        'tau_syn_d2': 10. * ms,
        'balance': .5,
    }
    syn_ii_gaba = Synapse(**syn_spec_gaba)

    source_ii_gaba = MFDynamicSource('II Gaba', pop_i, {
        SP.GM: params_standard['I']['gGABA'],
        SP.VE: -70 * mV,
        SP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, from_pop=pop_i, synapse=syn_ii_gaba)
    #source_ii_gaba.g_base = params_standard['I']['gGABA']
    #source_ii_gaba.g_dyn = lambda: syn_ii_gaba(pop_i.rate_ms) * pop_i.n
    #source_ii_gaba.E_rev = -70.

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba1 = MFDynamicSource('EI Gaba', pop_e1, {
        SP.GM: params_standard['E']['gGABA'],
        SP.VE: -70 * mV,
        SP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, from_pop=pop_i, synapse=syn_ei_gaba)
    #source_ei_gaba1.g_base = params_standard['E']['gGABA']
    #source_ei_gaba1.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    #source_ei_gaba1.E_rev = -70.

    source_ei_gaba2 = MFDynamicSource('EI Gaba', pop_e2, {
        SP.GM: params_standard['E']['gGABA'],
        SP.VE: -70 * mV,
        SP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, from_pop=pop_i, synapse=syn_ei_gaba)
    #source_ei_gaba2.g_base = params_standard['E']['gGABA']
    #source_ei_gaba2.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    #source_ei_gaba2.E_rev = -70.

    return system




def setup_brunel99(w_plus_val=2.5):

    # brunel 1999 system, one up state pop
    ff = .1
    w_plus = w_plus_val
    w_min = 1. - ff * (w_plus - 1.) / (1. - ff)
    initials = {'nu_plus': 25 * Hz, 'nu_min': 1.5 * Hz, 'nu_i': 9 * Hz}

    system = MFSystem('Brunel')
    pop_e1 = MFLinearPop('Eup', 800 * ff , {
        NP.GM: params_standard['E']['g_L'],
        NP.VL: params_standard['E']['V_L'],
        NP.CM: params_standard['E']['C_m'],
        NP.VTHR: params_standard['E']['V_th'],
        NP.VRES: params_standard['E']['V_reset'],
        NP.TAU_RP: params_standard['E']['tau_AMPA']
    })
    #pop_e1.n = 800
    pop_e2 = MFLinearPop('Edown', 800 * (1. - ff), {
        NP.GM: params_standard['E']['g_L'],
        NP.VL: params_standard['E']['V_L'],
        NP.CM: params_standard['E']['C_m'],
        NP.VTHR: params_standard['E']['V_th'],
        NP.VRES: params_standard['E']['V_reset'],
        NP.TAU_RP: params_standard['E']['tau_AMPA']
    })
    #pop_e2.n = 800
    pop_i = MFLinearPop('I', 200, {
        NP.GM: params_standard['I']['g_L'],
        NP.VL: params_standard['I']['V_L'],
        NP.CM: params_standard['I']['C_m'],
        NP.VTHR: params_standard['I']['V_th'],
        NP.VRES: params_standard['I']['V_reset'],
        NP.TAU_RP: params_standard['I']['tau_AMPA']
    })
    #pop_i.n = 200
    pop_e1.rate = initials['nu_plus']
    pop_e1.v_mean = -51. * mV
    pop_e2.rate = initials['nu_min']
    pop_e2.v_mean = -55. * mV
    pop_i.rate = initials['nu_i']
    pop_i.v_mean = -55. * mV

    system.pops += [pop_e1, pop_e2, pop_i]

    # noise pops

    source_e_noise1 = MFStaticSource('E_noise1', pop_e1, {
        SP.GM: params_standard['E']['gAMPA'],
        SP.VE: 0 * mV,
        SP.TAU: params_standard['E']['tau_AMPA'],
    }, params_standard['E']['nu_ext'], params_standard['E']['Cext'])
    #source_e_noise1.g_base = params_standard['E']['gAMPA']
    #source_e_noise1.g_dyn = lambda: params_standard['E']['nu_ext'] * params_standard['E']['Cext'] * params_standard['E']['tau_AMPA']
    #source_e_noise1.noise_tau = params_standard['E']['tau_AMPA']
    pop_e1.noise = source_e_noise1

    source_e_noise2 = MFStaticSource('E_noise2', pop_e2, {
        SP.GM: params_standard['E']['gAMPA'],
        SP.VE: 0 * mV,
        SP.TAU: params_standard['E']['tau_AMPA'],
    }, params_standard['E']['nu_ext'], params_standard['E']['Cext'])
    #source_e_noise2.g_base = params_standard['E']['gAMPA']
    #source_e_noise2.g_dyn = lambda: params_standard['E']['nu_ext'] * params_standard['E']['Cext'] * params_standard['E']['tau_AMPA']
    #source_e_noise2.noise_tau = params_standard['E']['tau_AMPA']
    pop_e2.noise = source_e_noise2

    source_i_noise = MFStaticSource('I_noise', pop_i, {
        SP.GM: params_standard['I']['gAMPA'],
        SP.VE: 0 * mV,
        SP.TAU: params_standard['I']['tau_AMPA'],
    }, params_standard['I']['nu_ext'],  params_standard['I']['Cext'])
    #source_i_noise.g_base = params_standard['I']['gAMPA']
    #source_i_noise.g_dyn = lambda: params_standard['I']['nu_ext'] * params_standard['I']['Cext'] * params_standard['I']['tau_AMPA']
    #source_i_noise.noise_tau = params_standard['I']['tau_AMPA']
    pop_i.noise = source_i_noise

    # E->E NMDA
    syn_spec_nmda = {
        'tau_syn_rise': 1. * ms,
        'tau_syn_d1': 100. * ms,
        'tau_syn_d2': 100. * ms,
        'balance': .5,
        'tau_x': 150. * ms,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda1 = MFDynamicSource('EE Nmda1', pop_e1, {
        SP.GM: params_standard['E']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: w_plus
    }, from_pop=pop_e1, synapse=syn_ee_nmda)
    source_ee_nmda12 = MFDynamicSource('EE Nmda12', pop_e1, {
        SP.GM: params_standard['E']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: w_min
    }, from_pop=pop_e2, synapse=syn_ee_nmda)
    #source_ee_nmda1.g_base = params_standard['E']['gNMDA']
    #source_ee_nmda1.g_dyn = lambda: (
    #        ff * w_plus * syn_ee_nmda(pop_e1.rate_ms) + (1. - ff) * w_min * syn_ee_nmda(pop_e2.rate_ms)
    #    ) * pop_e1.n
    #source_ee_nmda1.is_nmda = True

    source_ee_nmda2 = MFDynamicSource('EE Nmda2', pop_e2, {
        SP.GM: params_standard['E']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: w_min
    }, from_pop=pop_e1, synapse=syn_ee_nmda)
    source_ee_nmda22 = MFDynamicSource('EE Nmda22', pop_e2, {
        SP.GM: params_standard['E']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: (w_plus * ff + (1. - 2. * ff) * w_min) / (1 - ff)
    }, from_pop=pop_e2, synapse=syn_ee_nmda)
    #source_ee_nmda2.g_base = params_standard['E']['gNMDA']
    #source_ee_nmda2.g_dyn = lambda: (
    #        ff * w_min * syn_ee_nmda(pop_e1.rate_ms) + (ff * w_plus + (1. - 2.*ff) * w_min) * syn_ee_nmda(pop_e2.rate_ms)
    #    ) * pop_e2.n
    #source_ee_nmda2.is_nmda = True

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFDynamicSource('IE Nmda', pop_i, {
        SP.GM: params_standard['I']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: 1
    }, from_pop=pop_e1, synapse=syn_ie_nmda)
    source_ie_nmda2 = MFDynamicSource('IE Nmda2', pop_i, {
        SP.GM: params_standard['I']['gNMDA'],
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],
        SP.W: 1
    }, from_pop=pop_e2, synapse=syn_ie_nmda)
    #source_ie_nmda.g_base = params_standard['I']['gNMDA']
    #source_ie_nmda.g_dyn = lambda: (
    #        ff * syn_ie_nmda(pop_e1.rate_ms) + (1. - ff) * syn_ie_nmda(pop_e2.rate_ms)
    #    ) * pop_e1.n
    #source_ie_nmda.is_nmda = True

    # I->I GABA
    syn_spec_gaba = {
        'tau_syn_rise': 0. * ms,
        'tau_syn_d1': 10. * ms,
        'tau_syn_d2': 10. * ms,
        'balance': .5,
    }
    syn_ii_gaba = Synapse(**syn_spec_gaba)

    source_ii_gaba = MFDynamicSource('II Gaba', pop_i, {
        SP.GM: params_standard['I']['gGABA'],
        SP.VE: -70 * mV,
        SP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, from_pop=pop_i, synapse=syn_ii_gaba)
    #source_ii_gaba.g_base = params_standard['I']['gGABA']
    #source_ii_gaba.g_dyn = lambda: syn_ii_gaba(pop_i.rate_ms) * pop_i.n
    #source_ii_gaba.E_rev = -70.

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba1 = MFDynamicSource('EI Gaba', pop_e1, {
        SP.GM: params_standard['E']['gGABA'],
        SP.VE: -70 * mV,
        SP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, from_pop=pop_i, synapse=syn_ei_gaba)
    #source_ei_gaba1.g_base = params_standard['E']['gGABA']
    #source_ei_gaba1.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    #source_ei_gaba1.E_rev = -70.

    source_ei_gaba2 = MFDynamicSource('EI Gaba', pop_e2, {
        SP.GM: params_standard['E']['gGABA'],
        SP.VE: -70 * mV,
        SP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, from_pop=pop_i, synapse=syn_ei_gaba)
    #source_ei_gaba2.g_base = params_standard['E']['gGABA']
    #source_ei_gaba2.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    #source_ei_gaba2.E_rev = -70.

    return system


def setup_EI(has_nmda=False):

    # brunel 1999 system, one up state pop
    initials = {'nu_e': 3 * Hz, 'nu_i': 9 * Hz}

    mult = 0.25
    if has_nmda:
        mult = 1.

    system = MFSystem('EI')
    pop_e = MFLinearPop('E', 800, {
        NP.GM: params_standard['E']['g_L'],
        NP.VL: params_standard['E']['V_L'],
        NP.CM: params_standard['E']['C_m'],
        NP.VTHR: params_standard['E']['V_th'],
        NP.VRES: params_standard['E']['V_reset'],
        NP.TAU_RP: params_standard['E']['tau_AMPA']
    })
    #pop_e.n = 800
    pop_i = MFLinearPop('I', 200, {
        NP.GM: params_standard['I']['g_L'],
        NP.VL: params_standard['I']['V_L'],
        NP.CM: params_standard['I']['C_m'],
        NP.VTHR: params_standard['I']['V_th'],
        NP.VRES: params_standard['I']['V_reset'],
        NP.TAU_RP: params_standard['I']['tau_AMPA']
    })
    #pop_i.n = 200
    pop_e.rate = initials['nu_e']
    pop_e.v_mean = -55. * mV
    pop_i.rate = initials['nu_i']
    pop_i.v_mean = -55. * mV

    system.pops += [pop_e, pop_i]

    # noise pops
    source_e_noise = MFStaticSource('E_noise', pop_e, {
        SP.GM: params_standard['E']['gAMPA'],
        SP.VE: 0 * mV,
        SP.TAU: params_standard['E']['tau_AMPA'],
    }, params_standard['E']['nu_ext'], params_standard['E']['Cext'])
    #source_e_noise.g_base = params_standard['E']['gAMPA']
    #source_e_noise.g_dyn = lambda: params_standard['E']['nu_ext'] * params_standard['E']['Cext'] * params_standard['E']['tau_AMPA']
    #source_e_noise.noise_tau = params_standard['E']['tau_AMPA']
    pop_e.noise = source_e_noise

    source_i_noise = MFStaticSource('I_noise', pop_i, {
        SP.GM: params_standard['I']['gAMPA'],
        SP.VE: 0 * mV,
        SP.TAU: params_standard['I']['tau_AMPA'],
    }, params_standard['I']['nu_ext'], params_standard['I']['Cext'])
    #source_i_noise.g_base = params_standard['I']['gAMPA']
    #source_i_noise.g_dyn = lambda: params_standard['I']['nu_ext'] * params_standard['I']['Cext'] * params_standard['I']['tau_AMPA']
    #source_i_noise.noise_tau = params_standard['I']['tau_AMPA']
    pop_i.noise = source_i_noise

    # E->E NMDA
    syn_spec_nmda = {
        'tau_syn_rise': 1. * ms,
        'tau_syn_d1': 100. * ms,
        'tau_syn_d2': 100. * ms,
        'balance': .5,
        'tau_x': 150. * ms,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda = MFDynamicSource('EE Nmda', pop_e, {
        SP.GM: params_standard['E']['gNMDA'] * mult,
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],  # not sure
    }, from_pop=pop_e, synapse=syn_ee_nmda)
    #source_ee_nmda.g_base = params_standard['E']['gNMDA'] * mult
    #source_ee_nmda.g_dyn = lambda: syn_ee_nmda(pop_e.rate_ms) * pop_e.n
    #source_ee_nmda.is_nmda = has_nmda

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFDynamicSource('IE Nmda', pop_i, {
        SP.GM: params_standard['I']['gNMDA'] * mult,
        SP.VE: 0 * mV,
        SP.TAU: syn_spec_nmda['tau_syn_d1'],  # not sure
    }, from_pop=pop_e, synapse=syn_ie_nmda)
    #source_ie_nmda.g_base = params_standard['I']['gNMDA'] * mult
    #source_ie_nmda.g_dyn = lambda: syn_ie_nmda(pop_e.rate_ms) * pop_e.n
    #source_ie_nmda.is_nmda = has_nmda

    # I->I GABA
    syn_spec_gaba = {
        'tau_syn_rise': 0. * ms,
        'tau_syn_d1': 10. * ms,
        'tau_syn_d2': 10. * ms,
        'balance': .5,
    }
    syn_ii_gaba = Synapse(**syn_spec_gaba)
    source_ii_gaba = MFDynamicSource('II Gaba', pop_i, {
        SP.GM: params_standard['I']['gGABA'],
        SP.VE: -70 * mV,
        SP.TAU: syn_spec_gaba['tau_syn_d1'],  # not sure
    }, from_pop=pop_i, synapse=syn_ii_gaba)
    #source_ii_gaba.g_base = params_standard['I']['gGABA']
    #source_ii_gaba.g_dyn = lambda: syn_ii_gaba(pop_i.rate_ms) * pop_i.n
    #source_ii_gaba.E_rev = -70.

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba = MFDynamicSource('EI Gaba', pop_e, {
        SP.GM: params_standard['E']['gGABA'],
        SP.VE: -70 * mV,
        SP.TAU: syn_spec_gaba['tau_syn_d1'],  # not sure
    }, from_pop=pop_i, synapse=syn_ei_gaba)
    #source_ei_gaba.g_base = params_standard['E']['gGABA']
    #source_ei_gaba.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    #source_ei_gaba.E_rev = -70.

    return system
