from brian2.units import *

from meanfield.MFSystem import MFSystem
from meanfield.parameters import PP, IP
from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.inputs.MFStaticInput import MFStaticInput
from meanfield.inputs.Synapses import Synapse

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
    pop_e1 = MFLinearPopulation('Eup', 800 * ff, {
        PP.GM: params_standard['E']['g_L'],
        PP.VL: params_standard['E']['V_L'],
        PP.CM: params_standard['E']['C_m'],
        PP.VTHR: params_standard['E']['V_th'],
        PP.VRES: params_standard['E']['V_reset'],
        PP.TAU_RP: params_standard['E']['tau_AMPA']
    })
    pop_e2 = MFLinearPopulation('Edown', 800 * (1. - ff), {
        PP.GM: params_standard['E']['g_L'],
        PP.VL: params_standard['E']['V_L'],
        PP.CM: params_standard['E']['C_m'],
        PP.VTHR: params_standard['E']['V_th'],
        PP.VRES: params_standard['E']['V_reset'],
        PP.TAU_RP: params_standard['E']['tau_AMPA']
    })
    pop_i = MFLinearPopulation('I', 200, {
        PP.GM: params_standard['I']['g_L'],
        PP.VL: params_standard['I']['V_L'],
        PP.CM: params_standard['I']['C_m'],
        PP.VTHR: params_standard['I']['V_th'],
        PP.VRES: params_standard['I']['V_reset'],
        PP.TAU_RP: params_standard['I']['tau_AMPA']
    })
    pop_e1.rate = initials['nu_plus']
    pop_e1.v_mean = -51. * mV
    pop_e2.rate = initials['nu_min']
    pop_e2.v_mean = -55. * mV
    pop_i.rate = initials['nu_i']
    pop_i.v_mean = -55. * mV

    system.populations += [pop_e1, pop_e2, pop_i]

    source_e_noise1 = MFStaticInput('E_noise1', pop_e1, params_standard['E']['Cext'], params_standard['E']['nu_ext'], {
        IP.GM: params_standard['E']['gAMPA'],
        IP.VREV: 0 * mV,
        IP.TAU: params_standard['E']['tau_AMPA'],
    })
    #source_e_noise1.g_base = params_standard['E']['gAMPA']
    #source_e_noise1.g_dyn = lambda: params_standard['E']['nu_ext'] * params_standard['E']['Cext'] * params_standard['E']['tau_AMPA']
    #source_e_noise1.noise_tau = params_standard['E']['tau_AMPA']
    pop_e1.noise = source_e_noise1

    source_e_noise2 = MFStaticInput('E_noise2', pop_e2, params_standard['E']['Cext'], params_standard['E']['nu_ext'], {
        IP.GM: params_standard['E']['gAMPA'],
        IP.VREV: 0 * mV,
        IP.TAU: params_standard['E']['tau_AMPA'],
    })
    #source_e_noise2.g_base = params_standard['E']['gAMPA']
    #source_e_noise2.g_dyn = lambda: params_standard['E']['nu_ext'] * params_standard['E']['Cext'] * params_standard['E']['tau_AMPA']
    #source_e_noise2.noise_tau = params_standard['E']['tau_AMPA']
    pop_e2.noise = source_e_noise2

    source_i_noise = MFStaticInput('I_noise', pop_i, params_standard['I']['Cext'], params_standard['I']['nu_ext'], {
        IP.GM: params_standard['I']['gAMPA'],
        IP.VREV: 0 * mV,
        IP.TAU: params_standard['I']['tau_AMPA'],
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

    source_ee_nmda1 = MFLinearInput('EE Nmda1', pop_e1, {
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: w_plus
    }, origin=pop_e1)#, synapse=syn_ee_nmda)
    source_ee_nmda12 = MFLinearInput('EE Nmda12', pop_e1, {
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: w_min
    }, origin=pop_e2)#, synapse=syn_ee_nmda)
    #source_ee_nmda1.g_base = params_standard['E']['gNMDA']
    #source_ee_nmda1.g_dyn = lambda: (
    #        ff * w_plus * syn_ee_nmda(pop_e1.rate_ms) + (1. - ff) * w_min * syn_ee_nmda(pop_e2.rate_ms)
    #    ) * pop_e1.n
    #source_ee_nmda1.is_nmda = True

    source_ee_nmda2 = MFLinearInput('EE Nmda2', pop_e2, {
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: w_min
    }, origin=pop_e1)#, synapse=syn_ee_nmda)
    source_ee_nmda22 = MFLinearInput('EE Nmda22', pop_e2, {
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: (w_plus * ff + (1. - 2. * ff) * w_min) / (1 - ff)
    }, origin=pop_e2)#, synapse=syn_ee_nmda)
    #source_ee_nmda2.g_base = params_standard['E']['gNMDA']
    #source_ee_nmda2.g_dyn = lambda: (
    #        ff * w_min * syn_ee_nmda(pop_e1.rate_ms) + (ff * w_plus + (1. - 2.*ff) * w_min) * syn_ee_nmda(pop_e2.rate_ms)
    #    ) * pop_e2.n
    #source_ee_nmda2.is_nmda = True

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFLinearInput('IE Nmda', pop_i, {
        IP.GM: params_standard['I']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: 1
    }, origin=pop_e1)#, synapse=syn_ie_nmda)
    source_ie_nmda2 = MFLinearInput('IE Nmda2', pop_i, {
        IP.GM: params_standard['I']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: 1
    }, origin=pop_e2)#, synapse=syn_ie_nmda)
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

    source_ii_gaba = MFLinearInput('II Gaba', pop_i, {
        IP.GM: params_standard['I']['gGABA'],
        IP.VREV: -70 * mV,
        IP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, origin=pop_i)#, synapse=syn_ii_gaba)
    #source_ii_gaba.g_base = params_standard['I']['gGABA']
    #source_ii_gaba.g_dyn = lambda: syn_ii_gaba(pop_i.rate_ms) * pop_i.n
    #source_ii_gaba.E_rev = -70.

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba1 = MFLinearInput('EI Gaba', pop_e1, {
        IP.GM: params_standard['E']['gGABA'],
        IP.VREV: -70 * mV,
        IP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, origin=pop_i)#, synapse=syn_ei_gaba)
    #source_ei_gaba1.g_base = params_standard['E']['gGABA']
    #source_ei_gaba1.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    #source_ei_gaba1.E_rev = -70.

    source_ei_gaba2 = MFLinearInput('EI Gaba', pop_e2, {
        IP.GM: params_standard['E']['gGABA'],
        IP.VREV: -70 * mV,
        IP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, origin=pop_i)#, synapse=syn_ei_gaba)
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
    pop_e1 = MFLinearPopulation('Eup', 800 * ff, {
        PP.GM: params_standard['E']['g_L'],
        PP.VL: params_standard['E']['V_L'],
        PP.CM: params_standard['E']['C_m'],
        PP.VTHR: params_standard['E']['V_th'],
        PP.VRES: params_standard['E']['V_reset'],
        PP.TAU_RP: params_standard['E']['tau_AMPA']
    })
    #pop_e1.n = 800
    pop_e2 = MFLinearPopulation('Edown', 800 * (1. - ff), {
        PP.GM: params_standard['E']['g_L'],
        PP.VL: params_standard['E']['V_L'],
        PP.CM: params_standard['E']['C_m'],
        PP.VTHR: params_standard['E']['V_th'],
        PP.VRES: params_standard['E']['V_reset'],
        PP.TAU_RP: params_standard['E']['tau_AMPA']
    })
    #pop_e2.n = 800
    pop_i = MFLinearPopulation('I', 200, {
        PP.GM: params_standard['I']['g_L'],
        PP.VL: params_standard['I']['V_L'],
        PP.CM: params_standard['I']['C_m'],
        PP.VTHR: params_standard['I']['V_th'],
        PP.VRES: params_standard['I']['V_reset'],
        PP.TAU_RP: params_standard['I']['tau_AMPA']
    })
    #pop_i.n = 200
    pop_e1.rate = initials['nu_plus']
    pop_e1.v_mean = -51. * mV
    pop_e2.rate = initials['nu_min']
    pop_e2.v_mean = -55. * mV
    pop_i.rate = initials['nu_i']
    pop_i.v_mean = -55. * mV

    system.populations += [pop_e1, pop_e2, pop_i]

    # noise pops

    source_e_noise1 = MFStaticInput('E_noise1', pop_e1, params_standard['E']['Cext'], params_standard['E']['nu_ext'], {
        IP.GM: params_standard['E']['gAMPA'],
        IP.VREV: 0 * mV,
        IP.TAU: params_standard['E']['tau_AMPA'],
    })
    #source_e_noise1.g_base = params_standard['E']['gAMPA']
    #source_e_noise1.g_dyn = lambda: params_standard['E']['nu_ext'] * params_standard['E']['Cext'] * params_standard['E']['tau_AMPA']
    #source_e_noise1.noise_tau = params_standard['E']['tau_AMPA']
    pop_e1.add_noise(source_e_noise1)

    source_e_noise2 = MFStaticInput('E_noise2', pop_e2, params_standard['E']['Cext'], params_standard['E']['nu_ext'], {
        IP.GM: params_standard['E']['gAMPA'],
        IP.VREV: 0 * mV,
        IP.TAU: params_standard['E']['tau_AMPA'],
    })
    #source_e_noise2.g_base = params_standard['E']['gAMPA']
    #source_e_noise2.g_dyn = lambda: params_standard['E']['nu_ext'] * params_standard['E']['Cext'] * params_standard['E']['tau_AMPA']
    #source_e_noise2.noise_tau = params_standard['E']['tau_AMPA']
    pop_e2.add_noise(source_e_noise2)

    source_i_noise = MFStaticInput('I_noise', pop_i, params_standard['I']['Cext'], params_standard['I']['nu_ext'], {
        IP.GM: params_standard['I']['gAMPA'],
        IP.VREV: 0 * mV,
        IP.TAU: params_standard['I']['tau_AMPA'],
    })
    #source_i_noise.g_base = params_standard['I']['gAMPA']
    #source_i_noise.g_dyn = lambda: params_standard['I']['nu_ext'] * params_standard['I']['Cext'] * params_standard['I']['tau_AMPA']
    #source_i_noise.noise_tau = params_standard['I']['tau_AMPA']
    pop_i.add_noise(source_i_noise)

    # E->E NMDA
    syn_spec_nmda = {
        'tau_syn_rise': 1. * ms,
        'tau_syn_d1': 100. * ms,
        'tau_syn_d2': 100. * ms,
        'balance': .5,
        'tau_x': 150. * ms,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda1 = MFLinearInput('EE Nmda1', pop_e1, {
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: w_plus
    }, origin=pop_e1)#, synapse=syn_ee_nmda)
    source_ee_nmda12 = MFLinearInput('EE Nmda12', pop_e1, {
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: w_min
    }, origin=pop_e2)#, synapse=syn_ee_nmda)
    #source_ee_nmda1.g_base = params_standard['E']['gNMDA']
    #source_ee_nmda1.g_dyn = lambda: (
    #        ff * w_plus * syn_ee_nmda(pop_e1.rate_ms) + (1. - ff) * w_min * syn_ee_nmda(pop_e2.rate_ms)
    #    ) * pop_e1.n
    #source_ee_nmda1.is_nmda = True

    source_ee_nmda2 = MFLinearInput('EE Nmda2', pop_e2, {
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: w_min
    }, origin=pop_e1)#, synapse=syn_ee_nmda)
    source_ee_nmda22 = MFLinearInput('EE Nmda22', pop_e2, {
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: (w_plus * ff + (1. - 2. * ff) * w_min) / (1 - ff)
    }, origin=pop_e2)#, synapse=syn_ee_nmda)
    #source_ee_nmda2.g_base = params_standard['E']['gNMDA']
    #source_ee_nmda2.g_dyn = lambda: (
    #        ff * w_min * syn_ee_nmda(pop_e1.rate_ms) + (ff * w_plus + (1. - 2.*ff) * w_min) * syn_ee_nmda(pop_e2.rate_ms)
    #    ) * pop_e2.n
    #source_ee_nmda2.is_nmda = True

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFLinearInput('IE Nmda', pop_i, {
        IP.GM: params_standard['I']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: 1
    }, origin=pop_e1)#, synapse=syn_ie_nmda)
    source_ie_nmda2 = MFLinearInput('IE Nmda2', pop_i, {
        IP.GM: params_standard['I']['gNMDA'],
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],
        IP.W: 1
    }, origin=pop_e2)#, synapse=syn_ie_nmda)
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

    source_ii_gaba = MFLinearInput('II Gaba', pop_i, {
        IP.GM: params_standard['I']['gGABA'],
        IP.VREV: -70 * mV,
        IP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, origin=pop_i)#, synapse=syn_ii_gaba)
    #source_ii_gaba.g_base = params_standard['I']['gGABA']
    #source_ii_gaba.g_dyn = lambda: syn_ii_gaba(pop_i.rate_ms) * pop_i.n
    #source_ii_gaba.E_rev = -70.

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba1 = MFLinearInput('EI Gaba', pop_e1, {
        IP.GM: params_standard['E']['gGABA'],
        IP.VREV: -70 * mV,
        IP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, origin=pop_i)#, synapse=syn_ei_gaba)
    #source_ei_gaba1.g_base = params_standard['E']['gGABA']
    #source_ei_gaba1.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    #source_ei_gaba1.E_rev = -70.

    source_ei_gaba2 = MFLinearInput('EI Gaba', pop_e2, {
        IP.GM: params_standard['E']['gGABA'],
        IP.VREV: -70 * mV,
        IP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, origin=pop_i)#, synapse=syn_ei_gaba)
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
    pop_e = MFLinearPopulation('E', 800, {
        PP.GM: params_standard['E']['g_L'],
        PP.VL: params_standard['E']['V_L'],
        PP.CM: params_standard['E']['C_m'],
        PP.VTHR: params_standard['E']['V_th'],
        PP.VRES: params_standard['E']['V_reset'],
        PP.TAU_RP: params_standard['E']['tau_AMPA']
    })
    #pop_e.n = 800
    pop_i = MFLinearPopulation('I', 200, {
        PP.GM: params_standard['I']['g_L'],
        PP.VL: params_standard['I']['V_L'],
        PP.CM: params_standard['I']['C_m'],
        PP.VTHR: params_standard['I']['V_th'],
        PP.VRES: params_standard['I']['V_reset'],
        PP.TAU_RP: params_standard['I']['tau_AMPA']
    })
    #pop_i.n = 200
    pop_e.rate = initials['nu_e']
    pop_e.v_mean = -55. * mV
    pop_i.rate = initials['nu_i']
    pop_i.v_mean = -55. * mV

    system.populations += [pop_e, pop_i]

    # noise pops
    source_e_noise = MFStaticInput('E_noise', pop_e, params_standard['E']['Cext'], params_standard['E']['nu_ext'], {
        IP.GM: params_standard['E']['gAMPA'],
        IP.VREV: 0 * mV,
        IP.TAU: params_standard['E']['tau_AMPA'],
    })
    #source_e_noise.g_base = params_standard['E']['gAMPA']
    #source_e_noise.g_dyn = lambda: params_standard['E']['nu_ext'] * params_standard['E']['Cext'] * params_standard['E']['tau_AMPA']
    #source_e_noise.noise_tau = params_standard['E']['tau_AMPA']
    pop_e.add_noise(source_e_noise)

    source_i_noise = MFStaticInput('I_noise', pop_i, params_standard['I']['Cext'], params_standard['I']['nu_ext'], {
        IP.GM: params_standard['I']['gAMPA'],
        IP.VREV: 0 * mV,
        IP.TAU: params_standard['I']['tau_AMPA'],
    })
    #source_i_noise.g_base = params_standard['I']['gAMPA']
    #source_i_noise.g_dyn = lambda: params_standard['I']['nu_ext'] * params_standard['I']['Cext'] * params_standard['I']['tau_AMPA']
    #source_i_noise.noise_tau = params_standard['I']['tau_AMPA']
    pop_i.add_noise(source_i_noise)

    # E->E NMDA
    syn_spec_nmda = {
        'tau_syn_rise': 1. * ms,
        'tau_syn_d1': 100. * ms,
        'tau_syn_d2': 100. * ms,
        'balance': .5,
        'tau_x': 150. * ms,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda = MFLinearInput('EE Nmda', pop_e, {
        IP.GM: params_standard['E']['gNMDA'] * mult,
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],  # not sure
    }, origin=pop_e)#, synapse=syn_ee_nmda)
    #source_ee_nmda.g_base = params_standard['E']['gNMDA'] * mult
    #source_ee_nmda.g_dyn = lambda: syn_ee_nmda(pop_e.rate_ms) * pop_e.n
    #source_ee_nmda.is_nmda = has_nmda

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFLinearInput('IE Nmda', pop_i, {
        IP.GM: params_standard['I']['gNMDA'] * mult,
        IP.VREV: 0 * mV,
        IP.TAU: syn_spec_nmda['tau_syn_d1'],  # not sure
    }, origin=pop_e)#, synapse=syn_ie_nmda)
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
    source_ii_gaba = MFLinearInput('II Gaba', pop_i, {
        IP.GM: params_standard['I']['gGABA'],
        IP.VREV: -70 * mV,
        IP.TAU: syn_spec_gaba['tau_syn_d1'],  # not sure
    }, origin=pop_i)#, synapse=syn_ii_gaba)
    #source_ii_gaba.g_base = params_standard['I']['gGABA']
    #source_ii_gaba.g_dyn = lambda: syn_ii_gaba(pop_i.rate_ms) * pop_i.n
    #source_ii_gaba.E_rev = -70.

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba = MFLinearInput('EI Gaba', pop_e, {
        IP.GM: params_standard['E']['gGABA'],
        IP.VREV: -70 * mV,
        IP.TAU: syn_spec_gaba['tau_syn_d1'],  # not sure
    }, origin=pop_i)#, synapse=syn_ei_gaba)
    #source_ei_gaba.g_base = params_standard['E']['gGABA']
    #source_ei_gaba.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    #source_ei_gaba.E_rev = -70.

    return system
