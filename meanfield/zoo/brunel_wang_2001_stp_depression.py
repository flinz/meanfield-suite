from brian2.units import *

from meanfield.inputs.MFLinearNMDAInput import MFLinearNMDAInput
from meanfield.MFSystem import MFSystem
from meanfield.parameters import PP, IP
from meanfield.populations.MFLinearPopulation import MFLinearPopulation
from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.inputs.MFStaticInput import MFStaticInput
from meanfield.inputs.Synapses import Synapse

# BRUNEL, Nicolas et WANG, Xiao-Jing. Effects of neuromodulation in a cortical network model of object working memory dominated by recurrent inhibition. Journal of computational neuroscience, 2001, vol. 11, no 1, p. 63-85.
# TSODYKS, Misha, PAWELZIK, Klaus, et MARKRAM, Henry. Neural networks with dynamic synapses. Neural computation, 1998, vol. 10, no 4, p. 821-835.

params_standard = {
    'NMDA': {
        'gamma': 0.280112,
        'beta': 0.062,
    },
    'E': {
        'gamma': 0.280112,
        'beta': 0.062,
        'VE': 0. * mV,
        'VI': -70. * mV,
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
        'VE': 0. * mV,
        'VI': -70. * mV,
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

def no_subpopulation(has_nmda=True):

    # brunel 1999 system, one up state pop
    initials = {'nu_e': 3 * Hz, 'nu_i': 9 * Hz}

    mult = 0.25
    if has_nmda:
        mult = 1.

    ctor = MFLinearNMDAInput if has_nmda else MFLinearInput

    pop_e = MFLinearPopulation(800, {
        PP.GM: params_standard['E']['g_L'],
        PP.VL: params_standard['E']['V_L'],
        PP.CM: params_standard['E']['C_m'],
        PP.VTHR: params_standard['E']['V_th'],
        PP.VRES: params_standard['E']['V_reset'],
        PP.TAU_RP: params_standard['E']['t_ref']
    }, name='E')

    pop_i = MFLinearPopulation(200, {
        PP.GM: params_standard['I']['g_L'],
        PP.VL: params_standard['I']['V_L'],
        PP.CM: params_standard['I']['C_m'],
        PP.VTHR: params_standard['I']['V_th'],
        PP.VRES: params_standard['I']['V_reset'],
        PP.TAU_RP: params_standard['I']['t_ref']
    }, name='I')

    pop_e.rate = initials['nu_e']
    pop_e.v_mean = -51. * mV
    pop_i.rate = initials['nu_i']
    pop_i.v_mean = -55. * mV

    system = MFSystem(pop_e, pop_i, name='EI')

    # noise pops
    source_e_noise = MFStaticInput(params_standard['E']['Cext'], params_standard['E']['nu_ext'], pop_e, {
        IP.GM: params_standard['E']['gAMPA'],
        IP.VREV: params_standard['E']['VE'],
        IP.TAU: params_standard['E']['tau_AMPA'],
    }, name='E_noise')

    source_i_noise = MFStaticInput(params_standard['I']['Cext'], params_standard['I']['nu_ext'], pop_i, {
        IP.GM: params_standard['I']['gAMPA'],
        IP.VREV: params_standard['I']['VE'],
        IP.TAU: params_standard['I']['tau_AMPA'],
    }, name='I_noise')

    # E->E NMDA
    syn_spec_nmda = {
        'tau_syn_rise': 1. * ms,
        'tau_syn_d1': 100. * ms,
        'tau_syn_d2': 100. * ms,
        'balance': .5,
        'tau_x': 150. * ms,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda = ctor(pop_e, pop_e, {
        IP.BETA: params_standard['E']['beta'],
        IP.GAMMA: params_standard['E']['gamma'],
        IP.GM: params_standard['E']['gNMDA'] * mult,
        IP.VREV: params_standard['E']['VE'],
        IP.TAU: syn_spec_nmda['tau_syn_d1'],  # not sure
        IP.W: 1
    }, name='EE Nmda', synapse=syn_ee_nmda)

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = ctor(pop_e, pop_i, {
        IP.BETA: params_standard['I']['beta'],
        IP.GAMMA: params_standard['I']['gamma'],
        IP.GM: params_standard['I']['gNMDA'] * mult,
        IP.VREV: params_standard['I']['VE'],
        IP.TAU: syn_spec_nmda['tau_syn_d1'],  # not sure
        IP.W: 1
    }, name='IE Nmda', synapse=syn_ie_nmda)

    # I->I GABA
    syn_spec_gaba = {
        'tau_syn_rise': 0. * ms,
        'tau_syn_d1': 10. * ms,
        'tau_syn_d2': 10. * ms,
        'balance': .5,
    }
    syn_ii_gaba = Synapse(**syn_spec_gaba)
    source_ii_gaba = MFLinearInput(pop_i, pop_i, {
        IP.GM: params_standard['I']['gGABA'],
        IP.VREV:params_standard['I']['VI'],
        IP.TAU: syn_spec_gaba['tau_syn_d1'],  # not sure
    }, name='II Gaba', synapse=syn_ii_gaba)

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba = MFLinearInput(pop_i, pop_e, {
        IP.GM: params_standard['E']['gGABA'],
        IP.VREV: params_standard['E']['VI'],
        IP.TAU: syn_spec_gaba['tau_syn_d1'],  # not sure
    }, name='EI Gaba', synapse=syn_ei_gaba)

    return system


def one_subpopulation(w_plus_val=2.5):

    # brunel 1999 system, one up state pop
    ff = .1
    w_plus = w_plus_val
    w_min = 1. - ff * (w_plus - 1.) / (1. - ff)
    initials = {'nu_plus': 25 * Hz, 'nu_min': 1.5 * Hz, 'nu_i': 9 * Hz}

    pop_e1 = MFLinearPopulation(800 * ff, {
        PP.GM: params_standard['E']['g_L'],
        PP.VL: params_standard['E']['V_L'],
        PP.CM: params_standard['E']['C_m'],
        PP.VTHR: params_standard['E']['V_th'],
        PP.VRES: params_standard['E']['V_reset'],
        PP.TAU_RP: params_standard['E']['t_ref']
    }, name='Eup')

    pop_e2 = MFLinearPopulation(800 * (1. - ff), {
        PP.GM: params_standard['E']['g_L'],
        PP.VL: params_standard['E']['V_L'],
        PP.CM: params_standard['E']['C_m'],
        PP.VTHR: params_standard['E']['V_th'],
        PP.VRES: params_standard['E']['V_reset'],
        PP.TAU_RP: params_standard['E']['t_ref']
    }, name='Edown')

    pop_i = MFLinearPopulation(200, {
        PP.GM: params_standard['I']['g_L'],
        PP.VL: params_standard['I']['V_L'],
        PP.CM: params_standard['I']['C_m'],
        PP.VTHR: params_standard['I']['V_th'],
        PP.VRES: params_standard['I']['V_reset'],
        PP.TAU_RP: params_standard['I']['t_ref']
    }, name='I')

    pop_e1.rate = initials['nu_plus']
    pop_e1.v_mean = -51. * mV
    pop_e2.rate = initials['nu_min']
    pop_e2.v_mean = -55. * mV
    pop_i.rate = initials['nu_i']
    pop_i.v_mean = -55. * mV

    system = MFSystem(pop_e1, pop_e2, pop_i, name='Brunel')

    # noise pops

    source_e_noise1 = MFStaticInput(params_standard['E']['Cext'], params_standard['E']['nu_ext'], pop_e1, {
        IP.GM: params_standard['E']['gAMPA'],
        IP.VREV: params_standard['E']['VE'],
        IP.TAU: params_standard['E']['tau_AMPA'],
    }, name='E_noise1')

    source_e_noise2 = MFStaticInput(params_standard['E']['Cext'], params_standard['E']['nu_ext'], pop_e2, {
        IP.GM: params_standard['E']['gAMPA'],
        IP.VREV: params_standard['E']['VE'],
        IP.TAU: params_standard['E']['tau_AMPA'],
    }, name='E_noise2')

    source_i_noise = MFStaticInput(params_standard['I']['Cext'], params_standard['I']['nu_ext'], pop_i, {
        IP.GM: params_standard['I']['gAMPA'],
        IP.VREV: params_standard['I']['VE'],
        IP.TAU: params_standard['I']['tau_AMPA'],
    }, name='I_noise')

    # E->E NMDA
    syn_spec_nmda = {
        'tau_syn_rise': 1. * ms,
        'tau_syn_d1': 100. * ms,
        'tau_syn_d2': 100. * ms,
        'balance': .5,
        'tau_x': 150. * ms,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda1 = MFLinearNMDAInput(pop_e1, pop_e1, {
        IP.BETA: params_standard['E']['beta'],
        IP.GAMMA: params_standard['E']['gamma'],
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: params_standard['E']['VE'],
        IP.TAU: params_standard['E']['tau_AMPA'],
        IP.W: w_plus
    }, name='EE Nmda1', synapse=syn_ee_nmda)

    source_ee_nmda12 = MFLinearNMDAInput(pop_e2, pop_e1, {
        IP.BETA: params_standard['E']['beta'],
        IP.GAMMA: params_standard['E']['gamma'],
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: params_standard['E']['VE'],
        IP.TAU: params_standard['E']['tau_AMPA'],
        IP.W: w_min
    }, name='EE Nmda12', synapse=syn_ee_nmda)

    source_ee_nmda2 = MFLinearNMDAInput(pop_e1, pop_e2, {
        IP.BETA: params_standard['E']['beta'],
        IP.GAMMA: params_standard['E']['gamma'],
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: params_standard['E']['VE'],
        IP.TAU: params_standard['E']['tau_AMPA'],
        IP.W: w_min
    }, name='EE Nmda2', synapse=syn_ee_nmda)

    source_ee_nmda22 = MFLinearNMDAInput(pop_e2, pop_e2, {
        IP.BETA: params_standard['E']['beta'],
        IP.GAMMA: params_standard['E']['gamma'],
        IP.GM: params_standard['E']['gNMDA'],
        IP.VREV: params_standard['E']['VE'],
        IP.TAU: params_standard['E']['tau_AMPA'],
        IP.W: (w_plus * ff + (1. - 2. * ff) * w_min) / (1 - ff)
    }, name='EE Nmda22', synapse=syn_ee_nmda)

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFLinearNMDAInput(pop_e1, pop_i, {
        IP.BETA: params_standard['I']['beta'],
        IP.GAMMA: params_standard['I']['gamma'],
        IP.GM: params_standard['I']['gNMDA'],
        IP.VREV: params_standard['I']['VE'],
        IP.TAU: params_standard['I']['tau_AMPA'],
        IP.W: 1
    }, name='IE Nmda', synapse=syn_ie_nmda)

    source_ie_nmda2 = MFLinearNMDAInput(pop_e2, pop_i, {
        IP.BETA: params_standard['I']['beta'],
        IP.GAMMA: params_standard['I']['gamma'],
        IP.GM: params_standard['I']['gNMDA'],
        IP.VREV: params_standard['I']['VE'],
        IP.TAU: params_standard['I']['tau_AMPA'],
        IP.W: 1
    }, name='IE Nmda2', synapse=syn_ie_nmda)

    # I->I GABA
    syn_spec_gaba = {
        'tau_syn_rise': 0. * ms,
        'tau_syn_d1': 10. * ms,
        'tau_syn_d2': 10. * ms,
        'balance': .5,
    }
    syn_ii_gaba = Synapse(**syn_spec_gaba)

    source_ii_gaba = MFLinearInput(pop_i, pop_i, {
        IP.GM: params_standard['I']['gGABA'],
        IP.VREV: params_standard['I']['VI'],
        IP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, name='II Gaba', synapse=syn_ii_gaba)

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba1 = MFLinearInput(pop_i, pop_e1, {
        IP.GM: params_standard['E']['gGABA'],
        IP.VREV: params_standard['E']['VI'],
        IP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, name='EI Gaba', synapse=syn_ei_gaba)

    source_ei_gaba2 = MFLinearInput(pop_i, pop_e2, {
        IP.GM: params_standard['E']['gGABA'],
        IP.VREV: params_standard['E']['VI'],
        IP.TAU: syn_spec_gaba['tau_syn_d1'],
    }, name='EI Gaba', synapse=syn_ei_gaba)

    return system
