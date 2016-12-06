from functools import partial

from MFConstraint import MFConstraint
from MFPop import MFOldPop
from MFSolver import MFSolver
from MFSource import MFSource
from MFState import MFState
from MFSystem import MFSystem

params_standard = {
    "NMDA": {
        "gamma": 0.280112,
        "beta": 0.062,
    },
    "E": {
        "gamma": 0.280112,
        "beta": 0.062,
        "VE": 0.,
        "E_L": -70.,
        "VI": -70.,
        "V_th": -50.,
        "V_reset": -60.,
        "tau_AMPA": 2.,
        "t_ref": 2.,
        "C_m": 500.,
        "g_L": 25.,
        "Cext": 1000,
        "nu_ext": 0.0024,
        "gAMPA": 2.08,
        "gNMDA": 0.327,
        "gGABA": 1.25
    },
    "I": {
        "gamma": 0.280112,
        "beta": 0.062,
        "VE": 0.,
        "E_L": -70.,
        "VI": -70.,
        "V_th": -50.,
        "V_reset": -60.,
        "tau_AMPA": 2.,
        "t_ref": 1.,
        "C_m": 200.,
        "g_L": 20.,
        "Cext": 1000,
        "nu_ext": 0.0024,
        "gAMPA": 1.62,
        "gNMDA": 0.258,
        "gGABA": 0.973
    }
}


def setup_brunel99(w_plus_val=2.5):

    # brunel 1999 system, one up state pop
    ff = .1
    w_plus = w_plus_val
    w_min = 1. - ff*(w_plus - 1.)/(1. - ff)
    initials = {'nu_plus': .025, 'nu_min': .0015, 'nu_i': .009}

    system = MFSystem("Brunel")
    pop_e1 = MFOldPop("Eup", params_standard["E"])
    pop_e1.n = 800
    pop_e2 = MFOldPop("Edown", params_standard["E"])
    pop_e2.n = 800
    pop_i = MFOldPop("I", params_standard["I"])
    pop_i.n = 200
    pop_e1.rate_ms = initials["nu_plus"]
    pop_e1.v_mean = -51.
    pop_e2.rate_ms = initials["nu_min"]
    pop_e2.v_mean = -55.
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
        "tau_syn_rise": 1.,
        "tau_syn_d1": 100.,
        "tau_syn_d2": 100.,
        "balance": .5,
        "tau_x": 150.,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda1 = MFSource('EE Nmda1', pop_e1)
    source_ee_nmda1.g_base = params_standard["E"]["gNMDA"]
    source_ee_nmda1.g_dyn = lambda: (
            ff * w_plus * syn_ee_nmda(pop_e1.rate_ms) + (1. - ff) * w_min * syn_ee_nmda(pop_e2.rate_ms)
        ) * pop_e1.n
    source_ee_nmda1.is_nmda = True

    source_ee_nmda2 = MFSource('EE Nmda2', pop_e2)
    source_ee_nmda2.g_base = params_standard["E"]["gNMDA"]
    source_ee_nmda2.g_dyn = lambda: (
            ff * w_min * syn_ee_nmda(pop_e1.rate_ms) + (ff * w_plus + (1. - 2.*ff) * w_min) * syn_ee_nmda(pop_e2.rate_ms)
        ) * pop_e2.n
    source_ee_nmda2.is_nmda = True

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFSource('IE Nmda', pop_i)
    source_ie_nmda.g_base = params_standard["I"]["gNMDA"]
    source_ie_nmda.g_dyn = lambda: (
            ff * syn_ie_nmda(pop_e1.rate_ms) + (1. - ff) * syn_ie_nmda(pop_e2.rate_ms)
        ) * pop_e1.n
    source_ie_nmda.is_nmda = True

    # I->I GABA
    syn_spec_gaba = {
        "tau_syn_rise": 0.,
        "tau_syn_d1": 10.,
        "tau_syn_d2": 10.,
        "balance": .5,
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
    source_ei_gaba1.E_rev = -70.
    source_ei_gaba2 = MFSource('EI Gaba', pop_e2)
    source_ei_gaba2.g_base = params_standard["E"]["gGABA"]
    source_ei_gaba2.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    source_ei_gaba2.E_rev = -70.

    return system


def setup_EI(has_nmda=True):

    # brunel 1999 system, one up state pop
    initials = {'nu_e': .003, 'nu_i': .009}

    mult = 0.25
    if has_nmda:
        mult = 1.

    system = MFSystem("EI")
    pop_e = MFOldPop("E", params_standard["E"])
    pop_e.n = 800
    pop_i = MFOldPop("I", params_standard["I"])
    pop_i.n = 200
    pop_e.rate_ms = initials["nu_e"]
    pop_e.v_mean = -51.
    pop_i.rate_ms = initials["nu_i"]
    pop_i.v_mean = -55.

    system.pops += [pop_e, pop_i]

    # noise pops
    source_e_noise = MFSource("E_noise", pop_e)
    source_e_noise.g_base = params_standard["E"]["gAMPA"]
    source_e_noise.g_dyn = lambda: params_standard["E"]["nu_ext"] * params_standard["E"]["Cext"] * params_standard["E"]["tau_AMPA"]
    source_e_noise.noise_tau = params_standard["E"]["tau_AMPA"]
    pop_e.noise = source_e_noise

    source_i_noise = MFSource("I_noise", pop_i)
    source_i_noise.g_base = params_standard["I"]["gAMPA"]
    source_i_noise.g_dyn = lambda: params_standard["I"]["nu_ext"] * params_standard["I"]["Cext"] * params_standard["I"]["tau_AMPA"]
    source_i_noise.noise_tau = params_standard["I"]["tau_AMPA"]
    pop_i.noise = source_i_noise

    # E->E NMDA
    syn_spec_nmda = {
        "tau_syn_rise": 1.,
        "tau_syn_d1": 100.,
        "tau_syn_d2": 100.,
        "balance": .5,
        "tau_x": 150.,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda = MFSource('EE Nmda', pop_e)
    source_ee_nmda.g_base = params_standard["E"]["gNMDA"] * mult
    source_ee_nmda.g_dyn = lambda: syn_ee_nmda(pop_e.rate_ms) * pop_e.n
    source_ee_nmda.is_nmda = has_nmda

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFSource('IE Nmda', pop_i)
    source_ie_nmda.g_base = params_standard["I"]["gNMDA"] * mult
    source_ie_nmda.g_dyn = lambda: syn_ie_nmda(pop_e.rate_ms) * pop_e.n
    source_ie_nmda.is_nmda = has_nmda

    # I->I GABA
    syn_spec_gaba = {
        "tau_syn_rise": 0.,
        "tau_syn_d1": 10.,
        "tau_syn_d2": 10.,
        "balance": .5,
    }
    syn_ii_gaba = Synapse(**syn_spec_gaba)
    source_ii_gaba = MFSource('II Gaba', pop_i)
    source_ii_gaba.g_base = params_standard["I"]["gGABA"]
    source_ii_gaba.g_dyn = lambda: syn_ii_gaba(pop_i.rate_ms) * pop_i.n
    source_ii_gaba.E_rev = -70.

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba = MFSource('EI Gaba', pop_e)
    source_ei_gaba.g_base = params_standard["E"]["gGABA"]
    source_ei_gaba.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    source_ei_gaba.E_rev = -70.

    return system

import numpy as np
import matplotlib.pyplot as plt

def stp_dgl_u(U, tau_f, tau_x, rate_ms):
    """ Differential equation equilibrium solution of short-term pltastic synaptic input: facilitation variable u."""
    return float(U) * (1. + rate_ms * tau_f) * (1. + float(U) * rate_ms * tau_f) ** (-1.)

def stp_dgl_x(U, tau_f, tau_x, rate_ms):
    """ Differential equation equilibrium solution of short-term pltastic synaptic input: depression variable x."""
    return (1. + float(U) * rate_ms * tau_f) * (1. + float(U) * rate_ms * (tau_f + tau_x + rate_ms * tau_f * tau_x)) ** (-1.)


class Synapse(object):
    """
        Synapse with 3 time constants and possibly depression/facilitation.
    """
    fun = None

    def __init__(self, tau_syn_rise=1., tau_syn_d1=100., tau_syn_d2=100., balance=.5, U=1., tau_f=1000., tau_x=0.):

        # facilitation & depression parameters
        self.tau_f = tau_f
        self.tau_x = tau_x
        self.U = U

        # synaptic timeconstants parameters
        self.tau_syn_rise = tau_syn_rise
        self.tau_syn_d1 = tau_syn_d1
        self.tau_syn_d2 = tau_syn_d2
        self.balance = balance

    @property
    def taus(self):
        return np.array([
            self.balance * self.tau_syn_d1,
            (1.-self.balance) * self.tau_syn_d2,
            - self.balance * self.tau_syn_d1 * self.tau_syn_rise / (self.tau_syn_d1 + self.tau_syn_rise),
            - (1.-self.balance) * self.tau_syn_d2 * self.tau_syn_rise / (self.tau_syn_d2 + self.tau_syn_rise)
        ])

    def __call__(self, x):
        if self.fun is None:
            self.fun = self.make_fun()
        return self.fun(x)

    def stp_u(self, x):
        return stp_dgl_u(self.U, self.tau_f, self.tau_x, x)

    def stp_x(self, x):
        return stp_dgl_x(self.U, self.tau_f, self.tau_x, x)

    def stp_ur(self, x):
        return self.stp_x(x) * self.stp_u(x)

    def plot(self, **param):
        try:
            param['label']
        except KeyError:
            param['label'] = self.__repr__()

        x = np.arange(1., 181., 1.) * 1e-3
        plt.plot(x*1e3, [self(xv) for xv in x], **param)
        plt.xlim((0, 180))
        plt.legend(loc="best")

        plt.xlabel("Input rate [Hz]")
        plt.ylabel("Channel activation [1]")

    def make_fun(self):
        return lambda x: np.sum(self.taus) * x * self.stp_ur(x)
