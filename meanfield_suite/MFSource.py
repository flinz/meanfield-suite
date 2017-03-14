from abc import abstractmethod

import numpy as np
import matplotlib.pyplot as plt
from brian2 import units, Equations, Synapses, check_units

from MFParams import MFParams
from Utils import lazy, name2identifier
from params import SP, NP

class MFSource(object):
    """Source: a synapse coupled to pops"""

    def __init__(self, name, pop, params):

        self.name = name
        self.ref = name2identifier(name)
        self.is_nmda = False
        self.g_base = 0.  # [nS]
        self.E_rev = 0.   # [mV] excitatory by default
        self.noise_tau = 0.   # [mV] excitatory by default

        defaults = {}
        expectations = {
            SP.GM: units.siemens,
            SP.VE: units.volt,
            SP.TAU_M: units.second,
        }
        # TODO : which is static/dynamic

        self.params = MFParams(params)
        self.params.fill(defaults)
        self.params.verify(expectations)

        self.g_base = params[SP.GM]
        self.E_rev = params[SP.VE]
        self.noise_tau = params[SP.TAU_M]

        # link to pop
        self.pop = pop
        self.pop.sources.append(self) # TODO consistent

    @property
    @check_units(result=units.siemens)
    def conductance(self):
        tmp = self.g_dyn() * self.g_base
        # print('__')
        #if self.is_nmda:
        #    J_ = 1./self.pop.J
        #    return tmp * J_ * (
        #        1. + (1.-J_) * self.pop.params[NP.BETA] * (self.pop.v_mean - self.E_rev)
        #    )
        return tmp

    @property
    def voltage_conductance(self):
        cond = self.g_dyn() * self.g_base
        tmp = cond * (self.E_rev - self.pop.params[NP.VL])
        # TODO should get order 1 nA
        #if self.is_nmda:
        #    J_ = 1. / self.pop.J
        #    return J_ * (
        #        tmp + (1.-J_) * cond * self.pop.params[NP.BETA] * (self.pop.v_mean - self.E_rev) * (self.pop.v_mean - self.pop.params["E_L"])
        #    )
        return tmp

    def print_sys(self):
        print("%s: %f" % (self, self.conductance))

    def __repr__(self):
        return "MFSource [%s] <%s, nmda: %s, E_rev: %.1f>" % (id(self),self.name, self.is_nmda, self.E_rev)

    @abstractmethod
    def g_dyn(self):
        pass

    @abstractmethod
    def brian2(self):
        pass

    @property
    def current_name(self):
        return 'I_' + self.ref

    @property
    def post_variable_name(self):
        return 's_' + self.ref

    def brian2_model(self):
        return Equations(
            '''
            I = g * (v - ve) * s : amp
            ds / dt = - s / tau : 1
            ''',
            s=self.post_variable_name,
            I=self.current_name,
            g=self.params[SP.GM],
            ve=self.params[SP.VE],
            tau=self.params[SP.TAU_M]
        )

        # TODO store post-corresponding var s_{}


class MFStaticSource(MFSource):

    @check_units(rate=units.hertz, n=1)
    def __init__(self, name, pop, params, rate, n):
        super().__init__(name, pop, params)
        self.rate = rate
        self.n = n
        # TODO real param

    @check_units(result=1)
    def g_dyn(self):
        return self.rate * self.n * self.params[SP.TAU_M]

    @lazy
    def brian2(self):
        return None


class MFDynamicSource(MFSource):

    def __init__(self, name, pop, params, from_pop, synapse=None):
        super().__init__(name, pop, params)
        self.from_pop = from_pop
        self.synapse = synapse
        defaults = {
            SP.W: 1.,
            SP.FRAC: 1.
        }
        expectations = {
            SP.W: 1., # unitless
            SP.FRAC: 1.
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    @check_units(result=1)
    def g_dyn(self):
        activation = self.synapse(self.from_pop.rate) if self.synapse else self.from_pop.rate * self.params[SP.TAU_M]
        return self.from_pop.n * activation * self.params[SP.W]

    @lazy
    def brian2(self, mode='i != j'):
        model = Equations('w : 1')
        eqs_pre = '''
        {} += w
        '''.format(self.post_variable_name)
        C = Synapses(self.from_pop.brian2, self.pop.brian2, method='euler', model=model, on_pre=eqs_pre)
        C.connect(mode)
        C.w[:] = 1
        return C







@check_units(result=1)
def stp_dgl_u(U, tau_f, tau_x, rate):
    """ Differential equation equilibrium solution of short-term pltastic synaptic input: facilitation variable u."""
    return float(U) * (1. + rate * tau_f) * (1. + float(U) * rate * tau_f) ** (-1.)

@check_units(result=1)
def stp_dgl_x(U, tau_f, tau_x, rate):
    """ Differential equation equilibrium solution of short-term pltastic synaptic input: depression variable x."""
    return (1. + float(U) * rate * tau_f) * (1. + float(U) * rate * (tau_f + tau_x + rate * tau_f * tau_x)) ** (-1.)


class Synapse(object):
    """
        Synapse with 3 time constants and possibly depression/facilitation.
    """
    fun = None

    @check_units(tau_syn_rise=units.second, tau_syn_d1=units.second, tau_syn_d2=units.second, tau_f=units.second, tau_x=units.second)
    def __init__(self, tau_syn_rise=1. * units.ms, tau_syn_d1=100. * units.ms, tau_syn_d2=100. * units.ms, balance=.5, U=1., tau_f=1000. * units.ms, tau_x=0. * units.ms):

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
    @check_units(result=units.second)
    def taus(self):
        return [
            self.balance * self.tau_syn_d1,
            (1.-self.balance) * self.tau_syn_d2,
            - self.balance * self.tau_syn_d1 * self.tau_syn_rise / (self.tau_syn_d1 + self.tau_syn_rise),
            - (1.-self.balance) * self.tau_syn_d2 * self.tau_syn_rise / (self.tau_syn_d2 + self.tau_syn_rise)
        ]

    def __call__(self, x):
        if self.fun is None:
            self.fun = self.make_fun()
        return self.fun(x)

    @check_units(result=1)
    def stp_u(self, x):
        return stp_dgl_u(self.U, self.tau_f, self.tau_x, x)

    @check_units(result=1)
    def stp_x(self, x):
        return stp_dgl_x(self.U, self.tau_f, self.tau_x, x)

    @check_units(result=1)
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
        plt.legend(loc='best')

        plt.xlabel('Input rate [Hz]')
        plt.ylabel('Channel activation [1]')

    def make_fun(self):
        return check_units(result=1)(lambda x: sum(self.taus) * x * self.stp_ur(x))

