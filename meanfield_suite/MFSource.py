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
        self.g_base = 0. * units.siemens  # [nS]

        defaults = {}
        expectations = {
            SP.GM: units.siemens,
            SP.VE: units.volt,
            SP.TAU: units.second,
        }

        self.params = MFParams(params)
        self.params.fill(defaults)
        self.params.verify(expectations)

        self.g_base = params[SP.GM]  # TODO : parameterize, solve for ?

        # link to pop
        self.pop = pop
        self.pop.sources.append(self)  # TODO consistent

    @property
    @check_units(result=units.siemens)
    def conductance(self):
        ret = self.g_dyn() * self.g_base
        if self.is_nmda:
            return ret / self.pop.J * (
                1. + (1. - 1 / self.pop.J) * self.pop.params[NP.BETA] * (self.pop.v_mean - self.params[SP.VE])
            )
        return ret

    @property
    def voltage_conductance(self):
        cond = self.g_dyn() * self.g_base
        ret = cond * (self.params[SP.VE] - self.pop.params[NP.VL])
        if self.is_nmda:
            return 1. / self.pop.J * (
                ret + (1. - 1. / self.pop.J) * cond * self.pop.params[NP.BETA] * (self.pop.v_mean - self.params[SP.VE]) * (self.pop.v_mean - self.pop.params[NP.VL])
            )
        return ret

    def print_sys(self):
        print("{}: {}".format(self, self.conductance))

    def __repr__(self):
        return "MFSource [{}] <{}, nmda: {}, E_rev: {}>".format(id(self),self.name, self.is_nmda, self.params[SP.VE])

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
            tau=self.params[SP.TAU]
        )

class MFStaticSource(MFSource):

    @check_units(rate=units.hertz, n=1)
    def __init__(self, name, pop, params, rate, n):
        super().__init__(name, pop, params)
        self.rate = rate
        self.n = n

    @check_units(result=1)
    def g_dyn(self):
        return self.rate * self.n * self.params[SP.TAU]

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
        activation = self.synapse(self.from_pop.rate) if self.synapse else self.from_pop.rate * self.params[SP.TAU]
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


class MFNMDASource(MFDynamicSource):

    def __init__(self, name, pop, params, from_pop, synapse=None):
        super().__init__(name, pop, params, from_pop, synapse)
        defaults = {

        }
        expectations = {
            SP.TAU_NMDA: 1.,  # unitless
            SP.TAU_NMDA_RISE: 1.
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    @check_units(result=1)
    def g_dyn(self):
        activation = self.synapse(self.from_pop.rate) if self.synapse else self.from_pop.rate * self.params[SP.TAU]
        return self.from_pop.n * activation * self.params[SP.W]
        # TODO : SAME ?

    def brian2_model(self):
        return Equations(
            '''
            I = g * (v - ve) / (1 + exp(-0.062 * v) / 3.57) * s_post : amp
            s_post: 1
            ''',
            s_post=self.post_variable_name + '_post',
            I=self.current_name,
            g=self.params[SP.GM],
            ve=self.params[SP.VE]
        )
    # TODO parameterize, / mV

    @property
    def post_nonlinear_name(self):
        return 'x_' + self.ref
    # TODO name

    @lazy
    def brian2(self, mode='i != j'):
        # weight
        model = Equations(
            '''
            w : 1
            s_post = w * s : 1 (summed)
            ds / dt = - s / tau_decay + alpha * x * (1 - s) : 1 (clock-driven)
            dx / dt = - x / tau_rise : 1 (clock-driven)
            ''',
            s_post=self.post_variable_name + '_post',
            s=self.post_variable_name,
            x=self.post_nonlinear_name,
            tau_decay=self.params[SP.TAU_NMDA],
            tau_rise=self.params[SP.TAU_NMDA_RISE],
            alpha=self.params[SP.ALPHA], # TODO ALPHA ?
        )
        eqs_pre = '''
        {} += 1
        '''.format(self.post_nonlinear_name)
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

