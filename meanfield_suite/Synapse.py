import numpy as np
import pylab as pl

def stp_dgl_u(U, tau_f, tau_x, rate_ms):
    """ Differential equation equilibrium solution of short-term plastic synaptic input: facilitation variable u."""
    return float(U) * (1. + rate_ms * tau_f) * (1. + float(U) * rate_ms * tau_f)**(-1.)

def stp_dgl_x(U, tau_f, tau_x, rate_ms):
    """ Differential equation equilibrium solution of short-term plastic synaptic input: depression variable x."""
    return (1. + float(U) * rate_ms * tau_f) * (1. + float(U) * rate_ms * (tau_f + tau_x + rate_ms * tau_f * tau_x))**(-1.)

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
        pl.plot(x*1e3, [self(xv) for xv in x], **param)
        pl.xlim((0, 180))
        pl.legend(loc="best")

        pl.xlabel("Input rate [Hz]")
        pl.ylabel("Channel activation [1]")

    def make_fun(self):
        return lambda x: np.sum(self.taus) * x * self.stp_ur(x)

