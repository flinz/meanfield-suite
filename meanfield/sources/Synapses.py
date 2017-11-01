from brian2 import check_units, units, plt, np


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
            (1. - self.balance) * self.tau_syn_d2,
            - self.balance * self.tau_syn_d1 * self.tau_syn_rise / (self.tau_syn_d1 + self.tau_syn_rise),
            - (1. - self.balance) * self.tau_syn_d2 * self.tau_syn_rise / (self.tau_syn_d2 + self.tau_syn_rise)
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

