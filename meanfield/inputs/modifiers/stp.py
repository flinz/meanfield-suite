from brian2 import check_units, units


@check_units(result=1)
def stp_dgl_u(U, tau_f, tau_x, rate):
    """
    Differential equation equilibrium solution of short-term pltastic synaptic input: facilitation variable u.

    :param U:
    :param tau_f:
    :param tau_x:
    :param rate:
    :return:
    """
    return float(U) * (1. + rate * tau_f) * (1. + float(U) * rate * tau_f) ** (-1.)


@check_units(result=1)
def stp_dgl_x(U, tau_f, tau_x, rate):
    """
    Differential equation equilibrium solution of short-term pltastic synaptic input: depression variable x.

    :param U:
    :param tau_f:
    :param tau_x:
    :param rate:
    :return:
    """
    return (1. + float(U) * rate * tau_f) * (1. + float(U) * rate * (tau_f + tau_x + rate * tau_f * tau_x)) ** (-1.)


class STP(object):
    """
    Synapse with 3 time constants and possibly depression/facilitation.
    """

    params = {
        'tau_syn_rise':  1. * units.ms,
        'tau_syn_d1':  100. * units.ms,
        'tau_syn_d2':  100. * units.ms,
        'balance':  .5,
        'U':  1.,
        'tau_f':  1000. * units.ms,
        'tau_x': 0. * units.ms
    }

    @check_units(tau_syn_rise=units.second, tau_syn_d1=units.second, tau_syn_d2=units.second, tau_f=units.second, tau_x=units.second)
    def __init__(self, tau_syn_rise, tau_syn_d1, tau_syn_d2, balance, U, tau_f, tau_x):

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

    @check_units(result=1)
    def __call__(self, x):
        stp_x = stp_dgl_x(self.U, self.tau_f, self.tau_x, x)
        stp_u = stp_dgl_u(self.U, self.tau_f, self.tau_x, x)
        return sum(self.taus) * x * stp_x * stp_u

