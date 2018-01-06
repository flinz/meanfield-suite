from types import MappingProxyType
from typing import Union, Dict

from brian2 import check_units, Equations, Synapses, BrianObject, units

from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.parameters import IP
from meanfield.parameters.MFParameters import MFParameters
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import lazyproperty


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


class MFLinearSTPInput(MFLinearInput):

    arguments = MappingProxyType({

    })

    defaults = MappingProxyType({})

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParameters], synapse, **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

        self.synapse = STP(**synapse)

    # Theory

    @check_units(result=1)
    def g_dyn(self):
        return self.connection.theory(self.origin.n) * self.synapse(self.origin.rate) * self[IP.W]

    # Simulation

    @lazyproperty
    def brian2(self) -> BrianObject:
        model = Equations('''
        x : 1
        u : 1
        w : 1
        ''')

        U = self.synapse.U
        tauf = self.synapse.tau_f
        taud = self.synapse.tau_x

        on_pre = f'''
        u = {U} + (u - {U}) * exp(- (t - lastupdate) / ({tauf / units.ms} * ms))
        x = 1 + (x - 1) * exp(- (t - lastupdate) / ({taud / units.ms} * ms))

        {self.post_variable_name} += w * u * x

        x *= (1 - u)
        u += {U} * (1 - u)
        '''
        syn = Synapses(
            source=self.origin.brian2,
            target=self.target.brian2,
            method='euler',
            model=model,
            on_pre=on_pre,
            name=self.ref
        )
        self.connection.simulation(syn)
        syn.w[:] = self[IP.W]
        syn.x[:] = 1
        syn.u[:] = 1
        return syn

