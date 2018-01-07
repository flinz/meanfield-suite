from types import MappingProxyType
from typing import Union, Dict

from brian2 import check_units, Equations, Synapses, BrianObject, units

from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.parameters import IP
from meanfield.parameters.MFParameters import MFParameters
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import lazyproperty


def stp_dgl_u(U, tau_f, tau_d, rate):
    """
    Differential equation equilibrium solution of short-term pltastic synaptic input: facilitation variable u.

    :param U:
    :param tau_f:
    :param tau_x:
    :param rate:
    :return:
    """
    return float(U) * (1. + rate * tau_f) * (1. + float(U) * rate * tau_f) ** (-1.)


def stp_dgl_x(U, tau_f, tau_d, rate):
    """
    Differential equation equilibrium solution of short-term pltastic synaptic input: depression variable x.

    :param U:
    :param tau_f:
    :param tau_d:
    :param rate:
    :return:
    """
    return (1. + float(U) * rate * tau_f) * (1. + float(U) * rate * (tau_f + tau_d + rate * tau_f * tau_d)) ** (-1.)


class MFLinearSTPInput(MFLinearInput):

    arguments = MappingProxyType({
        IP.GM: units.siemens,
        IP.VREV: units.volt,
        IP.TAU: units.second,
        IP.W: 1,
        IP.TAU_F: units.second,
        IP.TAU_D: units.second,
        IP.U: 1
    })

    defaults = MappingProxyType({
        IP.W: 1,
    })

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParameters], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    @check_units(result=1)
    def g_dyn(self):
        stp_x = stp_dgl_x(self[IP.U], self[IP.TAU_F], self[IP.TAU_D], self.origin.rate)
        stp_u = stp_dgl_u(self[IP.U], self[IP.TAU_F], self[IP.TAU_D], self.origin.rate)
        activation = self.origin.rate * stp_x * stp_u
        return self.connection.theory(self.origin.n) * self[IP.TAU] * activation * self[IP.W]

    # Simulation

    @lazyproperty
    def brian2(self) -> BrianObject:
        model = Equations('''
        x : 1
        u : 1
        w : 1
        ''')

        U = self[IP.U]
        tauf = self[IP.TAU_F]
        taud = self[IP.TAU_D]

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

