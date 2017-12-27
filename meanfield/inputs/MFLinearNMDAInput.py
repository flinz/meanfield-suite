from types import MappingProxyType
from typing import Union, Dict

from brian2 import units, Equations, check_units, np

from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.parameters import IP, PP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.parameters import Connection
from meanfield.parameters.Connection import ConnectionStrategy


class MFLinearNMDAInput(MFLinearInput):

    arguments = MappingProxyType({
        IP.TAU_NMDA: units.second,
        IP.ALPHA: 1,
        IP.BETA: 1,  # TODO beta 1/v ?
        IP.GAMMA: 1,
    })

    defaults = MappingProxyType({})

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParams], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    @property
    def J(self):
        """Linearization factor for NMDA"""
        return 1 + self.parameters[PP.GAMMA] * np.exp(-self.parameters[PP.BETA] * self.target.v_mean / units.volt)

    @property
    @check_units(result=1)
    def rho1(self):
        return 1 / self.J

    @property
    @check_units(result=1)
    def rho2(self):
        return (self.J - 1) / self.J ** 2 * self.parameters[PP.BETA] * (self.target.v_mean - self.parameters[IP.VREV]) / units.volt # TODO unitless?

    @property
    @check_units(result=units.siemens)
    def conductance(self):
        return self.g_dyn() * self.g_base * (self.rho1 + self.rho2)

    @property
    def voltage_conductance(self):
        return self.g_dyn() * self.g_base * (
                self.rho1 * (self.parameters[IP.VREV] - self.target.params[PP.VL]) +
                self.rho2 * (self.target.v_mean - self.target.params[PP.VL])
        )

    # Simulation

    def brian2_model(self):
        return Equations(
            '''
            I = g * (v - vrev ) / (1 + gamma * exp(- beta * v)) * s : amp
            ds / dt = - s / tau_decay : 1
            ''',
            I=self.current_name,
            g=self.parameters[IP.GM],
            s=self.post_variable_name,
            s_post=self.post_variable_name + '_post',
            vrev=self.parameters[IP.VREV],
            tau_decay=self.parameters[IP.TAU_NMDA],
            gamma=self.parameters[IP.GAMMA],
            beta=self.parameters[IP.BETA] / units.mV,
        )

