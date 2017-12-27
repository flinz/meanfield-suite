from typing import Union, Dict

from brian2 import units, Equations, check_units, np

from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.parameters import IP, PP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPop import MFPop
from meanfield.parameters import Connection
from meanfield.parameters.Connection import ConnectionStrategy


class MFLinearNMDAInput(MFLinearInput):

    def __init__(self, name: str, pop: MFPop, params: Union[Dict, MFParams], from_pop: MFPop, connection: ConnectionStrategy=Connection.all_to_all()):
        super().__init__(name, pop, params, from_pop)

        defaults = {}
        expectations = {
            IP.TAU_NMDA: units.second,
            IP.ALPHA: 1,
            IP.BETA: 1, # TODO beta 1/v ?
            IP.GAMMA: 1,
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    @property
    def J(self):
        """Linearization factor for NMDA"""
        return 1 + self.params[PP.GAMMA] * np.exp(-self.params[PP.BETA] * self.pop.v_mean / units.volt)

    @property
    @check_units(result=1)
    def rho1(self):
        return 1 / self.J

    @property
    @check_units(result=1)
    def rho2(self):
        return (self.J - 1) / self.J ** 2 * self.params[PP.BETA] * (self.pop.v_mean - self.params[IP.VREV]) / units.volt # TODO unitless?

    @property
    @check_units(result=units.siemens)
    def conductance(self):
        return self.g_dyn() * self.g_base * (self.rho1 + self.rho2)

    @property
    def voltage_conductance(self):
        return self.g_dyn() * self.g_base * (
                self.rho1 * (self.params[IP.VREV] - self.pop.params[PP.VL]) +
                self.rho2 * (self.pop.v_mean - self.pop.params[PP.VL])
        )

    def brian2_model(self):
        return Equations(
            '''
            I = g * (v - vrev ) / (1 + gamma * exp(- beta * v)) * s : amp
            ds / dt = - s / tau_decay : 1
            ''',
            I=self.current_name,
            g=self.params[IP.GM],
            s=self.post_variable_name,
            s_post=self.post_variable_name + '_post',
            vrev=self.params[IP.VREV],
            tau_decay=self.params[IP.TAU_NMDA],
            gamma=self.params[IP.GAMMA],
            beta=self.params[IP.BETA] / units.mV,
        )

