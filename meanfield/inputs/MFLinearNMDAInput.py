from types import MappingProxyType
from typing import Union, Dict

from brian2 import units, Equations, check_units, np

from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.parameters import IP, PP
from meanfield.parameters.MFParameters import MFParameters
from meanfield.populations.MFPopulation import MFPopulation


class MFLinearNMDAInput(MFLinearInput):

    arguments = MappingProxyType({
        IP.GM: units.siemens,
        IP.VREV: units.volt,
        IP.TAU: units.second,
        IP.W: 1,
        IP.BETA: 1,
        IP.GAMMA: 1,
    })

    defaults = MappingProxyType({
        IP.W: 1,
    })

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParameters], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    @property
    def J(self):
        """Linearization factor for NMDA"""

        # careful on mvolt conversion
        return 1 + self[IP.GAMMA] * np.exp(-self[IP.BETA] * self.target.v_mean / units.mV)

    @property
    @check_units(result=1)
    def rho1(self):
        return 1 / self.J

    @property
    @check_units(result=1)
    def rho2(self):
        return (self.J - 1) / (self.J ** 2) * self[IP.BETA] * (self.target.v_mean - self[IP.VREV]) / units.mV

    @property
    @check_units(result=units.siemens)
    def conductance(self) -> units.siemens:
        return self.g_dyn() * self[IP.GM] * (self.rho1 + self.rho2)

    @property
    @check_units(result=units.amp)
    def voltage_conductance(self) -> units.amp:
        return self.g_dyn() * self[IP.GM] * (
                self.rho1 * (self[IP.VREV] - self.target[PP.VL]) +
                self.rho2 * (self.target.v_mean - self.target[PP.VL])
        )

    # Simulation

    @property
    def brian2_model(self) -> Equations:
        return Equations(
            '''
            I = g * (v - vrev ) / (1 + gamma * exp(- beta * v)) * s : amp
            ds / dt = - s / tau_decay : 1
            ''',
            I=self.current_name,
            g=self[IP.GM],
            s=self.post_variable_name,
            s_post=self.post_variable_name + '_post',
            vrev=self[IP.VREV],
            tau_decay=self[IP.TAU],
            gamma=self[IP.GAMMA],
            beta=self[IP.BETA] / units.mV,
        )

