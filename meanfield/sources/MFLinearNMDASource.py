from brian2 import units, Equations, check_units, np

from .MFLinearSource import MFLinearSource
from ..params import SP, NP


class MFLinearNMDASource(MFLinearSource):

    def __init__(self, name, pop, params, from_pop):
        super().__init__(name, pop, params, from_pop)

        defaults = {}
        expectations = {
            SP.TAU_NMDA: units.second,
            SP.ALPHA: 1,
            SP.BETA: 1, # TODO beta 1/v ?
            SP.GAMMA: 1,
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    @property
    def J(self):
        """Linearization factor for NMDA"""
        return 1 + self.params[NP.GAMMA] * np.exp(-self.params[NP.BETA] * self.pop.v_mean / units.volt)

    @property
    @check_units(result=1)
    def rho1(self):
        return 1 / self.J

    @property
    @check_units(result=1)
    def rho2(self):
        return (self.J - 1) / self.J ** 2 * self.params[NP.BETA] * (self.pop.v_mean - self.params[SP.VREV]) / units.volt # TODO unitless?

    @property
    @check_units(result=units.siemens)
    def conductance(self):
        return self.g_dyn() * self.g_base * (self.rho1 + self.rho2)

    @property
    def voltage_conductance(self):
        return self.g_dyn() * self.g_base * (
            self.rho1 * (self.params[SP.VREV] - self.pop.params[NP.VL]) +
            self.rho2 * (self.pop.v_mean - self.pop.params[NP.VL])
        )

    def b2_dyn(self):
        return Equations(
            '''
            I = g * (v - vrev ) / (1 + gamma * exp(- beta * v)) * s : amp
            ds / dt = - s / tau_decay : 1
            ''',
            I=self.current_name,
            g=self.params[SP.GM],
            s=self.post_variable_name,
            s_post=self.post_variable_name + '_post',
            vrev=self.params[SP.VREV],
            tau_decay=self.params[SP.TAU_NMDA],
            gamma=self.params[SP.GAMMA],
            beta=self.params[SP.BETA] / units.mV,
        )

