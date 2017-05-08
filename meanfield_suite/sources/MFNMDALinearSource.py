from brian2 import units, Equations, Synapses, check_units

from MFLinearSource import MFLinearSource
from MFSource import MFSource
from Utils import lazy
from params import SP, NP

# FIXME triple time scale source : MFTLinearTTSSource


class MFNMDALinearSource(MFLinearSource):

    def __init__(self, name, pop, params, from_pop):
        super().__init__(name, pop, params, from_pop)
        defaults = {

        }
        # FIXME synapse
        #
        expectations = {
            SP.TAU_NMDA: 1.,  # unitless
            SP.TAU_NMDA_RISE: 1.,
            SP.ALPHA: 1.
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    @property
    def J(self):
        """Linearization factor for NMDA"""
        return 1 + self.params[NP.GAMMA] * np.exp(-self.params[NP.BETA] * self.v_mean)

    @property
    def rho1(self):
        return 1 / self.pop.J

    @property
    def rho2(self):
        return (self.J - 1) / self.J ** 2 * self.pop.params[NP.BETA] * (self.pop.v_mean - self.params[SP.VREV])

    @property
    @check_units(result=units.siemens)
    def conductance(self):
        return self.g_dyn() * self.g_base * (self.rho1 + self.rho2)

    @property
    def voltage_conductance(self):
        return self.g_dyn() * self.g_base * (
            self.rho1 * (self.params[SP.VREV ] - self.pop.params[NP.VL]) +
            self.rho2 * (self.pop.v_mean - self.pop.params[NP.VREV])
        )

    def brian2_model(self):
        return Equations(
            '''
            I = g * (v - VREV ) / (1 + gamma * exp(- beta * v)) * s : amp
            ds / dt = - s / tau_decay : 1
            ''',
            s_post=self.post_variable_name + '_post',
            I=self.current_name,
            g=self.params[SP.GM],
            VREV =self.params[SP.VREV ],
            tau_decay=self.params[SP.TAU_NMDA],
            gamma=self.params[0], # TODO
            beta=self.params[0] / units.mV
        )

    # simulate poisson spike train (100 neurons) connect all with same source, record s (should be linear function of input), 10 rates * 10 reps
    # s = tau * nu

    ## s ~= g_dyn (theory)
    # TTS different slope but still linear


    def __repr__(self):
        return "MFSource [{}] <{}, nmda: {}, E_rev: {}>".format(id(self), self.name, self.is_nmda, self.params[SP.VREV ])
