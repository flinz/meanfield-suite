from brian2 import Equations, Synapses, check_units

from MFSource import MFSource
from Utils import lazy
from params import SP


class MFNMDASource(MFSource):

    def __init__(self, name, pop, params, from_pop, synapse=None):
        super().__init__(name, pop, params, from_pop, synapse)
        defaults = {

        }
        expectations = {
            SP.TAU_NMDA: 1.,  # unitless
            SP.TAU_NMDA_RISE: 1.
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    @check_units(result=1)
    def g_dyn(self):
        activation = self.synapse(self.from_pop.rate) if self.synapse else self.from_pop.rate * self.params[SP.TAU]
        return self.from_pop.n * activation * self.params[SP.W]

    def brian2_model(self):
        return Equations(
            '''
            I = g * (v - ve) / (1 + gamma * exp(- beta * v) ) * s_post : amp
            s_post: 1
            ''',
            s_post=self.post_variable_name + '_post',
            I=self.current_name,
            g=self.params[SP.GM],
            ve=self.params[SP.VE]
        )
    # TODO parameterize, / mV
    # Beta 1 / mV

    @property
    def post_nonlinear_name(self):
        return 'x_' + self.ref

    @lazy
    def brian2(self, mode='i != j'):
        # weight
        model = Equations(
            '''
            w : 1
            s_post = w * s : 1 (summed)
            ds / dt = - s / tau_decay + alpha * x * (1 - s) : 1 (clock-driven)
            dx / dt = - x / tau_rise : 1 (clock-driven)
            ''',
            s_post=self.post_variable_name + '_post',
            s=self.post_variable_name,
            x=self.post_nonlinear_name,
            tau_decay=self.params[SP.TAU_NMDA],
            tau_rise=self.params[SP.TAU_NMDA_RISE],
            alpha=self.params[SP.ALPHA], # TODO ALPHA ? 1 / ms
        )
        eqs_pre = '''
        {} += 1
        '''.format(self.post_nonlinear_name)
        C = Synapses(self.from_pop.brian2, self.pop.brian2, method='euler', model=model, on_pre=eqs_pre)
        C.connect(mode)
        C.w[:] = 1
        return C

    def __repr__(self):
        return "MFSource [{}] <{}, nmda: {}, E_rev: {}>".format(id(self), self.name, self.is_nmda, self.params[SP.VE])
