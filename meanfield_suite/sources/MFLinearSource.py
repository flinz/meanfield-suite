from brian2 import check_units, Equations, Synapses

from MFSource import MFSource
from Utils import lazy
from params import SP


class MFLinearSource(MFSource):

    def __init__(self, name, pop, params, from_pop):
        super().__init__(name, pop, params)
        self.from_pop = from_pop
        defaults = {
            SP.W: 1.,
        }
        expectations = {
            SP.W: 1.,
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    @check_units(result=1)
    def g_dyn(self):
        return self.from_pop.n * self.from_pop.rate * self.params[SP.TAU] * self.params[SP.W]

    @lazy
    def brian2(self, mode='i != j'):
        model = Equations('w : 1')
        eqs_pre = '''
        {} += w
        '''.format(self.post_variable_name)
        C = Synapses(self.from_pop.brian2, self.pop.brian2, method='euler', model=model, on_pre=eqs_pre)
        C.connect(mode)
        C.w[:] = 1
        return C

    def brian2_model(self):
        return Equations(
            ## FIXME split current outside and s also
            '''
            I = g * (v - vrev) * s : amp
            ds / dt = - s / tau : 1
            ''',
            s=self.post_variable_name,
            I=self.current_name,
            g=self.params[SP.GM],
            vrev=self.params[SP.VREV],
            tau=self.params[SP.TAU]
        )


    def __repr__(self):
        return "MFSource [{}] <{}, nmda: {}, E_rev: {}>".format(id(self), self.name, self.is_nmda, self.params[SP.VE])