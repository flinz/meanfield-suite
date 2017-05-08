from brian2 import units, check_units

from MFSource import MFSource
from Utils import lazy
from params import SP


class MFStaticSource(MFSource):

    @check_units(rate=units.hertz, n=1)
    def __init__(self, name, pop, params, rate, n):
        super().__init__(name, pop, params)

        self.rate = rate
        self.n = n

    @check_units(result=1)
    def g_dyn(self):
        return self.rate * self.n * self.params[SP.TAU]

    @lazy
    def b2_syn(self):
        # TODO add PoissonInput
        return None


