from brian2 import units, check_units, PoissonInput

from meanfield.sources.MFSource import MFSource
from meanfield.utils import lazyproperty
from meanfield.parameters import SP


class MFStaticSource(MFSource):

    @check_units(rate=units.hertz, n=1)
    def __init__(self, name, pop, params, rate, n):
        super().__init__(name, pop, params)

        self.rate = rate
        self.n = n

    @check_units(result=1)
    def g_dyn(self):
        return self.rate * self.n * self.params[SP.TAU]

    @lazyproperty
    def b2_syn(self):
        return PoissonInput(self.pop, self.post_variable_name, self.n, self.rate, 1)
