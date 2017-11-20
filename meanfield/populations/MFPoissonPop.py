from brian2 import check_units, units, PoissonGroup

from meanfield.populations.MFPop import MFPop
from meanfield.utils import lazyproperty
from meanfield.parameters import NP


class MFPoissonPop(MFPop):
    @check_units(rate=units.Hz)
    def __init__(self, name, n, rate, params):
        super().__init__(name, n, params)
        self.rate = rate

        defaults = {
        }
        expectations = {
            NP.GM: units.siemens,
            NP.VRES: units.volt,
            NP.TAU_RP: units.second
        }
        self.params.fill(defaults)
        self.params.verify(expectations)

    @lazyproperty
    def brian2(self):
        return PoissonGroup(self.n, self.rate)

    @property
    def rate_prediction(self):
        return self.rate

    @property
    def v_mean_prediction(self):
        return None
