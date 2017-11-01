from brian2 import check_units, units, PoissonGroup

from meanfield.populations.MFPop import MFPop
from meanfield.utils import lazyproperty


class MFPoissonSource(MFPop):
    @check_units(rate=units.Hz)
    def __init__(self, name, n, rate):
        super().__init__(name, n, {})
        self.rate = rate

    @lazyproperty
    def brian2(self):
        return PoissonGroup(self.n, self.rate)

    @property
    def rate_prediction(self):
        return self.rate

    @property
    def v_mean_prediction(self):
        return None
