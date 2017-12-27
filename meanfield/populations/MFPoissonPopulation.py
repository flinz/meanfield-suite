from types import MappingProxyType

from brian2 import check_units, units, PoissonGroup

from meanfield.parameters import PP
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import lazyproperty


class MFPoissonPopulation(MFPopulation):

    arguments = MappingProxyType({
        PP.GM: units.siemens,
        PP.VRES: units.volt,
        PP.TAU_RP: units.second
    })

    defaults = MappingProxyType({})

    @check_units(rate=units.Hz)
    def __init__(self, name, n, rate, parameters):
        super().__init__(name, n, parameters)
        self.rate = rate

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    @property
    def rate_prediction(self):
        return self.rate

    @property
    def v_mean_prediction(self):
        return None

    # Population

    @lazyproperty
    def brian2(self):
        return PoissonGroup(self.n, self.rate, name=self.name)

    def brian2_model(self):
        return None
