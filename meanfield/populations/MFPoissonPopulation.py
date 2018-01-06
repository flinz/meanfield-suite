from types import MappingProxyType
from typing import Union, Optional

from brian2 import check_units, units, PoissonGroup, BrianObject, Equations

from meanfield.parameters import PP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import lazyproperty


class MFPoissonPopulation(MFPopulation):

    arguments = MappingProxyType({})

    defaults = MappingProxyType({})

    @check_units(rate=units.Hz)
    def __init__(self, n: int, rate: units.Hz, **kwargs):
        super().__init__(n, {}, **kwargs)
        self.rate = rate

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    @property
    @check_units(result=units.Hz)
    def rate_prediction(self) -> units.Hz:
        return self.rate

    # Simulation

    @lazyproperty
    def brian2(self) -> BrianObject:
        return PoissonGroup(self.n, self.rate, name=self.name)

    @property
    def brian2_model(self) -> Optional[Equations]:
        return None
