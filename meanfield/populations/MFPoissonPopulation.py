from types import MappingProxyType
from typing import Union, Optional

from brian2 import check_units, units, PoissonGroup, BrianObject, Equations

from meanfield.parameters import PP
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import lazyproperty
from meanfield.parameters.MFParams import MFParams


class MFPoissonPopulation(MFPopulation):

    arguments = MappingProxyType({
        PP.GM: units.siemens,
        PP.VRES: units.volt,
        PP.TAU_RP: units.second
    })

    defaults = MappingProxyType({})

    @check_units(rate=units.Hz)
    def __init__(self, n: int, rate: units.Hz, parameters: Union[dict, MFParams], **kwargs):
        super().__init__(n, parameters, **kwargs)
        self.rate = rate

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory

    @property
    @check_units(result=units.Hz)
    def rate_prediction(self) -> units.Hz:
        return self.rate

    @property
    @check_units(result=units.volt)
    def v_mean_prediction(self) -> units.volt:
        raise NotImplementedError

    # Population

    @lazyproperty
    def brian2(self) -> BrianObject:
        return PoissonGroup(self.n, self.rate, name=self.name)

    def brian2_model(self) -> Optional[Equations]:
        return None
