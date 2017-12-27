from types import MappingProxyType
from typing import Dict, Union

from brian2 import units, check_units, PoissonInput

from meanfield.inputs.MFInput import MFInput
from meanfield.utils import lazyproperty
from meanfield.parameters import IP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.parameters import Connection
from meanfield.parameters.Connection import ConnectionStrategy


class MFStaticInput(MFInput):

    arguments = MappingProxyType({})

    defaults = MappingProxyType({})

    @check_units(rate=units.hertz, n=1)
    def __init__(self, n: int, rate: units.hertz, target: MFPopulation, parameters: Union[Dict, MFParams], **kwargs):
        super().__init__(target, parameters, **kwargs, add_as_input=False)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

        self.rate = rate
        self.n = n

    # Theory

    @check_units(result=1)
    def g_dyn(self):
        return self.rate * self.n * self.parameters[IP.TAU]

    # Simulation

    @lazyproperty
    def brian2(self):
        return PoissonInput(self.target.brian2, self.post_variable_name, self.n, self.rate, 1)

