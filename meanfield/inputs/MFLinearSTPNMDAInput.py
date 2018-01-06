from types import MappingProxyType
from typing import Union, Dict

from meanfield.inputs.MFLinearNMDAInput import MFLinearNMDAInput
from meanfield.inputs.MFLinearSTPInput import MFLinearSTPInput
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPopulation import MFPopulation


class MFLinearSTPNMDAInput(MFLinearSTPInput, MFLinearNMDAInput):

    arguments = MappingProxyType({})

    defaults = MappingProxyType({})

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParams], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

