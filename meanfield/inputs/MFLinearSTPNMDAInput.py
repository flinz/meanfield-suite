from types import MappingProxyType
from typing import Union, Dict

from brian2 import units

from meanfield.inputs.MFLinearNMDAInput import MFLinearNMDAInput
from meanfield.inputs.MFLinearSTPInput import MFLinearSTPInput
from meanfield.parameters.MFParameters import MFParameters
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.parameters import IP


class MFLinearSTPNMDAInput(MFLinearSTPInput, MFLinearNMDAInput):

    arguments = MappingProxyType({
        IP.GM: units.siemens,
        IP.VREV: units.volt,
        IP.TAU: units.second,
        IP.W: 1,
        IP.TAU_F: units.second,
        IP.TAU_D: units.second,
        IP.U: 1,
        IP.BETA: 1,
        IP.GAMMA: 1,
    })

    defaults = MappingProxyType({
        IP.W: 1,
    })

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParameters], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

