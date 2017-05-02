import numpy as np

from MFLinearPop import MFLinearPop
from params import NP


class MFNonLinearPop(MFLinearPop):
    """pop: similar neurons"""

    def __init__(self, name, n, params):
        super().__init__(name, n, params)

        defaults = {}
        expectations = {
            NP.GAMMA: 1,
            NP.BETA: 1,
        }

        self.params.fill(defaults)
        self.params.verify(expectations)

    @property
    def J(self):
        """Linearization factor for NMDA"""
        return 1 + self.params[NP.GAMMA] * np.exp(-self.params[NP.BETA] * self.v_mean)


