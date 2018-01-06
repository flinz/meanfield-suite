from types import MappingProxyType
from typing import Dict, Union

import numpy as np
from brian2 import Equations, Synapses, check_units, BrianObject, units
from scipy import optimize

from meanfield.inputs.MFLinearNMDAInput import MFLinearNMDAInput
from meanfield.parameters import IP
from meanfield.parameters.MFParameters import MFParameters
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import lazyproperty

# numerical fit
interpy = [0., 0.23969905, 0.38048725, 0.47689264, 0.56519249, 0.61862329, 0.65427461, 0.69026088, 0.72211861, 0.74889075, 0.76569837, 0.78312683, 0.79899099, 0.81322308, 0.82531825, 0.83315572, 0.843196, 0.84944546, 0.85572625, 0.86329107, 0.86957444, 0.87651456, 0.88072668, 0.88622306, 0.89035198, 0.89594234, 0.89912681, 0.90223962, 0.90720305, 0.90878449, 0.91204239, 0.91506795, 0.91678712, 0.91852342, 0.92235024, 0.92395159]
interpx = [0., 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175]


class Synapse(object):

    fits = {
        0: lambda x, a: a * x,
        1: lambda x, a, b, c, d, e: e + b / (c + a * np.exp(-x * d)),
        2: lambda x, a, b, c, d: a * x ** b / (c + x ** b) + d,
    }

    def __init__(self, fit_data_x, fit_data_y, fit_type=2):

        print("Fitting function to data, fit type {}".format(fit_type))
        fun = Synapse.fits[fit_type]

        popt, pcov = optimize.curve_fit(fun, fit_data_x, fit_data_y)
        self.fit = lambda x: fun(x, *popt)

    def __call__(self, x):
        return max(0, self.fit(x))


class MFNonLinearNMDAInput(MFLinearNMDAInput):

    arguments = MappingProxyType({
        IP.TAU_NMDA: units.second,
        IP.TAU_NMDA_RISE: units.second,
        IP.GAMMA: 1.,
        IP.BETA: 1.,
        IP.ALPHA: units.hertz,
        IP.MG: 1.,
    })

    defaults = MappingProxyType({})

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParameters], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

        self.synapse = Synapse(interpx, interpy)

    # Theory

    @check_units(result=1)
    def g_dyn(self):
        return self.connection.theory(self.origin.n) * self.synapse(self.origin.rate) * self[IP.W]

    # Simulation

    @property
    def post_nonlinear_name(self):
        return f'x_{self.ref}'

    @property
    def post_variable_name_tot(self):
        return f's_{self.ref}_tot'

    @lazyproperty
    def brian2(self) -> BrianObject:
        model = Equations(
            '''
            w : 1
            s_post = w * s : 1 (summed)
            ds / dt = - s / tau_decay + alpha * x * (1 - s) : 1 (clock-driven)
            dx / dt = - x / tau_rise : 1 (clock-driven)
            ''',
            s_post=f'{self.post_variable_name_tot}_post',
            s=self.post_variable_name,
            x=self.post_nonlinear_name,
            tau_decay=self[IP.TAU_NMDA],
            tau_rise=self[IP.TAU_NMDA_RISE],
            alpha=self[IP.ALPHA],
        )
        eqs_pre = f'''
        {self.post_nonlinear_name} += 1
        '''
        C = Synapses(
            self.origin.brian2,
            self.target.brian2,
            method='euler',
            model=model,
            on_pre=eqs_pre,
            name=self.ref,
        )
        C.connect()
        C.w[:] = 1
        return C

    @property
    def brian2_model(self) -> Equations:
        return Equations(
            '''
            I = g * (v - vrev) / (1 + mg * gamma * exp(- beta * v / mV) ) * s : amp
            s: 1
            ''',
            s=self.post_variable_name_tot,
            I=self.current_name,
            g=self[IP.GM],
            vrev=self[IP.VREV],
            beta=self[IP.BETA],
            gamma=self[IP.GAMMA],
            mg=self[IP.MG],
        )
