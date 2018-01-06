from types import MappingProxyType
from typing import Dict, Union

from brian2 import Equations, Synapses, check_units, BrianObject

from meanfield.inputs.MFLinearNMDAInput import MFLinearNMDAInput
from meanfield.utils import lazyproperty
from meanfield.parameters import IP
from meanfield.parameters.MFParams import MFParams
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.parameters import Connection
from meanfield.parameters.Connection import ConnectionStrategy
from scipy import optimize
import numpy as np

interpy = [0., 0.23969905, 0.38048725, 0.47689264, 0.56519249, 0.61862329, 0.65427461, 0.69026088, 0.72211861, 0.74889075, 0.76569837, 0.78312683, 0.79899099, 0.81322308, 0.82531825, 0.83315572, 0.843196, 0.84944546, 0.85572625, 0.86329107, 0.86957444, 0.87651456, 0.88072668, 0.88622306, 0.89035198, 0.89594234, 0.89912681, 0.90223962, 0.90720305, 0.90878449, 0.91204239, 0.91506795, 0.91678712, 0.91852342, 0.92235024, 0.92395159]

interpx = [0., 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175]

#fit_syn = Synapse(interpx,interpy)


class Synapse(object):

    fits = {
        0: lambda x, a: a * x,
        1: lambda x, a, b, c, d, e: e + b / (c + a * np.exp(-x * d)),
        2: lambda x, a, b, c, d: a * x ** b / (c + x ** b) + d,
    }

    fun = None

    def __init__(self, fit_data_x, fit_data_y, fit_type=2):
        self.fit_type = fit_type
        self.fit_data_x = np.array(fit_data_x)
        self.fit_data_y = np.array(fit_data_y)
        self.fit_p_opt = None

    def do_fit(self):
        print("Fitting function to data, fit type %i".format(self.fit_type))
        fun = Synapse.fits[self.fit_type]
        # n_args = len(inspect.getargspec(fun).args)
        # p_0 = [.1]*n_args
        popt, pcov = optimize.curve_fit(fun, self.fit_data_x, self.fit_data_y)
        self.fit_p_opt = popt

    def __call__(self, x):
        if self.fun == None:
            self.fun = self.make_fun()
        return self.fun(x)

    def make_fun(self):
        if self.fit_p_opt == None:
            self.do_fit()
        return lambda x: max(0, Synapse.fits[self.fit_type](x, *self.fit_p_opt))



class MFNonLinearNMDAInput(MFLinearNMDAInput):

    arguments = MappingProxyType({
        IP.TAU_NMDA: 1.,
        IP.TAU_NMDA_RISE: 1.,
        IP.ALPHA: 1.
    })

    defaults = MappingProxyType({})

    def __init__(self, origin: MFPopulation, target: MFPopulation, parameters: Union[Dict, MFParams], **kwargs):
        super().__init__(origin, target, parameters, **kwargs)

        self.parameters.fill(self.defaults)
        self.parameters.verify(self.arguments)

    # Theory





    # Simulation

    @property
    def post_nonlinear_name(self):
        return 'x_' + self.ref

    @lazyproperty
    def brian2(self) -> BrianObject:
        model = Equations(
            '''
            w : 1
            s_post = w * s : 1 (summed)
            ds / dt = - s / tau_decay + alpha * x * (1 - s) : 1 (clock-driven)
            dx / dt = - x / tau_rise : 1 (clock-driven)
            ''',
            s_post=self.post_variable_name + '_post',
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
            I = g * (v - ve) / (1 + gamma * exp(- beta * v) ) * s_post : amp
            s_post: 1
            ''',
            s_post=self.post_variable_name + '_post',
            I=self.current_name,
            g=self[IP.GM],
            ve=self[IP.VE],
            beta=self[IP.BETA]
        )
