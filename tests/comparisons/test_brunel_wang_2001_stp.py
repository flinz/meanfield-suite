from functools import partial

import numpy as np
import pytest
from brian2 import PopulationRateMonitor, units
from brian2.units import *

from meanfield.solvers.MFConstraint import MFConstraint
from meanfield.solvers.MFSolver import MFSolver, MFSolverRatesVoltages
from meanfield.solvers.MFState import MFState
from meanfield.zoo.brunel_wang_2001_stp_depression import one_subpopulation, no_subpopulation
from meanfield.parameters import IP
import matplotlib.pyplot as plt

class TestMF2(object):

    def setup_method(self):
        self.system = one_subpopulation()


    def test_MFSolver(self):
        """Solve for firing rates & voltages with explicit function"""
        up, down, inter = self.system.populations

        up.rate = 1 * Hz
        down.rate = 1 * Hz
        inter.rate = 4 * Hz

        solver = MFSolverRatesVoltages(self.system, solver='simplex', maxiter=5)
        sol = solver.run()
        print(sol)

    def test_graph(self):
        #print(self.system.graph().render('test.png', view=True))
        pass

    def test_sim(self):
        up, down, inter = self.system.populations

        # MFState<Eup-rate: 0.887, Edown-rate: 4.401, I-rate: 11.466, Eup-v_mean: -0.054, Edown-v_mean: -0.053, I-v_mean: -0.053>
        # MFState<Eup-rate: 10.724, Edown-rate: 3.154, I-rate: 11.041, Eup-v_mean: -0.053, Edown-v_mean: -0.053, I-v_mean: -0.053>

        rate_up = PopulationRateMonitor(up.brian2)
        rate_down = PopulationRateMonitor(down.brian2)
        rate_inter = PopulationRateMonitor(inter.brian2)

        net = self.system.collect_brian2_network(rate_up, rate_down, rate_inter)
        net.run(1 * units.second)

        plt.plot(rate_up.t / ms, rate_up.smooth_rate(width=25 * ms) / Hz, label='up')
        plt.plot(rate_down.t / ms, rate_down.smooth_rate(width=25 * ms) / Hz, label='down')
        plt.plot(rate_inter.t / ms, rate_inter.smooth_rate(width=25 * ms) / Hz, label='inter')
        plt.legend()
        plt.show()



