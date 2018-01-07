import matplotlib.pyplot as plt
import pytest
from brian2 import PopulationRateMonitor, units
from brian2.units import *

from meanfield.solvers.MFSolver import MFSolver
from meanfield.zoo.brunel_wang_2001_stp import one_subpopulation
from meanfield.utils import reset_brian2


class TestMF2(object):

    def setup_method(self):
        self.system = one_subpopulation()

    def test_MFSolver(self):
        """Solve for firing rates & voltages with explicit function"""
        up, down, inter = self.system.populations

        up.rate = 1 * Hz
        down.rate = 4 * Hz
        inter.rate = 12 * Hz

        solver = MFSolver.rates_voltages(self.system, solver='simplex', maxiter=1)
        sol = solver.run()
        print(sol)

    def test_graph(self):
        #print(self.system.graph().render('test.png', view=True))
        pass

    @pytest.mark.skip(reason='plotting')
    def test_sim(self):
        reset_brian2()
        up, down, inter = self.system.populations

        # MFState < Eup - rate: 1.055, Eup - v_mean: -0.054, Edown - rate: 4.285, Edown - v_mean: -0.053, I - rate: 11.463, I - v_mean: -0.053 >

        rate_up = PopulationRateMonitor(up.brian2)
        rate_down = PopulationRateMonitor(down.brian2)
        rate_inter = PopulationRateMonitor(inter.brian2)

        net = self.system.collect_brian2_network(rate_up, rate_down, rate_inter)
        net.run(3 * units.second)

        plt.plot(rate_up.t / ms, rate_up.smooth_rate(width=200 * ms) / Hz, label='up')
        plt.plot(rate_down.t / ms, rate_down.smooth_rate(width=200 * ms) / Hz, label='down')
        plt.plot(rate_inter.t / ms, rate_inter.smooth_rate(width=200 * ms) / Hz, label='inter')
        plt.legend()
        plt.show()



