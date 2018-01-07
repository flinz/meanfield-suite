from functools import partial

import numpy as np
import pytest
from brian2.units import *

from meanfield.solvers.MFConstraint import MFConstraint
from meanfield.solvers.MFSolver import MFSolver
from meanfield.solvers.MFState import MFState
from meanfield.zoo.brunel_wang_2001_stp import one_subpopulation, no_subpopulation
from meanfield.parameters import IP


class TestMF(object):

    def setup_method(self):
        self.system = one_subpopulation()

    def test_taus(self):
        """Effective timeconstants equal old implementation"""
        taus = [p.tau_eff / ms for p in self.system.populations]
        # taus_brunel = [11.309393346834508, 10.259016117146679, 5.2504978292166333]  # NMDA STP 3TS
        taus_brunel = [11.343448458849716, 10.277459899195716, 5.260438755634331]  # NMDA STP
        np.testing.assert_array_almost_equal(taus, taus_brunel)

    def test_mus(self):
        """Mean input predictions equal old implementation"""
        mus = [p.mu / mV for p in self.system.populations]
        # mus_brunel = [18.868912317468116, 16.462020986464438, 16.524683021279241]  # NMDA STP 3TS
        mus_brunel = [18.95643395445861, 16.512878547524384, 16.5783612995361]  # NMDA STP
        np.testing.assert_array_almost_equal(mus, mus_brunel)

    def test_sigma_square(self):
        """Sigma square predictions equal old implementation"""
        sigma_square = [p.sigma_square / (mV) ** 2 for p in self.system.populations]
        # sigma_square_brunel = [4.8869461761143898, 5.1557159873625888, 10.003849121175195]  # NMDA STP 3TS
        sigma_square_brunel = [4.901661863714598, 5.164984995315381, 10.022789711428175]  # NMDA STP
        np.testing.assert_array_almost_equal(sigma_square, sigma_square_brunel)

    def test_firing_rate(self):
        """Mean firing rate predictions equal old implementation"""
        rate_pred = [p.rate_prediction / Hz for p in self.system.populations]
        # rate_pred_brunel = [15.207823717059475, 1.3208593687499856, 6.8435294603920544]  # NMDA STP 3TS
        rate_pred_brunel = [15.978576065910952, 1.4435632224458281, 7.249972010396413]  # NMDA STP
        np.testing.assert_array_almost_equal(rate_pred, rate_pred_brunel)

    def test_MFConstraint(self):
        """Create constraint, and basic properties"""
        pop = self.system.populations[0]
        c = MFConstraint(
            'aa',
            partial(lambda x: x.rate, pop),
            partial(lambda x, val: setattr(x, "rate", val), pop),
            partial(lambda x: x.rate - x.rate_prediction, pop),
            0.,
            750.
        )
        c.free = 111. * Hz
        assert pop.rate == c.free == 111. * Hz
        assert c.error > 0.

    def test_MFState(self):
        """Create state from constraints, and check basic properties"""
        constraints = [
                          MFConstraint(
                              "%s-%s" % (p.name, "rate"),
                              partial(lambda x: x.rate, p),
                              partial(lambda x, val: setattr(x, "rate", val), p),
                              partial(lambda x: x.rate - x.rate_prediction, p),
                              0. * Hz, 750. * Hz
                          ) for p in self.system.populations
                      ] + [
                          MFConstraint(
                              "%s-%s" % (p.name, "v_mean"),
                              partial(lambda x: x.v_mean, p),
                              partial(lambda x, val: setattr(x, "v_mean", val), p),
                              partial(lambda x: x.v_mean-x.v_mean_prediction, p),
                              -80. * mV, 50. * mV
                          ) for p in self.system.populations
                      ]

        state = MFState(constraints)
        error = [c.error_fun() for c in state.constraints]
        state_ = [c.free for c in state.constraints]
        with pytest.raises(AttributeError):
            state.error = 1.
        assert state.error == error
        np.testing.assert_array_almost_equal(state.state, state_)
        return state

    def test_MFSolver_RatesVoltages(self):
        """Solve for firing rates & voltages with specialized subclass"""

        solver = MFSolver.rates_voltages(self.system)
        r1 = solver.run()

        # take old implementation and compare
        state = self.test_MFState()
        solver = MFSolver(state)
        r2 = solver.run()

        for key in [c.name for c in r1.constraints]:
            assert r1[key] == r2[key]

    def test_MFState_DictionarylikeAccess(self):
        """Dictionarylike access to MFState"""
        state = self.test_MFState()
        for i in range(len(state.constraints)):
            name = state.constraints[i].name
            val = state.state[i]
            assert val == state[name]

    def test_MFSolver(self):
        """Solve for firing rates & voltages with explicit function"""
        state = self.test_MFState()
        solver = MFSolver(state)
        solver.run()

    def test_GradientMinimization_RatesVoltages(self):
        """Solve for firing rates & voltages with specialized subclass"""

        system = no_subpopulation()

        solver = MFSolver.rates_voltages(system, solver="mse", maxiter=1)
        sol1 = solver.run()

        solver = MFSolver.rates_voltages(system, solver="gradient", maxiter=1)
        sol2 = solver.run()

        np.testing.assert_array_almost_equal(sol1.state, sol2.state, 5)

    def test_ConductanceMinimization(self):
        """Solve for NMDA conductances with constrained firing rates"""

        system = no_subpopulation()

        constraints = [
                          MFConstraint(
                              "%s-%s" % (p.name, "gNMDA"),
                              partial(lambda x: x.inputs[0][IP.GM], p),
                              partial(lambda x, val: x.inputs[0].__setitem__(IP.GM, val), p),
                              partial(lambda x: x.rate - x.rate_prediction, p),
                              0. * nsiemens, 500. * nsiemens
                          ) for p in system.populations
                      ] + [
                          MFConstraint(
                              "%s-%s" % (p.name, "v_mean"),
                              partial(lambda x: x.v_mean, p),
                              partial(lambda x, val: setattr(x, "v_mean", val), p),
                              partial(lambda x: x.v_mean - x.v_mean_prediction, p),
                              -80. * mV, -50. * mV
                          ) for p in system.populations
                      ]

        state = MFState(constraints, bounds_check=True)
        solver = MFSolver(state, solver='hybr')
        solver.run()

        for p in system.populations:
            assert p.rate_prediction.has_same_dimensions(p.rate)
            np.testing.assert_almost_equal(np.array(p.rate_prediction), np.array(p.rate))

    def test_ConductanceMinimizationRatio(self):
        """Solve for NMDA & Gaba conductances with constrained firing rates & EI fixed ratio"""

        system = no_subpopulation()
        ratio = 4.

        def e_setter(p, val):
            p.inputs[0].__setitem__(IP.GM, val)
            p.inputs[1].__setitem__(IP.GM, ratio * val)

        constraints = [
                          MFConstraint(
                              "%s-%s" % (print(p.inputs), "gNMDA"),
                              partial(lambda x: x.inputs[0][IP.GM], p),
                              partial(e_setter, p),
                              partial(lambda x: x.rate - x.rate_prediction, p),
                              0. * nsiemens, 500. * nsiemens
                          ) for p in system.populations
                      ] + [
                          MFConstraint(
                              "%s-%s" % (p.name, "v_mean"),
                              partial(lambda x: x.v_mean, p),
                              partial(lambda x, val: setattr(x, "v_mean", val), p),
                              partial(lambda x: x.v_mean-x.v_mean_prediction, p),
                              -80. * mV, 50. * mV
                          ) for p in system.populations
                      ]

        state = MFState(constraints)
        solver = MFSolver(state, maxiter=1)
        solver.run()

        for p in system.populations:
            assert p.rate_prediction.has_same_dimensions(p.rate)
            np.testing.assert_almost_equal(np.array(p.rate_prediction), np.array(p.rate))

