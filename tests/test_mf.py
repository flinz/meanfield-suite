from functools import partial

import numpy as np
import pytest
from brian2.units import *

from meanfield.meanfield_classes import setup_brunel99, setup_EI
from meanfield.solvers.MFConstraint import MFConstraint
from meanfield.solvers.MFSolver import MFSolver, MFSolverRatesVoltages
from meanfield.solvers.MFState import MFState

@pytest.mark.skip(reason="need nmda impl")
class TestMF(object):

    def setup_method(self):
        self.system = setup_brunel99()

    def test_taus(self):
        """Effective timeconstants equal old implementation"""
        taus = [p.tau_eff / ms for p in self.system.pops]
        #taus_brunel = [11.309393346834508, 10.259016117146679, 5.2504978292166333] # NMDA
        taus_brunel = [4.334601,  5.106385,  2.545614]
        np.testing.assert_array_almost_equal(taus, taus_brunel)
        #assert np.allclose(taus, taus_brunel), "Values not equal: {} != {}".format(taus, taus_brunel)

    def test_mus(self):
        """Mean input predictions equal old implementation"""
        mus = [p.mu for p in self.system.pops]
        print(mus)
        #mus_brunel = [18.868912317468116, 16.462020986464438, 16.524683021279241] # NMDA
        #mus_brunel = [41.174904,  36.042538,  36.576346]
        mus_brunel = [0.041175,  0.036043,  0.036576]
        np.testing.assert_array_almost_equal(mus, mus_brunel)

    def test_sigma_square(self):
        """Sigma square predictions equal old implementation"""
        sigma_square = [p.sigma_square / (mV) ** 2 for p in self.system.pops]
        #sigma_square_brunel = np.array([4.8869461761143898, 5.1557159873625888, 10.003849121175195]) # NMDA
        sigma_square_brunel = [1.873041,  2.566238,  4.850195]
        np.testing.assert_array_almost_equal(sigma_square, sigma_square_brunel)

    def test_firing_rate(self):
        """Mean firing rate predictions equal old implementation"""
        rate_pred = np.array([p.rate_prediction for p in self.system.pops])
        rate_pred_brunel = np.array([0.015207823717059475, 0.0013208593687499856, 0.0068435294603920544])
        diff = (rate_pred - 1e3 * rate_pred_brunel)
        assert all(diff < 1e-8), "Values not equal: %s" % diff

    def test_mean_voltage(self):
        """Mean voltage predictions equal old implementation"""
        v_mean_prediction = np.array([p.v_mean_prediction for p in self.system.pops])
        v_mean_prediction_brunel = np.array([4.8869461761143898, 5.1557159873625888, 10.003849121175195])
        diff = (v_mean_prediction - v_mean_prediction_brunel)
        assert all(diff < 1e-8), "Values not equal: %s" % diff

    def test_MFConstraint(self):
        """Create constraint, and basic properties"""
        pop = self.system.pops[0]
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
                          ) for p in self.system.pops
                      ] + [
                          MFConstraint(
                              "%s-%s" % (p.name, "v_mean"),
                              partial(lambda x: x.v_mean, p),
                              partial(lambda x, val: setattr(x, "v_mean", val), p),
                              partial(lambda x: x.v_mean-x.v_mean_prediction, p),
                              -80. * mV, 50. * mV
                          ) for p in self.system.pops
                      ]

        state = MFState(constraints)
        error = [c.error_fun() for c in state.constraints]
        state_ = [c.free for c in state.constraints]
        with pytest.raises(AttributeError):
            state.error = 1.
        assert state.error == error
        np.testing.assert_array_almost_equal(state.state, state_)
        print(state)
        return state

    def test_MFSolver_RatesVoltages(self):
        """Solve for firing rates & voltages with specialized subclass"""
        solver = MFSolverRatesVoltages(self.system)
        r1 = solver.run()

        # take old implementation and compare
        state = self.testMFState()
        solver = MFSolver(state, solver='gradient')
        r2 = solver.run()

        for key in r1.names:
            np.testing.assert_almost_equal(r1[key], r2[key])

    def test_MFState_DictionarylikeAccess(self):
        """Dictionarylike access to MFState"""
        state = self.testMFState()
        for i in range(len(state.names)):
            name = state.names[i]
            val = state.state[i]
            assert val == state[name]

    def test_MFSolver(self):
        """Solve for firing rates & voltages with explicit function"""
        state = self.testMFState()
        solver = MFSolver(state)
        solver.run()

    def test_MFSolver_RatesVoltagesNoNmda(self):
        """Solve for firing rates only with specialized subclass"""

        system1 = setup_EI(has_nmda=False)
        solver = MFSolverRatesVoltages(system1, maxiter=100)
        r1 = solver.run()

        system2 = setup_EI(has_nmda=False)
        solver = MFSolverRatesVoltages(system2, force_nmda=True, maxiter=100)
        r2 = solver.run()

        # test rates
        for key in r1.names:
            np.testing.assert_almost_equal(r1[key], r2[key])

        # test voltages. these are set implicitly in r1 by the dependent_function
        for key in [("E", "E-v_mean"), ("I", "I-v_mean")]:
            np.testing.assert_almost_equal(system1[key[0]].v_mean_prediction, r2[key[1]], 5)
            np.testing.assert_almost_equal(system1[key[0]].v_mean, r2[key[1]], 5)

    def test_GradientMinimization_RatesVoltages(self):
        """Solve for firing rates & voltages with specialized subclass"""

        system = setup_EI()

        solver = MFSolverRatesVoltages(system, solver="gradient")
        sol1 = solver.run()

        solver = MFSolverRatesVoltages(system)
        sol2 = solver.run()

        np.testing.assert_array_almost_equal(sol1.state, sol2.state, 5)

    def test_ConductanceMinimization(self):
        """Solve for NMDA conductances with constrained firing rates"""

        system = setup_EI()

        constraints = [
                          MFConstraint(
                              "%s-%s" % (p.name, "gNMDA"),
                              partial(lambda x: x.sources[1].g_base, p),
                              partial(lambda x, val: setattr(x.sources[1], "g_base", val), p),
                              partial(lambda x: x.rate - x.rate_prediction, p),
                              0. * nsiemens, 500. * nsiemens
                          ) for p in system.pops
                      ] + [
                          MFConstraint(
                              "%s-%s" % (p.name, "v_mean"),
                              partial(lambda x: x.v_mean, p),
                              partial(lambda x, val: setattr(x, "v_mean", val), p),
                              partial(lambda x: x.v_mean - x.v_mean_prediction, p),
                              -80. * mV, -50. * mV
                          ) for p in system.pops
                      ]

        state = MFState(constraints, bounds_check=True)
        solver = MFSolver(state, solver='hybr')
        solver.run()

        for p in system.pops:
            assert p.rate_prediction.has_same_dimensions(p.rate)
            np.testing.assert_almost_equal(np.array(p.rate_prediction), np.array(p.rate))

    def test_ConductanceMinimizationRatio(self):
        """Solve for NMDA & Gaba conductances with constrained firing rates & EI fixed ratio"""

        print('start', flush=True)
        system = setup_EI()
        ratio = 4.

        print('setup', flush=True)

        def e_setter(p, val):
            setattr(p.sources[1], "g_base", val)
            setattr(p.sources[2], "g_base", ratio * val)

        constraints = [
                          MFConstraint(
                              "%s-%s" % (p.name, "gNMDA"),
                              partial(lambda x: x.sources[1].g_base, p),
                              partial(e_setter, p),
                              partial(lambda x: x.rate - x.rate_prediction, p),
                              0. * nsiemens, 500. * nsiemens
                          ) for p in system.pops
                      ] + [
                          MFConstraint(
                              "%s-%s" % (p.name, "v_mean"),
                              partial(lambda x: x.v_mean, p),
                              partial(lambda x, val: setattr(x, "v_mean", val), p),
                              partial(lambda x: x.v_mean-x.v_mean_prediction, p),
                              -80. * mV, 50. * mV
                          ) for p in system.pops
                      ]

        print('state', flush=True)
        state = MFState(constraints)
        print('solve', flush=True)

        solver = MFSolver(state, solver='gradient')
        solver.run()

        for p in system.pops:
            assert p.rate_prediction.has_same_dimensions(p.rate)
            np.testing.assert_almost_equal(np.array(p.rate_prediction), np.array(p.rate))
            assert p.sources[1].g_base == p.sources[1].g_base

