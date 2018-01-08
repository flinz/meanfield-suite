import logging
import signal
import sys
from functools import partial
from timeit import default_timer as timer

import matplotlib.pylab as plt
import numpy as np
from brian2 import units
from brian2.units import fundamentalunits
from scipy.optimize import root, minimize

from meanfield.inputs.MFLinearNMDAInput import MFLinearNMDAInput
from meanfield.populations.MFPoissonPopulation import MFPoissonPopulation
from meanfield.solvers.MFConstraint import MFConstraint
from meanfield.solvers.MFState import MFState
from meanfield.solvers.algorithms import custom_gradient_solver

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("solver")

class MFSolver(object):

    @staticmethod
    def rates_voltages(system, force_nmda=False, *args, **kwargs):
        # create constraints on the firing rates and mean voltages

        constraints = []
        functions = []

        for p in system.populations:

            if not isinstance(p, MFPoissonPopulation):
                constraints.append(
                    MFConstraint(
                        "%s-%s" % (p.name, "rate"),
                        free_get=partial(lambda x: x.rate, p),
                        free_set=partial(lambda x, val: setattr(x, "rate", val), p),
                        error_fun=partial(lambda x: x.rate - x.rate_prediction, p),
                        bound_down=0. * units.Hz,
                        bound_up=1000. * units.Hz
                    )
                )

                if hasattr(p, 'v_mean'):

                    has_nmda = any(isinstance(i, MFLinearNMDAInput) for i in p.inputs)

                    if has_nmda or force_nmda:
                        print("Population %s has NMDA -> solving for voltages" % p.name)
                        constraints.append(
                            MFConstraint(
                                "%s-%s" % (p.name, "v_mean"),
                                partial(lambda x: x.v_mean, p),
                                partial(lambda x, val: setattr(x, "v_mean", val), p),
                                partial(lambda x: x.v_mean - x.v_mean_prediction, p),
                                -80. * units.mV, 50. * units.mV
                            )
                        )
                    else:
                        functions.append(
                            partial(lambda x: setattr(x, "v_mean", x.v_mean_prediction), p),
                        )

        state = MFState(constraints, dependent_functions=functions)
        return MFSolver(state, *args, **kwargs)

    def __init__(self, mfstate, maxiter=5, solver="hybr", crit=1e-4, tol=1e-12, fail_val=None):

        self.mfstate = mfstate
        self.solver = solver
        self.maxiter = maxiter
        self.crit = crit
        self.tol = tol
        self.it = 0
        self.error = 1e10
        self.state = "initialized"
        self.fail_val = fail_val

    def run(self, state_0=None, noise_percent=.1):

        # set solver state to the passed state if none was explicitly given
        if not state_0:
            state_0 = self.mfstate.state

        self.it = 0
        abs_err = 1e10
        min_sol = None
        min_abs_err = 1e10
        crit = self.crit
        tol = self.tol

        logger.info(f"initializing minimization with {self.solver}")

        # allow interruption of minimization loop, did not work otherwise.
        signal_handler = lambda signal, frame: sys.exit(0)
        signal.signal(signal.SIGINT, signal_handler)

        # do one pass to check unit, then disable unit checking to speed up optimization
        self.mfstate(state_0)
        fundamentalunits.unit_checking = False

        start = timer()
        steps = []

        # loop over tries to solve the system
        while abs_err > crit:

            self.it += 1

            # interrupt upon reaching the maxiter.
            if self.it > self.maxiter:
                logger.info('maximum iterations reached')
                self.state = "EXCEEDED"
                return self.finalize(min_sol)

            # get stochastic initial state
            up_dist = [c.bound_up - state_0[i] for i, c in enumerate(self.mfstate.constraints)]
            down_dist = [state_0[i] - c.bound_down for i, c in enumerate(self.mfstate.constraints)]
            max_dist = [min(up_dist[i], down_dist[i]) for i in range(len(self.mfstate.constraints))]

            p_0 = state_0 + (2. * np.random.rand(len(state_0)) - 1.) * np.array(max_dist) * noise_percent

            bounds = [[c.bound_down, c.bound_up] for c in self.mfstate.constraints]

            plotting = False

            # solve
            if self.solver == "gradient":  # own implementation
                sol = custom_gradient_solver(self.mfstate, p_0)
                abs_err = max(np.abs(sol.fun))

            elif self.solver == 'simplex':
                sq = lambda y: np.sum(np.array(y) ** 2)

                def bounded_f(x):
                    v = self.mfstate(x, fun=sq)

                    overbound = max(max(0, c.bound_down - xi, xi - c.bound_up) for xi, c in zip(x, self.mfstate.constraints))
                    #print(overbound)

                    return v + np.expm1(overbound)

                sol = minimize(bounded_f, p_0, method='Nelder-Mead')
                abs_err = np.sqrt(sol.fun)
                logger.debug(sol)

                if not sol.success:
                    logger.warning('optimizer failed')


            elif self.solver == "mse":
                sq = lambda y: np.sum(np.array(y) ** 2)

                def f(x):
                    v = self.mfstate(x, fun=sq)
                    return v

                sol = minimize(f, p_0, bounds=bounds, method='L-BFGS-B')
                #, options={
                #'disp': None,
                #'maxls': 20,
                #'iprint': -1,
                #'gtol': 1e-05,
                #    'eps': 0.00001,
                #'maxiter': 15000,
                #'ftol': 2.220446049250313e-09,
                #'maxcor': 10,
                #'maxfun': 15000
                #})
                logger.debug(sol)

                if not sol.success:
                    logger.warning('optimizer failed')

                abs_err = np.sqrt(sol.fun)

            else:  # scipy solvers
                sol = root(self.mfstate, p_0, jac=None, method=self.solver, tol=tol)

                abs_err = max(abs(sol.fun))


            # calculate the abs err, store minimal solution
            if abs_err < min_abs_err:
                min_abs_err = abs_err
                min_sol = sol

            if self.it > 1 and self.it % 50 == 0 and self.it < self.maxiter:
                sys.stdout.write('\n ')
            sys.stdout.flush()

            end = timer()
            steps.append(end - start)
            start = end

        fundamentalunits.unit_checking = True
        logger.info('finished successfully')
        logger.info(f'solver took {np.sum(steps).round(3)}s')
        self.state = "SUCCESS"
        return self.finalize(sol)

    def finalize(self, sol):
        self.mfstate.state = sol.x

        # set the state value to the fail_val if we did not converge
        if (self.fail_val is not None) and self.state != 'SUCCESS':
            self.mfstate.state = sol.x * self.fail_val

        self.error = sol.fun
        return self.mfstate
