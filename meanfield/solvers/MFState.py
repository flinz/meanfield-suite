import numpy as np
from brian2 import Quantity

from meanfield.solvers import OutOfBoundsError


class MFState(object):

    def __init__(self, constraints, dependent_functions=None, bounds_check=False, bounds_correct=False, bounds_escalate=True, bounds_interrupt=False):

        if dependent_functions is None:
            dependent_functions = []

        self.constraints = constraints
        self.dependent_functions = dependent_functions
        self.bounds_check = bounds_check
        self.bounds_correct = bounds_correct
        self.bounds_escalate = bounds_escalate
        self.bounds_interrupt = bounds_interrupt

    def __getitem__(self, key):
        """Dictionary-like access to state values."""
        try:
            idx = self.names.index(key)
        except ValueError:
            raise KeyError
        return self.state[idx]

    @property
    def names(self):
        return [c.name for c in self.constraints]

    @property
    def state(self):
        return [c.free for c in self.constraints]

    @state.setter
    def state(self, value):
        assert len(value) == len(self.constraints)
        for i, c in enumerate(self.constraints):
            dim = c.free.dimensions
            c.free = Quantity(value[i], dim=dim)
        for f in self.dependent_functions:
            f()

    @property
    def error(self):
        """Error vector of state."""
        return [c.error for c in self.constraints]

    @property
    def error_checked(self):
        """Bounds checked error."""
        # do bounds checking
        if self.bounds_check:
            bstate, is_out, err_states = self.checkbounds()
            if is_out:
                if self.bounds_interrupt:
                    raise OutOfBoundsError
                if self.bounds_correct:
                    self.state = bstate
                if self.bounds_escalate:
                    return err_states * 1e5  # return more detailed error to solver, i.e. which direction was wrong.

        return self.error

    def checkbounds(self, hard=True):
        """
        Check bounds and set variables to bounds
        -   bounds: array of (low,up)
        -   hard: whether set to hard bounds.
            False to reflect at bounds. Watch out, causes errors for too large reflections
        -   change: True to set state to the bounded values. False for only detection.
        """
        bstate = np.array(self.state).copy()
        err_states = np.zeros(len(self.state))

        for i, c in enumerate(self.constraints):

            s = self.state[i]

            # lower bound check
            if s < c.bound_down:
                if not c.bound_hard:  # reflect at bound
                    bstate[i] = 2 * c.bound_down - s
                else:  # set hard bound
                    bstate[i] = c.bound_down
                err_states[i] = 1

            # same for upper bound
            elif s > c.bound_up:
                if not c.bound_hard:
                    bstate[i] = 2 * c.bound_up - s
                else:
                    bstate[i] = c.bound_up
                err_states[i] = 1
            else:
                bstate[i] = self.state[i]
        return bstate, any(err_states), err_states

    def __repr__(self):
        return "MFState<" + ", ".join(["%s: %.3f" % (c.name, c.free) for c in self.constraints]) + ">"

    def __call__(self, value, fun=lambda x: x):
        self.state = value
        #print(value)
        error = self.error_checked
        return fun(error)
