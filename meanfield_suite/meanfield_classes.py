import numpy as np
from scipy.integrate import quad
from scipy.optimize import root
import pylab as pl
from math import erf
import sys
import signal
from functools import partial

params_standard = {
    "NMDA": {
        "gamma": 0.280112,
        "beta": 0.062,
    },
    "E": {
        "gamma": 0.280112,
        "beta": 0.062,
        "VE": 0.,
        "E_L": -70.,
        "VI": -70.,
        "V_th": -50.,
        "V_reset": -60.,
        "tau_AMPA": 2.,
        "t_ref": 2.,
        "C_m": 500.,
        "g_L": 25.,
        "Cext": 1000,
        "nu_ext": 0.0024,
        "gAMPA": 2.08,
        "gNMDA": 0.327,
        "gGABA": 1.25
    },
    "I": {
        "gamma": 0.280112,
        "beta": 0.062,
        "VE": 0.,
        "E_L": -70.,
        "VI": -70.,
        "V_th": -50.,
        "V_reset": -60.,
        "tau_AMPA": 2.,
        "t_ref": 1.,
        "C_m": 200.,
        "g_L": 20.,
        "Cext": 1000,
        "nu_ext": 0.0024,
        "gAMPA": 1.62,
        "gNMDA": 0.258,
        "gGABA": 0.973
    }
}


class OutOfBoundsError(Exception):
    """Raised when solver out of bounds."""
    pass


class MFSource(object):
    """Source: a synapse coupled to pops"""

    def __init__(self, name, pop):

        self.name = name
        self.g_dyn = lambda: 0.  # this should be given as a function describing the synaptic conductance
        self.is_nmda = False
        self.g_base = 0.  # [nS]
        self.E_rev = 0.   # [mV] excitatory by default
        self.noise_tau = 0.   # [mV] excitatory by default

        # link to pop
        pop.sources.append(self)
        self.pop = pop

    @property
    def conductance(self):
        tmp = self.g_dyn() * self.g_base
        if self.is_nmda:
            J_ = 1./self.pop.J
            return tmp * J_ * (
                1. + (1.-J_) * self.pop.params["beta"] * (self.pop.v_mean - self.E_rev)
            )
        return tmp

    @property
    def voltage_conductance(self):
        cond = self.g_dyn() * self.g_base
        tmp = cond * (self.E_rev - self.pop.params["E_L"])
        if self.is_nmda:
            J_ = 1./self.pop.J
            return J_ * (
                tmp + (1.-J_) * cond * self.pop.params["beta"] * (self.pop.v_mean - self.E_rev) * (self.pop.v_mean - self.pop.params["E_L"])
            )
        return tmp

    def print_sys(self):
        print("%s: %f" % (self, self.conductance))

    def __repr__(self):
        return "MFSource [%s] <%s, nmda: %s, E_rev: %.1f>" % (id(self),self.name, self.is_nmda, self.E_rev)


class MFpop(object):
    """pop: similar neurons"""

    def __init__(self, name, params):
        self.name = name
        self.sources = []
        self.params = params
        self.noise = None

        self.n = 1
        self.rate_ = 0.  # spikes/s
        self.v_mean = -60.        # pop mean voltage
        self.is_adapting = False

    @property
    def rate_ms(self):
        return self.rate_/1e3

    @property
    def rate_hz(self):
        return self.rate_

    @rate_hz.setter
    def rate_hz(self, value):
        self.rate_ = value

    @rate_ms.setter
    def rate_ms(self, value):
        self.rate_ = value*1e3

    @property
    def has_nmda(self):
        return any([s.is_nmda for s in self.sources])

    @property
    def J(self):
        """Linearization factor for NMDA"""
        return 1 + self.params["gamma"] * np.exp(-self.params["beta"]*self.v_mean)

    @property
    def total_cond(self):
        """Gm * SE in [1]"""
        return self.params["g_L"] + np.sum(s.conductance for s in self.sources)

    @property
    def tau_eff(self):
        return self.params["C_m"]/self.total_cond

    @property
    def mu(self):
        return 1./self.total_cond * np.sum(s.voltage_conductance for s in self.sources)

    @property
    def sigma_square(self):
        if not self.noise:
            return 0.
        return (self.noise.g_base / self.params["C_m"] * (self.v_mean - self.noise.E_rev))**2 * self.tau_eff * self.noise.g_dyn() * self.noise.noise_tau

    @property
    def rate_pred(self):
        return 1e3 * self.phi_firing_func()

    @property
    def v_mean_prediction(self):
        return self.params["E_L"] + self.mu - (self.params["V_th"]-self.params["V_reset"]) * self.rate_ms * self.tau_eff

    def phi_firing_func(self):

        sigma = np.sqrt(self.sigma_square)
        tau_eff = self.tau_eff

        def beta():
            return (self.params["V_reset"]-self.params["E_L"]-self.mu)/sigma

        def alpha():
            tmp = -0.5*self.noise.noise_tau/tau_eff \
                    + 1.03*np.sqrt(self.noise.noise_tau/tau_eff) \
                    + (-self.mu - self.params["E_L"] + self.params["V_th"])*(1. + (0.5 * self.noise.noise_tau / tau_eff))/sigma
            return tmp

        def integrand(x):
            if x < -10.:
                return np.exp(10.**2)*(1.+erf(10.))
            if x > 10.:
                return 0.
            return np.exp(x**2)*(1.+erf(x))
        return 1./(self.params["t_ref"] + tau_eff * np.sqrt(np.pi) * quad(integrand, beta(), alpha())[0])

    def __repr__(self):
        return "MFpop [%s] <%s (%i sources, n: %i, rate: %.4f, v_mean: %.4f)>" % (id(self), self.name, len(self.sources), self.n, self.rate_hz, self.v_mean)

    def print_sys(self, mf=False):
            print("\t%s - tau_eff: %.1fms, mu: %.4f, sig^2: %.4f, rate_pred: %.4f, v_mean_pred: %.4f" % (self, self.tau_eff, self.mu, self.sigma_square, self.rate_pred, self.v_mean_prediction))
            for s in self.sources:
                print("\t\t",
                s.print_sys())


class MFSystem(object):

    def __init__(self, name):
        self.name = name
        self.pops = []

    def __getitem__(self, key):
        """Dictionarylike access to state values."""
        names = [x.name for x in self.pops]
        try:
            idx = names.index(key)
        except ValueError:
            raise KeyError
        return self.pops[idx]

    def print_sys(self, mf=False):
        print("%s" % self)
        for p in self.pops:
            p.print_sys(mf)

    def __repr__(self):
        return "MFSystem <%s (%i pops)>" % (self.name, len(self.pops))


def stp_dgl_u(U, tau_f, tau_x, rate_ms):
    """ Differential equation equilibrium solution of short-term plastic synaptic input: facilitation variable u."""
    return float(U) * (1. + rate_ms * tau_f) * (1. + float(U) * rate_ms * tau_f)**(-1.)


def stp_dgl_x(U, tau_f, tau_x, rate_ms):
    """ Differential equation equilibrium solution of short-term plastic synaptic input: depression variable x."""
    return (1. + float(U) * rate_ms * tau_f) * (1. + float(U) * rate_ms * (tau_f + tau_x + rate_ms * tau_f * tau_x))**(-1.)


class Synapse(object):
    """
        Synapse with 3 time constants and possibly depression/facilitation.
    """
    fun = None

    def __init__(self, tau_syn_rise=1., tau_syn_d1=100., tau_syn_d2=100., balance=.5, U=1., tau_f=1000., tau_x=0.):

        # facilitation & depression parameters
        self.tau_f = tau_f
        self.tau_x = tau_x
        self.U = U

        # synaptic timeconstants parameters
        self.tau_syn_rise = tau_syn_rise
        self.tau_syn_d1 = tau_syn_d1
        self.tau_syn_d2 = tau_syn_d2
        self.balance = balance

    @property
    def taus(self):
        return np.array([
            self.balance * self.tau_syn_d1,
            (1.-self.balance) * self.tau_syn_d2,
            - self.balance * self.tau_syn_d1 * self.tau_syn_rise / (self.tau_syn_d1 + self.tau_syn_rise),
            - (1.-self.balance) * self.tau_syn_d2 * self.tau_syn_rise / (self.tau_syn_d2 + self.tau_syn_rise)
        ])

    def __call__(self, x):
        if self.fun is None:
            self.fun = self.make_fun()
        return self.fun(x)

    def stp_u(self, x):
        return stp_dgl_u(self.U, self.tau_f, self.tau_x, x)

    def stp_x(self, x):
        return stp_dgl_x(self.U, self.tau_f, self.tau_x, x)

    def stp_ur(self, x):
        return self.stp_x(x) * self.stp_u(x)

    def plot(self, **param):
        try:
            param['label']
        except KeyError:
            param['label'] = self.__repr__()

        x = np.arange(1., 181., 1.)*1e-3
        pl.plot(x*1e3, [self(xv) for xv in x], **param)
        pl.xlim((0, 180))
        pl.legend(loc="best")

        pl.xlabel("Input rate [Hz]")
        pl.ylabel("Channel activation [1]")

    def make_fun(self):
        return lambda x: np.sum(self.taus) * x * self.stp_ur(x)


class MFConstraint(object):

    def __init__(self, name, free_get, free_set, error_fun=lambda x: 0., bound_low=-1e5, bound_up=1e5, bound_hard=False):
        self.name = name
        self.free_set = free_set
        self.free_get = free_get
        self.error_fun = error_fun
        self.bound_low = bound_low
        self.bound_up = bound_up
        self.bound_hard = bound_hard

    @property
    def free(self):
        return self.free_get()

    @free.setter
    def free(self, value):
        self.free_set(value)

    @property
    def error(self):
        return self.error_fun()

    def __repr__(self):
        return "MFConstraint <%s>" % (self.name)


class MFState(object):

    def __init__(self, constraints, dependent_functions=[], bounds_check=True, bounds_correct=False, bounds_escalate=True, bounds_interrupt=False):
        self.constraints = constraints
        self.dependent_functions = dependent_functions
        self.bounds_check = bounds_check
        self.bounds_correct = bounds_correct
        self.bounds_escalate = bounds_escalate
        self.bounds_interrupt = bounds_interrupt

    def __getitem__(self, key):
        """Dictionarylike access to state values."""
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
        assert (len(value) == len(self.constraints))
        for i, c in enumerate(self.constraints):
            c.free = value[i]
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
                    return err_states*1e5  # return more detailed error to solver, i.e. which direction was wrong.
                

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
            if s < c.bound_low:
                if not c.bound_hard:  # reflect at bound
                    bstate[i] = 2 * c.bound_low - s
                else:  # set hard bound
                    bstate[i] = c.bound_low
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
        return "MFState<" + ", ".join(["%s: %.1f" % (c.name, c.free) for c in self.constraints]) + ">"

    def __call__(self, value, fun=lambda x: x):
        self.state = value
        error = self.error_checked
        return fun(error)


def gradient_solver(mfstate, p_0, dt=.1, tmax=30.):

    """Simple gradient descent along the error."""

    t = 0.
    state = np.array(p_0)
    states = [state]

    while t < tmax:
        state -= dt * np.array(mfstate(state))
        t += dt
        states.append(list(state))

    import pylab as pl
    
    states = np.array(states).T
    nspl = states.shape[0]
    for i in range(nspl):
        pl.subplot(nspl, 1, i+1)
        pl.plot(states[i, :])
    pl.show()

    # minimal solution object to return
    class sol:
        x = np.array(state)
        fun = np.array(mfstate.error)

    return sol()


class MFSolver(object):

    def __init__(self, mfstate, maxiter=20, solver="hybr", print_status=False, crit=1e-5, tol=1e-12, fail_val=None):

        self.mfstate = mfstate
        self.solver = solver
        self.maxiter = maxiter
        self.print_status = print_status
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

        print("\n-------------------")
        print("[%s] initializing minimization: %s" % (self.__class__.__name__, self.solver))
        print("[")

        # allow interruption of minimization loop, did not work otherwise.
        signal_handler = lambda signal, frame: sys.exit(0)
        signal.signal(signal.SIGINT, signal_handler)

        # loop over tries to solve the system
        while abs_err > crit:

            self.it += 1

            # interrupt upon reaching the maxiter.
            if self.it > self.maxiter:
                print("]\n[%s] maximum iterations reached" % self.__class__.__name__)
                self.state = "EXCEEDED"
                return self.finalize(min_sol)

            # get stochastic initial state
            up_dist = [c.bound_up - state_0[i] for i, c in enumerate(self.mfstate.constraints)]
            down_dist = [state_0[i] - c.bound_low for i, c in enumerate(self.mfstate.constraints)]
            max_dist = [min(up_dist[i], down_dist[i]) for i in range(len(self.mfstate.constraints))]
            p_0 = state_0 + (2.*np.random.rand(len(state_0))-1.) * np.array(max_dist) * noise_percent

            # solve
            if self.solver == "gradient":  # own implementation
                sol = gradient_solver(self.mfstate, p_0)
            else:  # scipy solvers
                sol = root(self.mfstate, p_0, jac=None, method=self.solver, tol=tol)

            # calculate the abs err, store minimal solution
            abs_err = max(abs(sol.fun))
            if abs_err < min_abs_err:
                min_abs_err = abs_err
                min_sol = sol
                sys.stdout.write('X')
            else:
                # display update
                sys.stdout.write('.')

            if self.it > 1 and self.it % 50 == 0 and self.it < self.maxiter:
                sys.stdout.write('\n ')
            sys.stdout.flush()

        print("]\n[%s] finished successfully" % self.__class__.__name__)
        self.state = "SUCCESS"
        return self.finalize(sol)

    def finalize(self, sol):
        self.mfstate.state = sol.x

        # set the state value to the fail_val if we did not converge
        if (self.fail_val is not None) and self.state != 'SUCCESS':
            self.mfstate.state = sol.x * self.fail_val

        self.error = sol.fun
        print(self.mfstate)
        print("-------------------\n")
        return self.mfstate


class MFSolver_RatesVoltages(MFSolver):

    def __init__(self, system, force_nmda=False, *args, **kwargs):

        # create constraints on the firing rates and mean voltages

        constraints = []
        functions = []

        for p in system.pops:
            constraints.append(
                MFConstraint(
                    "%s-%s" % (p.name, "rate_hz"),
                    partial(lambda x: x.rate_hz, p),
                    partial(lambda x, val: setattr(x, "rate_hz", val), p),
                    partial(lambda x: x.rate_hz-x.rate_pred, p),
                    0., 750.
                )
            )

            if p.has_nmda or force_nmda:
                print("Population %s has NMDA -> solving for voltages" % p.name)
                constraints.append(
                    MFConstraint(
                        "%s-%s" % (p.name, "v_mean"),
                        partial(lambda x: x.v_mean, p),
                        partial(lambda x, val: setattr(x, "v_mean", val), p),
                        partial(lambda x: x.v_mean-x.v_mean_prediction, p),
                        -80., 50.
                    )
                )
            else:
                functions.append(
                    partial(lambda x: setattr(x, "v_mean", x.v_mean_prediction), p),
                )

        state = MFState(constraints, dependent_functions=functions)
        super(MFSolver_RatesVoltages, self).__init__(state, *args, **kwargs)


def setup_brunel99(w_plus_val=2.5):

    # brunel 1999 system, one up state pop
    ff = .1
    w_plus = w_plus_val
    w_min = 1. - ff*(w_plus - 1.)/(1. - ff)
    initials = {'nu_plus': .025, 'nu_min': .0015, 'nu_i': .009}

    system = MFSystem("Brunel")
    pop_e1 = MFpop("Eup", params_standard["E"])
    pop_e1.n = 800
    pop_e2 = MFpop("Edown", params_standard["E"])
    pop_e2.n = 800
    pop_i = MFpop("I", params_standard["I"])
    pop_i.n = 200
    pop_e1.rate_ms = initials["nu_plus"]
    pop_e1.v_mean = -51.
    pop_e2.rate_ms = initials["nu_min"]
    pop_e2.v_mean = -55.
    pop_i.rate_ms = initials["nu_i"]
    pop_i.v_mean = -55.

    system.pops += [pop_e1, pop_e2, pop_i]

    # noise pops
    source_e_noise1 = MFSource("E_noise1", pop_e1)
    source_e_noise1.g_base = params_standard["E"]["gAMPA"]
    source_e_noise1.g_dyn = lambda: params_standard["E"]["nu_ext"] * params_standard["E"]["Cext"] * params_standard["E"]["tau_AMPA"]
    source_e_noise1.noise_tau = params_standard["E"]["tau_AMPA"]
    pop_e1.noise = source_e_noise1

    source_e_noise2 = MFSource("E_noise2", pop_e2)
    source_e_noise2.g_base = params_standard["E"]["gAMPA"]
    source_e_noise2.g_dyn = lambda: params_standard["E"]["nu_ext"] * params_standard["E"]["Cext"] * params_standard["E"]["tau_AMPA"]
    source_e_noise2.noise_tau = params_standard["E"]["tau_AMPA"]
    pop_e2.noise = source_e_noise2

    source_i_noise = MFSource("I_noise", pop_i)
    source_i_noise.g_base = params_standard["I"]["gAMPA"]
    source_i_noise.g_dyn = lambda: params_standard["I"]["nu_ext"] * params_standard["I"]["Cext"] * params_standard["I"]["tau_AMPA"]
    source_i_noise.noise_tau = params_standard["I"]["tau_AMPA"]
    pop_i.noise = source_i_noise

    # E->E NMDA
    syn_spec_nmda = {
        "tau_syn_rise": 1.,
        "tau_syn_d1": 100.,
        "tau_syn_d2": 100.,
        "balance": .5,
        "tau_x": 150.,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda1 = MFSource('EE Nmda1', pop_e1)
    source_ee_nmda1.g_base = params_standard["E"]["gNMDA"]
    source_ee_nmda1.g_dyn = lambda: (
            ff * w_plus * syn_ee_nmda(pop_e1.rate_ms) + (1. - ff) * w_min * syn_ee_nmda(pop_e2.rate_ms)
        ) * pop_e1.n
    source_ee_nmda1.is_nmda = True

    source_ee_nmda2 = MFSource('EE Nmda2', pop_e2)
    source_ee_nmda2.g_base = params_standard["E"]["gNMDA"]
    source_ee_nmda2.g_dyn = lambda: (
            ff * w_min * syn_ee_nmda(pop_e1.rate_ms) + (ff * w_plus + (1. - 2.*ff) * w_min) * syn_ee_nmda(pop_e2.rate_ms)
        ) * pop_e2.n
    source_ee_nmda2.is_nmda = True

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFSource('IE Nmda', pop_i)
    source_ie_nmda.g_base = params_standard["I"]["gNMDA"]
    source_ie_nmda.g_dyn = lambda: (
            ff * syn_ie_nmda(pop_e1.rate_ms) + (1. - ff) * syn_ie_nmda(pop_e2.rate_ms)
        ) * pop_e1.n
    source_ie_nmda.is_nmda = True

    # I->I GABA
    syn_spec_gaba = {
        "tau_syn_rise": 0.,
        "tau_syn_d1": 10.,
        "tau_syn_d2": 10.,
        "balance": .5,
    }
    syn_ii_gaba = Synapse(**syn_spec_gaba)
    source_ii_gaba = MFSource('II Gaba', pop_i)
    source_ii_gaba.g_base = params_standard["I"]["gGABA"]
    source_ii_gaba.g_dyn = lambda: syn_ii_gaba(pop_i.rate_ms) * pop_i.n
    source_ii_gaba.E_rev = -70.

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba1 = MFSource('EI Gaba', pop_e1)
    source_ei_gaba1.g_base = params_standard["E"]["gGABA"]
    source_ei_gaba1.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    source_ei_gaba1.E_rev = -70.
    source_ei_gaba2 = MFSource('EI Gaba', pop_e2)
    source_ei_gaba2.g_base = params_standard["E"]["gGABA"]
    source_ei_gaba2.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    source_ei_gaba2.E_rev = -70.

    return system


def setup_EI(has_nmda=True):

    # brunel 1999 system, one up state pop
    initials = {'nu_e': .003, 'nu_i': .009}

    mult = 0.25
    if has_nmda:
        mult = 1.

    system = MFSystem("EI")
    pop_e = MFpop("E", params_standard["E"])
    pop_e.n = 800
    pop_i = MFpop("I", params_standard["I"])
    pop_i.n = 200
    pop_e.rate_ms = initials["nu_e"]
    pop_e.v_mean = -51.
    pop_i.rate_ms = initials["nu_i"]
    pop_i.v_mean = -55.

    system.pops += [pop_e, pop_i]

    # noise pops
    source_e_noise = MFSource("E_noise", pop_e)
    source_e_noise.g_base = params_standard["E"]["gAMPA"]
    source_e_noise.g_dyn = lambda: params_standard["E"]["nu_ext"] * params_standard["E"]["Cext"] * params_standard["E"]["tau_AMPA"]
    source_e_noise.noise_tau = params_standard["E"]["tau_AMPA"]
    pop_e.noise = source_e_noise

    source_i_noise = MFSource("I_noise", pop_i)
    source_i_noise.g_base = params_standard["I"]["gAMPA"]
    source_i_noise.g_dyn = lambda: params_standard["I"]["nu_ext"] * params_standard["I"]["Cext"] * params_standard["I"]["tau_AMPA"]
    source_i_noise.noise_tau = params_standard["I"]["tau_AMPA"]
    pop_i.noise = source_i_noise

    # E->E NMDA
    syn_spec_nmda = {
        "tau_syn_rise": 1.,
        "tau_syn_d1": 100.,
        "tau_syn_d2": 100.,
        "balance": .5,
        "tau_x": 150.,  # depressing
    }
    syn_ee_nmda = Synapse(**syn_spec_nmda)

    source_ee_nmda = MFSource('EE Nmda', pop_e)
    source_ee_nmda.g_base = params_standard["E"]["gNMDA"] * mult
    source_ee_nmda.g_dyn = lambda: syn_ee_nmda(pop_e.rate_ms) * pop_e.n
    source_ee_nmda.is_nmda = has_nmda

    # E->I NMDA
    syn_ie_nmda = Synapse(**syn_spec_nmda)
    source_ie_nmda = MFSource('IE Nmda', pop_i)
    source_ie_nmda.g_base = params_standard["I"]["gNMDA"] * mult
    source_ie_nmda.g_dyn = lambda: syn_ie_nmda(pop_e.rate_ms) * pop_e.n
    source_ie_nmda.is_nmda = has_nmda

    # I->I GABA
    syn_spec_gaba = {
        "tau_syn_rise": 0.,
        "tau_syn_d1": 10.,
        "tau_syn_d2": 10.,
        "balance": .5,
    }
    syn_ii_gaba = Synapse(**syn_spec_gaba)
    source_ii_gaba = MFSource('II Gaba', pop_i)
    source_ii_gaba.g_base = params_standard["I"]["gGABA"]
    source_ii_gaba.g_dyn = lambda: syn_ii_gaba(pop_i.rate_ms) * pop_i.n
    source_ii_gaba.E_rev = -70.

    # I->E GABA
    syn_ei_gaba = Synapse(**syn_spec_gaba)
    source_ei_gaba = MFSource('EI Gaba', pop_e)
    source_ei_gaba.g_base = params_standard["E"]["gGABA"]
    source_ei_gaba.g_dyn = lambda: syn_ei_gaba(pop_i.rate_ms) * pop_i.n
    source_ei_gaba.E_rev = -70.

    return system
