class MFSource(object):
    """Source: a synapse coupled to pops"""

    def __init__(self, name, pop):

        # TODO source vs synapse
        # TODO sometimes setup, sometimes constructed

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
