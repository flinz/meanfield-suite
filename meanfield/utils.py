from brian2 import NeuronGroup, Synapses, PoissonInput, magic_network, Equations


def name2identifier(name):
    return name.replace(' ', '_')


def lazyproperty(fun):
    attr = '_lazy_' + fun.__name__

    @property
    def wrapper(self):
        if not hasattr(self, attr):
            setattr(self, attr, fun(self))
        return getattr(self, attr)
    return wrapper


def brian2_introspect(net, globals):

    def introspect_population(pop):
        if False:
            print(pop.thresholder)
            print(pop.resetter)
            print(pop.events['spike'])
            print(pop.event_codes.get('spike'))

        print('{} [{}] n: {}'.format(
            pop.__class__.__name__,
            pop.name,
            pop.N
        ))
        print(Equations(str(pop.user_equations), **globals))

    def introspect_source(source):
        print("{} [{}]".format(
            source.__class__.__name__,
            source.name,
        ))

    sources = {}

    for source in net.objects:
        if isinstance(source, Synapses):
            if source.target.id not in sources:
                sources[source.target.id] = []
            sources[source.target.id].append(source)

        if isinstance(source, PoissonInput):
            if source._group.id not in sources:
                sources[source._group.id] = []
            sources[source._group.id].append(source)

    for population in net.objects:
        if isinstance(population, NeuronGroup):
            introspect_population(population)
            print()

            for source in sources.get(population.id, []):
                introspect_source(source)

            print()




