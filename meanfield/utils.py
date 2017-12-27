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

    def introspect_input(input):
        print("{} [{}]".format(
            input.__class__.__name__,
            input.name,
        ))

    inputs = {}

    for source in net.objects:
        if isinstance(source, Synapses):
            if source.target.id not in inputs:
                inputs[source.target.id] = []
            inputs[source.target.id].append(source)

        if isinstance(source, PoissonInput):
            if source._group.id not in inputs:
                inputs[source._group.id] = []
            inputs[source._group.id].append(source)

    for population in net.objects:
        if isinstance(population, NeuronGroup):
            introspect_population(population)
            print()

            for source in inputs.get(population.id, []):
                introspect_input(source)

            print()




