import brian2 as b2
import brian2.units
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import create_name, create_identifier, reset_lazyproperty


class MFSystem(object):

    def __init__(self, *populations: MFPopulation, name: str=None):
        self.name = name if name else create_name(self)
        self.ref = create_identifier(self.name)
        self.populations = set(populations)

    def add_population(self, *populations: MFPopulation) -> 'MFSystem':
        self.populations.update(populations)
        return self

    def __getitem__(self, population_name: str) -> MFPopulation:
        """Dictionary-like access to state values."""
        for p in self.populations:
            if p.name == population_name:
                return p
        raise KeyError

    def setup_brian2(self, dt=0.01 * b2.units.ms, codegen: str='cpp_standalone', build_on_run: bool=True):
        b2.set_device(codegen, build_on_run=build_on_run)
        b2.defaultclock.dt = dt

    def reset_brian2(self, **kwargs):
        b2.device.reinit()
        b2.device.activate()
        self.setup_brian2(**kwargs)

        for p in self.populations:
            reset_lazyproperty(p, 'brian2')

            for i in p.inputs:
                reset_lazyproperty(i, 'brian2')

            for n in p.noises:
                reset_lazyproperty(n, 'brian2')

    def collect_brian2_network(self):
        for net in b2.Network.__instances__():
            if self.ref == net().name:
                # FIXME warning about call to reset brian2 and probably slow as codegen recreated
                break

        net = b2.Network(name=self.ref)

        for p in self.populations:
            net.add(p.brian2)

            for i in p.inputs:
                net.add(i.brian2)

            for n in p.noises:
                net.add(n.brian2)

        return net

    def introspect(self, indent=0) -> str:
        builder = []
        spaces = ' ' * indent
        builder.append(str(self))

        for pop in self.populations:
            builder.append('\n{}  - {}'.format(spaces, pop.introspect(indent=indent + 2)))

        builder.append('')
        return '\n'.join(builder)

    def print_introspect(self) -> None:
        print(self.introspect(), flush=True)

    def __repr__(self):
        return "{} [{}] ({} pops)".format(self.__class__.__name__, self.name, len(self.populations))
