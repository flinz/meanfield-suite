import brian2 as b2
import brian2.units
from meanfield.populations.MFPopulation import MFPopulation
from meanfield.utils import create_name, create_identifier, reset_lazyproperty
from graphviz import Digraph
import logging

logger = logging.getLogger('system')

class MFSystem(object):

    def __init__(self, *populations: MFPopulation, name: str=None):
        self.name = name if name else create_name(self)
        self.ref = create_identifier(self.name)
        self.populations = list(populations)

    def add_population(self, *populations: MFPopulation) -> 'MFSystem':
        self.populations.extend(populations)
        return self

    def __getitem__(self, population_name: str) -> MFPopulation:
        """Dictionary-like access to state values."""
        for p in self.populations:
            if p.name == population_name:
                return p
        raise KeyError

    def reset_lazyness(self) -> None:
        for p in self.populations:
            reset_lazyproperty(p, 'brian2')

            for i in p.inputs:
                reset_lazyproperty(i, 'brian2')

            for n in p.noises:
                reset_lazyproperty(n, 'brian2')

    def collect_brian2_network(self, *more_objects: b2.BrianObject):
        for net in b2.Network.__instances__():
            if self.ref == net().name:
                logger.warning('before recollecting the same network, you need to call reset_lazyness to clean up')
                break

        net = b2.Network(*more_objects, name=self.ref)

        for p in self.populations:
            net.add(p.brian2)

            for i in p.inputs:
                net.add(i.brian2)

            for n in p.noises:
                net.add(n.brian2)

        return net

    def brian2_run(self, *more_objects: b2.BrianObject, t: b2.units.second, **kwargs):
        self.collect_brian2_network(*more_objects).run(t, **kwargs)

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

    def graph(self, engine='dot'):
        g = Digraph(name=self.name, engine=engine)
        g.graph_attr.update(
            overlap='false',
            splines='true',
            rankdir='LR',
            fontname='Helvetica Neue',
        )
        g.node_attr.update(
            fontsize='11',
            fontname='Helvetica Neue',
        )
        g.edge_attr.update(
            fontsize='9',
            fontname='Helvetica Neue',
        )

        def table(**kwargs):
            rows_html = ''.join(tr(k, v) for k, v in kwargs.items())
            return f'''<
                <table
                    cellspacing="0" 
                    border="0"
                    cellpadding="3"
                    cellborder="1" 
                >
                    {rows_html}
                </table>
            >'''

        def tr(*args):
            cells_html = ''.join(td(c) for c in args)
            return f'<tr>{cells_html}</tr>'

        def td(value):
            return f'<td padding="50px">{value}</td>'

        for p in self.populations:
            g.node(p.ref, table(**{p.__class__.__name__: p.name}, n=p.n, **p.parameters), shape='plaintext')

            for i in p.inputs:
                g.edge(i.origin.ref, p.ref, label=i.name)

            for n in p.noises:
                g.edge(n.ref, p.ref, style='dotted', arrowhead='empty')

        with g.subgraph(name='cluster_noises') as c:
            c.node_attr.update(style='filled')
            c.attr(label='External noises')
            c.attr(style='filled')
            c.attr(color='#efefef')

            for p in self.populations:
                for n in p.noises:
                    c.node(n.ref, table(**{n.__class__.__name__: n.name}, n=n.n, rate=n.rate), shape='plaintext')

        return g

    def __repr__(self):
        return "{} [{}] ({} pops)".format(self.__class__.__name__, self.name, len(self.populations))
