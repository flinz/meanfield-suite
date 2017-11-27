class MFSystem(object):

    def __init__(self, name):
        self.name = name
        self.pops = []

    def __getitem__(self, key):
        """Dictionary-like access to state values."""
        names = [x.name for x in self.pops]
        try:
            idx = names.index(key)
        except ValueError:
            raise KeyError
        return self.pops[idx]

    def introspect(self, indent=0) -> str:
        builder = []
        spaces = ' ' * indent
        builder.append(str(self))

        for pop in self.pops:
            builder.append('\n{}  - {}'.format(spaces, pop.introspect(indent=indent + 2)))

        builder.append('')
        return '\n'.join(builder)

    def print_introspect(self) -> None:
        print(self.introspect(), flush=True)

    def __repr__(self):
        return "{} [{}] ({} pops)".format(self.__class__.__name__, self.name, len(self.pops))
