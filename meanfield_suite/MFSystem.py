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
