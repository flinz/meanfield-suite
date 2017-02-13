
def name2identifier(name):
    return name.replace(' ', '_')

def lazy(f):
    attr = '_lazy_' + f.__name__

    @property
    def wrapper(self):
        if not hasattr(self, attr):
            setattr(self, attr, f(self))
        return getattr(self, attr)
    return wrapper



