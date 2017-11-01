

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

