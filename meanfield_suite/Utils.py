
from itertools import count
iid = count()


def lazy(f):
    attr = '_lazy_' + f.__name__

    @property
    def wrapper(self):
        if not hasattr(self, attr):
            setattr(self, attr, f(self))
        return getattr(self, attr)
    return wrapper



