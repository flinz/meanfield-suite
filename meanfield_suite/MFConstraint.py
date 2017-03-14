
class MFConstraint(object):

    def __init__(self, name, free_get, free_set, error_fun=lambda x: 0., bound_down=-1e5, bound_up=1e5, bound_hard=False):
        self.name = name
        self.free_set = free_set
        self.free_get = free_get
        self.error_fun = error_fun
        self.bound_down = bound_down
        self.bound_up = bound_up
        self.bound_hard = bound_hard

    @property
    def free(self):
        return self.free_get()

    @free.setter
    def free(self, value):
        self.free_set(value)

    @property
    def error(self):
        print(self, self.error_fun())
        return self.error_fun()

    def __repr__(self):
        return "MFConstraint <%s>" % self.name
