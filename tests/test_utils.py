from meanfield.utils import lazyproperty


class Counter(object):

    def __init__(self):
        self.n = 0

    @lazyproperty
    def m1(self):
        self.n += 1
        return self.n

    @lazyproperty
    def m2(self):
        self.n += 1
        return self.n


class TestUtils(object):

    def test_lazy_property(self):
        c = Counter()
        assert(c.n == 0)
        assert(c.m1 == 1)
        assert(c.n == 1)
        assert(c.m1 == 1)
        assert(c.n == 1)
        assert(c.m2 == 2)
        assert(c.n == 2)
        assert(c.m2 == 2)
        assert(c.n == 2)

    def test_shared_lazy_property(self):
        c1 = Counter()
        c2 = Counter()
        assert(c1.n == 0)
        assert(c2.n == 0)
        assert(c1.m1 == 1)
        assert(c1.n == 1)
        assert(c2.n == 0)
