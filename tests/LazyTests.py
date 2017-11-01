import unittest

from utils import lazy


class Counter(object):

    def __init__(self):
        self.n = 0

    @lazy
    def m1(self):
        self.n += 1
        return self.n

    @lazy
    def m2(self):
        self.n += 1
        return self.n


class LazyTests(unittest.TestCase):

    def testLazyness(self):
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

    def testSharedLazyStructure(self):
        c1 = Counter()
        c2 = Counter()
        assert(c1.n == 0)
        assert(c2.n == 0)
        assert(c1.m1 == 1)
        assert(c1.n == 1)
        assert(c2.n == 0)
