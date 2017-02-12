import unittest

from Utils import lazy


class LazyTests(unittest.TestCase):

    @lazy
    def compute1(self):
        self.counter += 1

    @lazy
    def compute2(self):
        self.counter += 1

    def testLazyness(self):
        self.counter = 0
        self.compute1()
        assert(self.counter == 1)
        self.compute1()
        assert(self.counter == 1)
        self.compute1()
        assert(self.counter == 1)
        self.compute2()
        assert(self.counter == 2)
        self.compute2()
        assert(self.counter == 2)
        self.compute2()
        assert(self.counter == 2)
