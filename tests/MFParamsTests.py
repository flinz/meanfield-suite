import unittest

from MFParams import MFParams

class MFParamsTests(unittest.TestCase):

    def testConstructor(self):
        assert MFParams().copy() == {}
        assert MFParams(a=1).copy() == {'a': 1}

    def testGetter(self):
        try:
            _ = MFParams().a
            assert False
        except KeyError:
            pass
        assert MFParams(a=1).a == 1

    def testSetter(self):
        p = MFParams(a=1)
        p.a = 2
        p.b = 3
        assert p.a == 2
        assert p.b == 3

    def testAllKeysConsumed(self):
        assert MFParams().all_keys_consumed()
        assert not MFParams(a=1).all_keys_consumed()
        p = MFParams(a=1)
        p.b = 2
        assert not p.all_keys_consumed()
        _ = p.a
        assert not p.all_keys_consumed()
        _ = p.b
        assert p.all_keys_consumed()

    def testDestructorCheck(self):

        params = MFParams()
