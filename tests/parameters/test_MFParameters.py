from meanfield.parameters.MFParameters import MFParameters

class TestMFParams(object):

    def test_constructor(self):
        assert MFParameters({}).underlying == {}
        assert MFParameters({'a': 1}).underlying == {'a': 1}

    def test_getter(self):
        try:
            _ = MFParameters({})['a']
            assert False
        except:
            pass
        assert MFParameters({'a': 1})['a'] == 1

#    def testSetter(self):
#        p = MFParams({'a': 1})
#        p.a = 2
#        p.b = 3
#        assert p.a == 2
#        assert p.b == 3

    def test_all_keys_consumed(self):
        assert MFParameters({}).all_keys_consumed()
        p = MFParameters({'a': 1})
        assert not p.all_keys_consumed()
        _ = p['a']
        assert p.all_keys_consumed()
