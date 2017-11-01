from meanfield.parameters.MFParams import MFParams

class TestMFParams(object):

    def test_constructor(self):
        assert MFParams({}).underlying == {}
        assert MFParams({'a': 1}).underlying == {'a': 1}

    def test_getter(self):
        try:
            _ = MFParams({})['a']
            assert False
        except:
            pass
        assert MFParams({'a': 1})['a'] == 1

#    def testSetter(self):
#        p = MFParams({'a': 1})
#        p.a = 2
#        p.b = 3
#        assert p.a == 2
#        assert p.b == 3

    def test_all_keys_consumed(self):
        assert MFParams({}).all_keys_consumed()
        p = MFParams({'a': 1})
        assert not p.all_keys_consumed()
        _ = p['a']
        assert p.all_keys_consumed()
