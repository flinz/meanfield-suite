
class MFParams(dict):

    access_log = '_access_log'

    def __init__(self, **kwargs):
        super(MFParams, self).__init__(**{MFParams.access_log: {MFParams.access_log}}, **kwargs)

    def __getattr__(self, key):
        value = self[key]
        if key == MFParams.access_log:
            return value
        if value is not None:
            self[MFParams.access_log].add(key)
        return value

    def __setattr__(self, key, value):
        self[key] = value

    def all_keys_consumed(self):
        return set(self.keys()) == self[MFParams.access_log]

    def copy(self):
        cpy = super(MFParams, self).copy()
        del cpy[MFParams.access_log]
        return cpy

