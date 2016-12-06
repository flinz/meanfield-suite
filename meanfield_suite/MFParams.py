import copy

class MFParams(object):

    def __init__(self, **kwargs):
        self.__dict__['underlying'] = dict(**kwargs)
        self.__dict__['access_log'] = set()

    def __getattr__(self, key):
        value = self.__dict__['underlying'][key]
        self.__dict__['access_log'].add(key)
        return value

    def __setattr__(self, key, value):
        self.__dict__['underlying'][key] = value

    def all_keys_consumed(self):
        return set(self.__dict__['underlying'].keys()) == self.__dict__['access_log']

    def to_dict(self):
        return copy.deepcopy(self.__dict__['underlying'])
