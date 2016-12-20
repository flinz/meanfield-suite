import copy

from brian2 import in_unit, have_same_dimensions, DimensionMismatchError, get_dimensions


class MFParams2(object):

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

class MissingParameterError(Exception):
    pass

class MFParams(object):

    def __init__(self, params):
        self.underlying = dict(params)
        self.accesses = set()

    def __getitem__(self, key):
        value = self.underlying[key]
        self.accesses.add(key)
        return value

    def all_keys_consumed(self):
        not_used = self.accesses.difference(self.underlying.keys())
        return len(not_used) == 0

    def verify(self, expects):
        for param, unit in expects.items():

            if param not in self.underlying:
                raise MissingParameterError('parameter {} ({}) missing'.format(param, unit))

            if not have_same_dimensions(self.underlying[param], unit):
                raise DimensionMismatchError(
                    'parameter {} ({}) does not match dimension {}'.format(
                        param,
                        get_dimensions(self.underlying[param]),
                        unit
                    )
                )

    def __repr__(self):
        return self.underlying.__repr__()

