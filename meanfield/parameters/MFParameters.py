from collections import Mapping

from brian2 import have_same_dimensions, DimensionMismatchError, get_dimensions


class MissingParameterError(Exception):
    pass


class MFParameters(Mapping):

    def __init__(self, params):

        if isinstance(params, MFParameters):
            params = params.underlying

        self.underlying = dict(params)
        self.accesses = set()

    def __getitem__(self, key):
        self.accesses.add(key)
        return self.underlying[key]

    def __setitem__(self, key, value):
        self.underlying[key] = value

    def __iter__(self):
        return self.underlying.__iter__()

    def __len__(self):
        return self.underlying.__len__()

    def __add__(self, other):
        other_set = other if isinstance(other, dict) else other.underlying
        params = dict(**self.underlying)
        params.update(**other_set)
        return MFParameters(params)

    def all_keys_consumed(self):
        return self.accesses == set(self.underlying.keys())

    def verify(self, expectations):
        for param, unit in expectations.items():

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

        return self

    def fill(self, defaults):
        for param, value in defaults.items():

            if param not in self.underlying:
                self.underlying[param] = value

        return self

    def __repr__(self):
        items = ', '.join('{}: {}'.format(k, v) for k, v in self.underlying.items())
        return '{} {{ {} }}'.format(self.__class__.__name__, items)

