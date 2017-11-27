from typing import Callable

from brian2 import Synapses


class ConnectionStrategy:

    def __init__(self, name: str, theory: Callable[[int], int], simulation: Callable[[Synapses], Synapses]):
        self.name = name
        self.theory = theory
        self.simulation = simulation

    def __repr__(self):
        return self.name


def random(p):
    def theory(n):
        return round(p * n)

    def simulation(s):
        s.connect(p=p)
        return s

    return ConnectionStrategy('random {}'.format(p), theory, simulation)


def all_to_all():
    def theory(n):
        return n

    def simulation(s):
        s.connect()
        return s

    return ConnectionStrategy('all-to-all', theory, simulation)


def all_to_others():
    def theory(n):
        return n - 1

    def simulation(s):
        s.connect(condition='j != i')
        return s

    return ConnectionStrategy('all-to-others', theory, simulation)


def one_to_one():
    def theory(n):
        return 1

    def simulation(s):
        s.connect(j='i')
        return s

    return ConnectionStrategy('one-to-one', theory, simulation)


