from enum import Enum

class ConnectionStrategy:

    def __init__(self, theory, simulation):
        self.theory = theory
        self.simulation = simulation


def random(p):
    def theory(n):
        return round(p * n)

    def simulation(s):
        s.connect(p=p)
        return s

    return ConnectionStrategy(theory, simulation)


def all_to_all():
    def theory(n):
        return n
    # FIXME n-1 self connection, else n

    def simulation(s):
        s.connect(condition='j != i')
        return s

    return ConnectionStrategy(theory, simulation)


def one_to_one():
    def theory(n):
        return 1

    def simulation(s):
        s.connect(j='i')
        return s

    return ConnectionStrategy(theory, simulation)


