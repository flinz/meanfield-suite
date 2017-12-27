import pytest
from brian2 import device


@pytest.yield_fixture(autouse=True)
def run_around_tests():
    yield
    print('REEESET')
    device.reinit()
    device.activate()
