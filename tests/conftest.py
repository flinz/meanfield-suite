from time import sleep

import pytest
from brian2 import device


@pytest.yield_fixture(autouse=True)
def run_around_tests():
    yield
    device.reinit()
    device.activate()
    # workaround for cpp generation errors when running multiple tests in a row (reinit might take some time?)
    sleep(0.2)
