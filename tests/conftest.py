from time import sleep

import pytest
from brian2 import device, gc


@pytest.yield_fixture(autouse=True)
def run_around_tests():
    yield
    device.reinit()
    device.activate()
    # workaround for cpp generation errors when running multiple tests in a row (reinit might take some time?)
    # Network.__instances__() shows multiple instance of the network containing same objects
    sleep(0.5)

