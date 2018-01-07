import matplotlib.pyplot as plt
from brian2 import SpikeMonitor, units, PopulationRateMonitor
from typing import List


def activities(spike_monitors: List[SpikeMonitor], n: int = 15) -> None:

    for i, sm in enumerate(spike_monitors):
        plt.plot(sm.t / units.ms, sm.i + (len(spike_monitors) - i - 1) * n, '.', markersize=2, label=sm.name)

    plt.title('Population activities ({} neurons/pop)'.format(n))
    plt.xlabel('ms')
    plt.yticks([])
    plt.legend()


def rates(rate_monitors: List[PopulationRateMonitor], smoothing: units.second = 25 * units.ms) -> None:

    for rm in rate_monitors:
        plt.plot(rm.t / units.ms, rm.smooth_rate(width=smoothing) / units.Hz, label=rm.name)

    plt.title('Population rates')
    plt.xlabel('ms')
    plt.ylabel('Hz')
    plt.legend()
