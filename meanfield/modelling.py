from typing import List, Union

from brian2 import device, SpikeMonitor, units, PopulationRateMonitor

import MFSystem
from meanfield.inputs.MFInput import MFInput
from meanfield.parameters import Connection
from meanfield.parameters.MFParameters import MFParameters
from meanfield.populations.MFPopulation import MFPopulation
import numpy as np

from utils import listify


def seed(seed: Union[int, None]) -> None:
    np.random.seed(seed)
    device.seed(seed)


def fully_connected(ctor, starts: Union[MFPopulation, List[MFPopulation]],
                    ends: Union[MFPopulation, List[MFPopulation]], parameters: Union[dict, MFParameters], *args,
                    name: str=None, self_loop_parameters: Union[None, Union[dict, MFParameters]] = None,
                    **kwargs) -> List[MFInput]:

    inputs = []

    if name is None:
        name = ''

    starts = listify(starts)
    ends = listify(ends)

    for start in starts:
        for end in ends:

            inp_name = f'{name} {start.name}-{end.name}'
            inp = Connection.all_to_others() if start == end else Connection.all_to_all()
            inp_parameters = self_loop_parameters if start == end and self_loop_parameters is not None else parameters

            inp = ctor(start, end, inp_parameters, *args, **kwargs, name=inp_name, connection=inp)
            inputs.append(inp)

    return inputs


def noise_connected(ctor, n: int, rate: units.Hz, ends: Union[MFPopulation, List[MFPopulation]],
                          parameters: Union[dict, MFParameters], *args, name: str=None, **kwargs) -> List[MFInput]:

    inputs = []

    if name is None:
        name = ''

    ends = listify(ends)

    for end in ends:

        inp_name = f'{name} {end.name}'

        inp = ctor(n, rate, end, parameters, *args, **kwargs, name=inp_name)
        inputs.append(inp)

    return inputs


def multiple_populations(n_pop: int, ctor, n: Union[int, List[int]], parameters: Union[dict, MFParameters],
                         *args, name: str=None, **kwargs) -> List[MFPopulation]:

    populations = []

    if not isinstance(n, list):
        n = [n] * n_pop

    if len(n) != n_pop:
        raise ValueError(f'populations sizes (array) should have length {n_pop}')

    for i, size in enumerate(n):

        pop_name = None if name is None else f'{name}-{i}'
        pop = ctor(size, parameters, *args, **kwargs, name=pop_name)
        populations.append(pop)

    return populations


def brian2_spike_monitors(populations: Union[MFPopulation, List[MFPopulation]], n: int = 15) -> List[SpikeMonitor]:

    populations = listify(populations)

    return [SpikeMonitor(p.brian2[:n], name=f'spike_{p.ref}') for p in populations]


def brian2_rate_monitors(populations: Union[MFPopulation, List[MFPopulation]]) -> List[PopulationRateMonitor]:

    populations = listify(populations)

    return [PopulationRateMonitor(p.brian2, name=f'rate_{p.ref}') for p in populations]

