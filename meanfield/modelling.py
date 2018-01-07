from typing import List, Union

import brian2

import MFSystem
from meanfield.inputs import MFInput
from meanfield.parameters import MFParameters, Connection
from meanfield.populations import MFPopulation
import numpy as np


def seed(seed: Union[int, None]) -> None:
    np.random.seed(seed)
    brian2.device.seed(seed)


def fully_connected(ctor: type, starts: Union[MFPopulation, List[MFPopulation]],
                    ends: Union[MFPopulation, List[MFPopulation]], parameters: Union[dict, MFParameters], *args,
                    name: str=None, self_loop_parameters: Union[None, Union[dict, MFParameters]] = None,
                    **kwargs) -> List[MFInput]:

    inputs = []

    if name is None:
        name = ''

    if not isinstance(starts, list):
        starts = list(starts)

    if not isinstance(ends, list):
        ends = list(ends)

    for start in starts:
        for end in ends:

            inp_name = f'{name}-{start.name}-{end.name}'
            inp = Connection.all_to_others() if start == end else Connection.all_to_all()
            inp_parameters = self_loop_parameters if start == end and self_loop_parameters is not None else parameters

            inp = ctor(start, end, inp_parameters, *args, **kwargs, name=inp_name, connection=inp)
            inputs.append(inp)

    return inputs


def multiple_populations(system: MFSystem, n_pop: int, ctor: type, n: Union[int, List[int]],
                         parameters: Union[dict, MFParameters], *args, name: str=None, **kwargs) -> List[MFPopulation]:

    populations = []

    if not isinstance(n, list):
        n = list(n) * n_pop

    if len(n) != n_pop:
        raise ValueError(f'populations sizes (array) should have length {n_pop}')

    for i, size in enumerate(n):

        pop_name = None if name is None else f'{name}-{i}'
        pop = ctor(size, parameters, *args, **kwargs, name=pop_name)
        populations.append(pop)
        system.add_population(pop)

    return populations


