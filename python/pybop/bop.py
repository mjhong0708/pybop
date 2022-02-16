from typing import Sequence, Union

import numpy as np
from ase import Atoms
from joblib import Parallel, delayed

from .neighbors import get_neighbor_positions
from .pybop import calculate_bop as _bop_rs


def calculate_bop(
    atoms: Atoms,
    r_cut: float,
    params_list: Union[Sequence[Sequence[float]], np.ndarray],
    n_jobs: int = 1,
) -> float:
    """Calculates BOP for an Atoms object.

    Args:
        atoms (Atoms): Atoms object
        r_cut (float): cutoff radius in Angstrom
        params_list (Union[Sequence[Sequence[float]], np.ndarray]): list of local BOP parameters.
            Each element is a list of 8 parameters. (N, 8) numpy array is also accepted.
        n_jobs (int, optional): The number of parallel jobs. Defaults to 1.

    Returns:
        float: energy of Atoms
    """
    center_positions = atoms.positions
    neighbors_positions = get_neighbor_positions(atoms, r_cut)
    energy = 0
    if n_jobs > 1:
        atomic_energies = Parallel(n_jobs=n_jobs)(
            delayed(_bop_rs)(center_positions[i], neighbors_positions[i], params_list[i], r_cut)
            for i in range(len(atoms))
        )
        energy = sum(atomic_energies)
    else:
        for center, neighbors, params in zip(center_positions, neighbors_positions, params_list):
            energy += _bop_rs(center, neighbors, params, r_cut)
    return energy
