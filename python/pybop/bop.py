from joblib import Parallel, delayed
from .pybop import calculate_bop as _bop_rs
from .neighbors import get_neighbor_positions


def calculate_bop(atoms, r_cut, params_list, n_jobs=1):
    center_positions = atoms.positions
    neighbors_positions = get_neighbor_positions(atoms, r_cut)
    energy = 0
    if n_jobs > 1:
        atomic_energies = Parallel(n_jobs=n_jobs)(
            delayed(_bop_rs)(center_positions[i], neighbors_positions[i], params_list[i], r_cut)
            for i in range(len(atoms))
        )
        # atomic_energies = (
        #     _bop_rs(center_positions[i], neighbors_positions[i], params_list[i], r_cut)
        #     for i in range(len(atoms))
        # )
        energy = sum(atomic_energies)
    else:
        for center, neighbors, params in zip(center_positions, neighbors_positions, params_list):
            energy += _bop_rs(center, neighbors, params, r_cut)
    return energy
