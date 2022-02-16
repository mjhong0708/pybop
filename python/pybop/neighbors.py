import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList, NewPrimitiveNeighborList


def get_neighbor_positions(atoms: Atoms, r_c: float) -> np.ndarray:
    """Get positions of neighbors in a single `Atoms` object.
    For homogeneity of array, neighborlist is padded with `1e10`.

    Args:
        atoms ([type]): [description]
        r_c ([type]): [description]

    Returns:
        [type]: [description]
    """
    calc = NeighborList(
        cutoffs=[r_c / 2.0] * len(atoms),
        self_interaction=False,
        bothways=True,
        skin=0.0,
        primitive=NewPrimitiveNeighborList,
    )
    calc.update(atoms)
    neighborlist = [calc.get_neighbors(index) for index in range(len(atoms))]
    max_neighbor_num = max(map(lambda x: len(x[0]), neighborlist))

    neighbor_positions = 1e10 * r_c * np.ones((len(atoms), max_neighbor_num, 3), dtype=np.float64)
    for i, atom in enumerate(atoms):
        indices, offset = calc.get_neighbors(atom.index)
        pos = atoms.positions[indices] + offset @ atoms.get_cell()
        neighbor_positions[i, : len(pos)] = pos

    return neighbor_positions
