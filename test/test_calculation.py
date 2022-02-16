import numpy as np
from ase.build import bulk
from pybop.bop import calculate_bop

atoms = bulk("Pt", "fcc", a=3.92, cubic=True)
atoms = atoms.repeat((3, 3, 3))  # Supercell

bop_params = 0.5 * np.ones((len(atoms), 8), dtype=np.float64)

energy = calculate_bop(atoms, r_cut=6.5, params_list=bop_params, n_jobs=8)
print(energy)
