from ase.build import bulk
from pybop.bop import calculate_bop

Pt = bulk("Pt", "fcc", a=3.92, cubic=True)
Pt = Pt.repeat((3, 3, 3))

bop_params = [[0.5] * 8] * len(Pt)

energy = calculate_bop(Pt, r_cut=6.5, params_list=bop_params, n_jobs=8)
print(energy)
