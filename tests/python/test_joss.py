import sys, os
import numpy as np
from ctypes import cast, py_object, CDLL, RTLD_GLOBAL
from ase.build import molecule
from asi4py import ASI_ASE_calculator
from ase import units

def write_input(asi):
  from ase.calculators.dftb import Dftb
  calc = Dftb(label='Some_cluster',
        Hamiltonian_SCC='Yes',
        Hamiltonian_MaxAngularMomentum_='',
        Hamiltonian_MaxAngularMomentum_O='"p"',
        Hamiltonian_MaxAngularMomentum_H='"s"')
  calc.write_input(asi.atoms, properties=['forces'])

atoms = molecule('H2O')

atoms.calc = ASI_ASE_calculator(os.environ['ASI_LIB_PATH'], write_input, None, atoms)
atoms.calc.asi.keep_density_matrix = True
atoms.calc.asi.keep_hamiltonian = True
atoms.calc.asi.keep_overlap = True

print(f'E = {atoms.get_potential_energy():.6f}')

S = atoms.calc.asi.overlap_storage[(1,1)]
H = atoms.calc.asi.hamiltonian_storage[(1,1)]
DM = atoms.calc.asi.dm_storage.get((1,1), None)
DM_cnt = atoms.calc.asi.dm_calc_cnt[(1,1)]

print(f'Number of electrons = {np.sum(S*DM):.6f}')
print(f'Sum of eigenvalues = {np.sum(H*DM):.6f}')
print(f'Number of iterations = {DM_cnt}')
