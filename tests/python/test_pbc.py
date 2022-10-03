import sys, os
import numpy as np
from mpi4py import MPI
from mpiprint import parprint, ordprint
from ctypes import cast, py_object, CDLL, RTLD_GLOBAL

from ase.build import molecule, bulk
from ase.io import read, write
from asi4py.asecalc import ASI_ASE_calculator
from ase import units

ASI_LIB_PATH = os.environ['ASI_LIB_PATH']
asilib = CDLL(ASI_LIB_PATH, mode=RTLD_GLOBAL)

if asilib.ASI_flavour() == 1:
  def init_via_ase(asi):
    from ase.calculators.aims import Aims
    calc = Aims(xc='pbe', 
      relativistic="atomic_zora scalar",
      occupation_type="gaussian 0.010",
      sc_accuracy_etot=1e-06,
      k_grid=(2,2,2),
      sc_accuracy_forces=1e-1, # enables force calculation
    )
    calc.write_input(asi.atoms)
else:
  def init_via_ase(asi):
    from ase.calculators.dftb import Dftb
    calc = Dftb(label='Si_crystall',
          kpts=(2,2,2),
          Hamiltonian_SCC='Yes',
          Hamiltonian_MaxAngularMomentum_='',
          Hamiltonian_MaxAngularMomentum_Si='"p"')
    calc.write_input(asi.atoms, properties=['forces'])


atoms = bulk('Si')
atoms.pbc=True

atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_via_ase, None, atoms)
print(f'atoms.calc created!')
E = atoms.get_potential_energy()
print(f'E = {E:.6f}')
F = atoms.get_forces()
print(F)
S = atoms.get_stress()
print(S)

