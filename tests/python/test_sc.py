from mpi4py import MPI

from ctypes import c_double, c_void_p, CFUNCTYPE, cast, py_object, c_int, POINTER
from pyasi.pyasi import DFT_C_API
import numpy as np
from numpy.testing import assert_allclose
from ase.build import molecule
from ase import units
import os, sys
from time import sleep
from scipy.spatial import KDTree
from scipy.interpolate import Rbf as interpolator
from mpiprint import parprint, ordprint

@CFUNCTYPE(None, c_void_p, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double))
def ext_pot_func(ptr, n, coords_ptr, potential_ptr, potential_grad_ptr):
  func, _, _ = cast(ptr, py_object).value
  coords = np.ctypeslib.as_array(coords_ptr, shape=(n, 3))
  if potential_ptr:
    potential = np.ctypeslib.as_array(potential_ptr, shape=(n, ))
    potential[:] = func(coords)
  if potential_grad_ptr:
    potential_grad = np.ctypeslib.as_array(potential_grad_ptr, shape=(n, 3))
    potential_grad[:] = 0
    

X, Y, Z  = np.mgrid[-3:3:5j, -3:3:5j, -3:3:5j]
shifts = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

def get_rdf_potential(asi, target_pos=None):
  coords = (target_pos[:, None, :] + shifts[None, :, :]).reshape((-1, 3))
  #coords = np.vstack([target_pos, target_pos + [1.5, 0, 0], target_pos + [0., 1.5, 0], target_pos + [0, 0, 1.5], target_pos - [1.5, 0, 0], target_pos - [0., 1.5, 0], target_pos - [0, 0, 1.5]])
  esp, _ = asi.calc_esp(coords) if asi is not None else (np.zeros((coords.shape[0], )), None)
  func = interpolator(coords, esp, neighbors=128, smoothing=0.0)
  return ext_pot_func, (func, esp, coords)

def mix_rdf_potential(pot1, pot2):
  _, (_, esp1, coords1) = pot1
  _, (_, esp2, coords2) = pot2
  assert_allclose(coords1, coords2)
  esp = esp1 * 0.5 + esp2 * 0.5
  func = interpolator(coords1, esp, neighbors=128, smoothing=0.0)
  return ext_pot_func, (func, esp, coords1)



lib_file, get_pot, mixer = os.environ['ASI_LIB_PATH'], get_rdf_potential, mix_rdf_potential

if "aims" in lib_file:
  def initializer(asi):
    from ase.calculators.aims import Aims
    calc = Aims(xc='pw-lda')
    calc.write_input(asi.atoms)
else:
  def initializer(asi):
    from ase.calculators.dftb import Dftb
    calc = Dftb(label='Some_cluster',
          Hamiltonian_SCC='Yes',
          Hamiltonian_MaxAngularMomentum_='',
          Hamiltonian_MaxAngularMomentum_O='"p"',
          Hamiltonian_MaxAngularMomentum_H='"s"')
    calc.write_input(asi.atoms, properties=['forces'])

d = 2.75

h2o1 = molecule('H2O')
h2o2 = h2o1.copy()
h2o2.translate([0, 0, d])

with DFT_C_API(lib_file, initializer, MPI.COMM_WORLD, h2o1 + h2o2, 'asi0') as asi0:
  asi0.run()
  E0 = asi0.total_energy
parprint (f'E0={E0:.6f}')

with DFT_C_API(lib_file, initializer, MPI.COMM_WORLD, h2o1, 'asi1') as asi1:
  asi1.run()
  E01 = asi1.total_energy
  pot1 = get_pot(asi1, h2o2.positions / units.Bohr)

parprint (f'E01={E01:.6f}')

with DFT_C_API(lib_file, initializer, MPI.COMM_WORLD, h2o2, 'asi2') as asi2:
  asi2.run()
  E02 = asi2.total_energy
  pot2 = get_pot(asi2, h2o1.positions / units.Bohr)
parprint (f'E02={E02:.6f}')
dE = E0 - (E01 + E02)
parprint (f'dE={dE:.6f}')

#pot1 = get_pot(None, h2o2.positions / units.Bohr)
#pot2 = get_pot(None, h2o1.positions / units.Bohr)



for i in range (3):

  with DFT_C_API(lib_file, initializer, MPI.COMM_WORLD, h2o2, 'asi2') as asi2:
    asi2.register_external_potential(*pot1)
    asi2.run()
    new_pot2 = get_pot(asi2, h2o1.positions / units.Bohr)
    E2 = asi2.total_energy

  with DFT_C_API(lib_file, initializer, MPI.COMM_WORLD, h2o1, 'asi1') as asi1:
    asi1.register_external_potential(*pot2)
    asi1.run()
    new_pot1 = get_pot(asi1, h2o2.positions / units.Bohr)
    E1 = asi1.total_energy

  pot1 = mixer(pot1, new_pot1)
  pot2 = mixer(pot2, new_pot2)
  parprint (f"{d:.2f} {E0 - E01 - E02:.6f} {E1 - E01:.6f} {E2 - E02:.6f}")


