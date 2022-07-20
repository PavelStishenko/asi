from mpi4py import MPI

from ctypes import c_double, c_void_p, CFUNCTYPE, cast, py_object, c_int, POINTER
from asi4py import DFT_C_API, init_dftbp, init_aims
import numpy as np
from numpy.testing import assert_allclose
from ase.build import molecule
from ase import units
import os, sys
from time import sleep
from scipy.spatial import KDTree
from scipy.interpolate import RBFInterpolator as interpolator

@CFUNCTYPE(None, c_void_p, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double))
def ext_pot_nearest(ptr, n, coords_ptr, potential_ptr, potential_grad_ptr):
    tree, esp, esp_grad = cast(ptr, py_object).value
    coords = np.ctypeslib.as_array(coords_ptr, shape=(n, 3))
    d,i = tree.query(coords, 1, eps=0.001)
    assert np.all(d < 1e-6), d 
    
    if potential_ptr:
      potential = np.ctypeslib.as_array(potential_ptr, shape=(n, ))
      potential[:] = esp[i]

    if potential_grad_ptr:
      potential_grad = np.ctypeslib.as_array(potential_grad_ptr, shape=(n, 3))
      potential_grad[:] = esp_grad[i]

def get_tree_potential(asi, target_pos=None):
  tree = KDTree(target_pos)
  esp, esp_grad = asi.calc_esp(target_pos)
  return ext_pot_nearest, (tree, esp.copy(), esp_grad.copy())

def mix_tree_potential(pot1, pot2):
  _, (tree, esp1, esp_grad1) = pot1
  _, (tree, esp2, esp_grad2) = pot2
  return ext_pot_nearest, (tree, (esp1 + esp2)*0.5, (esp_grad1 + esp_grad2)*0.5)
 
ext_pot_func_CNT = 0

@CFUNCTYPE(None, c_void_p, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double))
def ext_pot_func(ptr, n, coords_ptr, potential_ptr, potential_grad_ptr):
  global ext_pot_func_CNT
  print ('ext_pot_func', n, ext_pot_func_CNT)
  ext_pot_func_CNT += 1
  func, _, _ = cast(ptr, py_object).value
  coords = np.ctypeslib.as_array(coords_ptr, shape=(n, 3))
  potential = np.ctypeslib.as_array(potential_ptr, shape=(n, ))  
  potential[:] = func(coords)


initializer, lib_file, get_pot, mixer = init_dftbp, os.environ['ASI_LIB_PATH'], get_tree_potential, mix_tree_potential
d = 3.5

h2o1 = molecule('H2O')
h2o2 = h2o1.copy()
h2o2.translate([0, 0, d])

with DFT_C_API(lib_file, initializer, MPI.COMM_WORLD, h2o1 + h2o2, 'asi0') as asi0:
  asi0.run()
  E0 = asi0.total_energy
  ch0 = asi0.atomic_charges
  print ('F0=',np.round(asi0.total_forces,6) )
print ('E0=', E0)

with DFT_C_API(lib_file, initializer, MPI.COMM_WORLD, h2o1, 'asi1') as asi1:
  asi1.run()
  E01 = asi1.total_energy
  ch01 = asi1.atomic_charges
  pot1 = get_pot(asi1, h2o2.positions / units.Bohr)
  #print ('F1=',np.round(asi1.total_forces,6) )
print ('E01=', E01)

with DFT_C_API(lib_file, initializer, MPI.COMM_WORLD, h2o2, 'asi2') as asi2:
  asi2.run()
  E02 = asi2.total_energy
  ch02 = asi2.atomic_charges
  pot2 = get_pot(asi2, h2o1.positions / units.Bohr)
  #print ('F2=',np.round(asi2.total_forces,6) )
print ('E02=', E02)

#sys.exit(0)
for i in range (10):
  
  
  with DFT_C_API(lib_file, initializer, MPI.COMM_WORLD, h2o1, 'asi1') as asi1:
    asi1.register_external_potential(*pot2)
    asi1.run()
    new_pot1 = get_pot(asi1, h2o2.positions / units.Bohr)
    E1 = asi1.total_energy
    ch1 = asi1.atomic_charges
    f1=asi1.total_forces.copy()
   
  with DFT_C_API(lib_file, initializer, MPI.COMM_WORLD, h2o2, 'asi2') as asi2:
    asi2.register_external_potential(*pot1)
    asi2.run()
    new_pot2 = get_pot(asi2, h2o1.positions / units.Bohr)
    E2 = asi2.total_energy
    ch2 = asi2.atomic_charges
    f2=asi2.total_forces.copy()

  pot1 = mixer(pot1, new_pot1)
  pot2 = mixer(pot2, new_pot2)
  print (d, E0 - E01 - E02, E1 - E01, E2 - E02, 'raw:', E0, E01, E02, E1, E2, ch0, ch01, ch02, ch1, ch2)
  print (np.round(f1,6))
  print (np.round(f2,6))
  

