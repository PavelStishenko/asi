import sys, os
import numpy as np
from mpi4py import MPI
from scalapack4py import ScaLAPACK4py
from mpiprint import parprint, ordprint
from ctypes import cast, py_object, CDLL, RTLD_GLOBAL

from ase.build import molecule
from ase.io import read, write
from asi4py.asecalc import ASI_ASE_calculator
from ase import units

ASI_LIB_PATH = os.environ['ASI_LIB_PATH']
asilib = CDLL(ASI_LIB_PATH, mode=RTLD_GLOBAL)
sl = ScaLAPACK4py(asilib)


if asilib.ASI_flavour() == 1:
  def init_via_ase(asi):
    from ase.calculators.aims import Aims
    calc = Aims(xc='pw-lda', 
      sc_accuracy_eev=1E-4,
      sc_accuracy_rho=1E-4,
      sc_accuracy_etot=1E-06,
      sc_accuracy_forces=5e-2, 
      sc_iter_limit = 100,
      postprocess_anyway = True,
      species_dir=os.environ["AIMS_SPECIES_DIR"] + '/../light',
      compensate_multipole_errors=True,
      spin='collinear',
      multiplicity=3,
      final_forces_cleaned=False,
      output=['hirshfeld-I',],
      density_update_method="density_matrix",
    )
    asi.atoms.set_initial_magnetic_moments([2,0,0])
    calc.write_input(asi.atoms)
else:
  def init_via_ase(asi):
    from ase.calculators.dftb import Dftb
    calc = Dftb(label='Some_cluster',
          Hamiltonian_SCC='Yes',
          Hamiltonian_MaxAngularMomentum_='',
          Hamiltonian_MaxAngularMomentum_O='"p"',
          Hamiltonian_MaxAngularMomentum_H='"s"')
    calc.write_input(asi.atoms, properties=['forces'])

def pnt_charge_esp(coords):
  pnt_charge = 1.0
  pnt_charge_coord = np.array([1, 0, 2.0])
  return pnt_charge / np.linalg.norm(coords - pnt_charge_coord, axis=1)

def pot_grad(coords, pot_func, d):
  p0 = pot_func(coords)
  px = pot_func(coords + [d, 0, 0])
  py = pot_func(coords + [0, d, 0])
  pz = pot_func(coords + [0, 0, d])
  res = np.zeros(shape=(coords.shape[0], 3), dtype=p0.dtype)
  res[:,0] = (px - p0) / d
  res[:,1] = (py - p0) / d
  res[:,2] = (pz - p0) / d
  return res

def esp(aux, n, coords, potential, potential_grad):
  coords = np.ctypeslib.as_array(coords, shape=(n, 3))
  #print('esp coords:')
  #np.savetxt(sys.stdout, coords)
  if (potential):
    potential = np.ctypeslib.as_array(potential, shape=(n, ))
    potential[:] = pnt_charge_esp(coords)
    #print('esp potential:')
    #np.savetxt(sys.stdout, potential)
    
  if (potential_grad):
    potential_grad = np.ctypeslib.as_array(potential_grad, shape=(n, 3))
    potential_grad[:,:] = pot_grad(coords, pnt_charge_esp, 0.0001)
    #print('esp potential_grad:')
    #np.savetxt(sys.stdout, potential_grad)

atoms = molecule('H2O')
atoms.positions = [[1.0000000011, 0.0000000000, -1.0000000011],   [1.0000000011,   0.7830640008,   0.0000000000],   [1.0000000011,  -0.7830640008,   0.0000000000]] # postions from test_esp.cpp


atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_via_ase, None, atoms)
atoms.calc.asi.register_external_potential(esp, atoms)
parprint(f'E = {atoms.get_potential_energy():.6f} eV')
#parprint(f'E = {atoms.get_potential_energy()*0.9999999439805366:.6f}') # AIMS units conversion to correspond test_esp.cpp
parprint(f'E = {atoms.get_potential_energy()/units.Ha:.6f} Ha')
if MPI.COMM_WORLD.rank == 0:
  print('Forces [Ha/Bohr]:')
  np.savetxt(sys.stdout, atoms.get_forces()/units.Ha*units.Bohr, fmt='%10.6f')
parprint(f'atomic_charges = {atoms.calc.asi.atomic_charges}')

