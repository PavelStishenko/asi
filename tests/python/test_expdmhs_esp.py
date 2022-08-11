import sys, os
import numpy as np, scipy
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
    calc = Aims(xc='pbe', 
      relativistic="atomic_zora scalar",
      occupation_type="gaussian 0.010",
      sc_accuracy_eev=1E-6,
      sc_accuracy_rho=1e-05,
      sc_accuracy_etot=1e-06,
      #sc_accuracy_forces=1e-1, # just to enable force calculation
      sc_iter_limit = 20,
      postprocess_anyway = True,
      species_dir=os.environ["AIMS_SPECIES_DIR"] + '/../really_tight/',
      tier = [1,1],
      density_update_method="density_matrix",
      #homogeneous_field='1 1 1',
    )
    calc.write_input(asi.atoms)
else:
  def init_via_ase(asi):
    from ase.calculators.dftb import Dftb
    calc = Dftb(label='Some_cluster',
          Hamiltonian_SCC='Yes',
          #Hamiltonian_Charge='-2',
          Hamiltonian_MaxAngularMomentum_='',
          Hamiltonian_MaxAngularMomentum_O='"f"',
          Hamiltonian_MaxAngularMomentum_H='"p"',
          #Hamiltonian_ElectricField_='',
          #Hamiltonian_ElectricField_External_='',
          #Hamiltonian_ElectricField_External_Strength='100.0',
          #Hamiltonian_ElectricField_External_Direction='1.0 2.0 3.0',
          #Hamiltonian_ElectricField_PointCharges_='',
          #Hamiltonian_ElectricField_PointCharges_CoordsAndCharges_='',
          #Hamiltonian_ElectricField_PointCharges_CoordsAndCharges_empty='  2 0 0 1\n  0 2 0 1.1\n  0 0 2 1.2'
          )
    calc.write_input(asi.atoms, properties=['forces'])

def dm_calc(aux, iK, iS, descr, data):
  asi = cast(aux, py_object).value
  asi.scf_cnt += 1
  try:
    if descr:
      descr = sl.wrap_blacs_desc(descr)
      if descr.is_distributed:
        parprint("distributed case not implemented")
        return # TODO distributed case not implemented
      else:
        pass
    else:
      pass
    # single process case:
    #print (f"dm_calc invoked {asi.scf_cnt}")
    data = np.ctypeslib.as_array(data, shape=(asi.n_basis,asi.n_basis))
    E = atoms.calc.asi.total_energy if asilib.ASI_flavour() == 1 else 0.0 # total_energy calls DM calculation in DFTB+ causing infinite recursion
    parprint (f"{asi.scf_cnt} S*D = {np.sum(data * asi.overlap):.6f} E = {E * units.Hartree:.6f}")
    if hasattr(asi,"hamiltonian"):
      parprint (f"{asi.scf_cnt} H*D = {np.sum(data * asi.hamiltonian):.8f}")
    asi.dm = data.copy()
  except Exception as eee:
    print ("Something happened in dm_calc", eee)

def overlap_calc(aux, iK, iS, descr, data):
  asi = cast(aux, py_object).value
  parprint (f"overlap_calc invoked {asi.scf_cnt}")
  try:
    if descr:
      descr = sl.wrap_blacs_desc(descr)
      if descr.is_distributed:
        parprint("distributed case not implemented")
        return # TODO distributed case not implemented
      else:
        pass
    else:
      pass
    # single process case:
    data = np.ctypeslib.as_array(data, shape=(asi.n_basis,asi.n_basis))
    asi.overlap = data.copy()
  except Exception as eee:
    print ("Something happened in dm_calc", eee)

def hamiltonian_calc(aux, iK, iS, descr, data):
  asi = cast(aux, py_object).value
  try:
    if descr:
      descr = sl.wrap_blacs_desc(descr)
      if descr.is_distributed:
        parprint("distributed case not implemented")
        return # TODO distributed case not implemented
      else:
        pass
    else:
      pass
    # single process case:
    #print (f"dm_calc invoked {asi.scf_cnt}")
    data = np.ctypeslib.as_array(data, shape=(asi.n_basis,asi.n_basis))
    asi.hamiltonian = data.copy()
  except Exception as eee:
    print ("Something happened in dm_calc", eee)
  #return np.loadtxt(f"{cwd}/dm_1_1.init.full").T
  
def dm_init(aux, iK, iS, descr, data):
  asi = cast(aux, py_object).value
  parprint("dm_init")
  try:
    data = np.ctypeslib.as_array(data, shape=(asi.n_basis,asi.n_basis)).T
    data[:, :] = asi.predicted_dm
  except Exception as eee:
    print ("Something happened in dm_init", eee)

pnt_charges = np.array([0.1, 0.1, 0.1], dtype=np.float64)
def pnt_charge_esp(coords):
  pnt_charge_coord = np.array([[5, 0, 0.], [0, 5, 0.], [0, 0, 5.]], dtype=np.float64)
  res = pnt_charges[:, None] / np.linalg.norm(coords[None, :, :] - pnt_charge_coord[:,None,:], axis=-1)
  parprint ('res.shape', res.shape, 'coords.shape', coords.shape)
  return np.sum(res, axis=0)

def homog_esp(coords):
  direction = np.array([1, 0, 0], dtype=np.float64)
  strength = 1.0
  direction /= np.linalg.norm(direction)
  res = (coords @ direction) * strength
  parprint ('res.shape', res.shape, 'coords.shape', coords.shape)
  return res

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
  pot_func = pnt_charge_esp
  #pot_func = homog_esp
  if (potential):
    potential = np.ctypeslib.as_array(potential, shape=(n, ))
    potential[:] = pot_func(coords)
    #print('esp potential:')
    #np.savetxt(sys.stdout, potential)
    
  if (potential_grad):
    potential_grad = np.ctypeslib.as_array(potential_grad, shape=(n, 3))
    potential_grad[:,:] = pot_grad(coords, pot_func, 0.0001)
    #print('esp potential_grad:')
    #np.savetxt(sys.stdout, potential_grad)

mol_name = 'H2O' if len(sys.argv)<2 else sys.argv[1]
if '.traj' in mol_name:
  idx = int(sys.argv[2])
  atoms = read(mol_name, index=idx)
  mol_name = f'{mol_name}.{idx}'
else:
  atoms = molecule(mol_name)
#atoms1.rattle() = atoms.copy()
parprint('atoms.positions:\n', atoms.positions)
#atoms.positions = atoms.positions @ [[0, 0, 1],[0, 1, 0],[1, 0, 0]]
parprint('atoms.positions:\n', atoms.positions)
#atoms.info['homogeneous_field'] = [1.0, 1.0, 2.0]

calc = ASI_ASE_calculator(ASI_LIB_PATH, init_via_ase, None, atoms)
for i,q in enumerate([[0, 0.0, 0.0], [0, 0.0, 0.0],]):
  pnt_charges[:] = q[:]
  calc.atoms = None
  atoms.calc = calc
  atoms.calc.asi.scf_cnt = 0
  atoms.calc.asi.register_overlap_callback(overlap_calc, atoms.calc.asi)
  atoms.calc.asi.register_dm_callback(dm_calc, atoms.calc.asi)
  atoms.calc.asi.register_hamiltonian_callback(hamiltonian_calc, atoms.calc.asi)
  #atoms.calc.asi.register_DM_init(dm_init, atoms.calc.asi)
  atoms.calc.asi.register_external_potential(esp, atoms)
  #atoms.calc.asi.predicted_dm = np.loadtxt('dm2.txt')

  parprint(f'E = {atoms.get_potential_energy():.6f}')
  np.savetxt(f'D_{mol_name}.txt', atoms.calc.asi.dm)
  np.savetxt(f'H_{mol_name}.txt', atoms.calc.asi.hamiltonian)
  np.savetxt(f'S_{mol_name}.txt', atoms.calc.asi.overlap)
  L,C = scipy.linalg.eigh(atoms.calc.asi.hamiltonian, atoms.calc.asi.overlap)
  np.savetxt(f'L_{mol_name}_{i}.txt', L)
  np.savetxt(f'C_{mol_name}_{i}.txt', C)

