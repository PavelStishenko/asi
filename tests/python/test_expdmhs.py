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
    calc = Aims(xc='pbe', 
      relativistic="atomic_zora scalar",
      occupation_type="gaussian 0.010",
      sc_accuracy_eev=1E-3,
      sc_accuracy_rho=1e-05,
      sc_accuracy_etot=1e-06,
      #sc_accuracy_forces=1e-1, # just to enable force calculation
      sc_iter_limit = 3,
      postprocess_anyway = True,
      #species_dir=os.environ["AIMS_SPECIES_DIR"],
      tier = [1, 2],
      density_update_method="density_matrix"
    )
    calc.write_input(asi.atoms)
else:
  def init_via_ase(asi):
    from ase.calculators.dftb import Dftb
    calc = Dftb(label='Some_cluster',
          Hamiltonian_SCC='Yes',
          Hamiltonian_MaxAngularMomentum_='',
          Hamiltonian_MaxAngularMomentum_O='"p"',
          Hamiltonian_MaxAngularMomentum_H='"s"')
    calc.write_input(asi.atoms) # , properties=['forces']

def utriang2herm(X):
  i_lower = np.tril_indices(X.shape[0], -1)
  X = X.copy()
  X[i_lower] = X.conj().T[i_lower]
  return X

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
    asi.dm = utriang2herm(data)
    #parprint("asi.dm")
    #np.savetxt(sys.stdout, asi.dm, fmt='%10.5f')
    E = atoms.calc.asi.total_energy if asilib.ASI_flavour() == 1 else 0.0 # total_energy calls DM calculation in DFTB+ causing infinite recursion
    parprint (f"{asi.scf_cnt} S*D = {np.sum(asi.dm * asi.overlap):.6f} E = {E * units.Hartree:.6f}")
    if hasattr(asi,"hamiltonian"):
      parprint (f"{asi.scf_cnt} H*D = {np.sum(asi.dm * asi.hamiltonian):.8f}")
  except Exception as eee:
    print ("Something happened in dm_calc", eee)

def overlap_calc(aux, iK, iS, descr, data):
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
    asi.overlap = data.copy()
    #parprint("asi.overlap")
    #np.savetxt(sys.stdout, asi.overlap, fmt='%10.5f')
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
    #parprint("asi.hamiltonian")
    #np.savetxt(sys.stdout, asi.hamiltonian, fmt='%10.5f')
  except Exception as eee:
    print ("Something happened in dm_calc", eee)

atoms = molecule('H2O')

atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_via_ase, None, atoms)
atoms.calc.asi.scf_cnt = 0
atoms.calc.asi.register_overlap_callback(overlap_calc, atoms.calc.asi)
atoms.calc.asi.register_dm_callback(dm_calc, atoms.calc.asi)
atoms.calc.asi.register_hamiltonian_callback(hamiltonian_calc, atoms.calc.asi)

parprint(f'E = {atoms.get_potential_energy():.6f}')

