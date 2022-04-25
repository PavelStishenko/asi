import sys, os
import numpy as np
from mpi4py import MPI
from scalapack4py import ScaLAPACK4py
from mpiprint import parprint, ordprint
from ctypes import cast, py_object, CDLL

from ase.build import molecule
from ase.io import read, write
from pyasi.asecalc import ASI_ASE_calculator
from ase import units

#os.environ["OMP_NUM_THREADS"]="1"
#os.environ["ASE_AIMS_COMMAND"]=f"mpiexec -n 4 {os.environ['HOME']}/prg/aims/build-exe/aims.220309.scalapack.mpi.x | tee -a aims.total.log.3"
#os.environ["AIMS_SPECIES_DIR"]=f"{os.environ['HOME']}/prg/aims/species_defaults/defaults_2010/tight"
#aims_lib_path = f"{os.environ['HOME']}/prg/aims/build-so-3/libaims.220309.scalapack.mpi.so"
ASI_LIB_PATH = os.environ['ASI_LIB_PATH']

sl = ScaLAPACK4py(CDLL(ASI_LIB_PATH))


def init_aims(asi):
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
    species_dir=os.environ["AIMS_SPECIES_DIR"],
    tier = [1, 2],
    density_update_method="density_matrix"
  )
  calc.write_input(asi.atoms)

def dm_init(aux, iK, iS, descr, data):
  asi = cast(aux, py_object).value
  try:
    is_distributed=False
    if descr:
      descr = sl.wrap_blacs_desc(descr)
      if descr.is_distributed:
        is_distributed = True
        is_root = (descr.myrow==0 and descr.mycol==0)
    
    locshape = (descr.locrow, descr.loccol) if is_distributed else (asi.n_basis,asi.n_basis)

    data = np.ctypeslib.as_array(data, shape=locshape).T
    
    #if not is_distributed or is_root: TODO 
    predicted_dm = predict_dm(asi.atoms)
    if is_distributed and not is_root:
      predicted_dm = None

    if is_distributed:
      sl.scatter(predicted_dm, descr, data)
    else:
      data[:, :] = predicted_dm

    #parprint ("dm predict done")
  except Exception as eee:
    print ("Something happened in dm_init", eee)


def dm_calc(aux, iK, iS, descr, data):
  asi = cast(aux, py_object).value
  asi.scf_cnt += 1
  ordprint(f"asi.scf_cnt = {asi.scf_cnt}")
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
    parprint (f"S*D = {np.sum(data * asi.overlap):.6f} E = {atoms.calc.asi.total_energy * units.Hartree:.6f}")
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
  except Exception as eee:
    print ("Something happened in dm_calc", eee)


atoms = molecule('H2O')

atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_aims, None, atoms)
atoms.calc.asi.scf_cnt = 0
atoms.calc.asi.register_overlap_callback(overlap_calc, atoms.calc.asi)
atoms.calc.asi.register_dm_callback(dm_calc, atoms.calc.asi)

parprint(f'E = {atoms.get_potential_energy():.6f}')

