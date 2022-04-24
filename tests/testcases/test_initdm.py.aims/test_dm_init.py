import sys, os
import numpy as np
from mpi4py import MPI
from scalapack4py import ScaLAPACK4py
from mpiprint import parprint, ordprint
from ctypes import cast, py_object, CDLL

from ase.build import molecule
from ase.io import read, write
from pyasi.asecalc import ASI_ASE_calculator
from ase.optimize.lbfgs import LBFGS
from ase import units

cwd = os.getcwd()


try:
  import dmpredict
  def predict_water_dm(mol):
    dmpredict_path = os.environ.get("DMPREDICT_PATH", cwd)
    return dmpredict.predict(dmpredict_path, mol, )
except ModuleNotFoundError:
  parprint("dmpredict package is not available. Loading ground state DM as a stub.")
  ground_state_dm = np.loadtxt(f"{cwd}/dm_1_1.init")
  def predict_water_dm(mol):
    return ground_state_dm
  

#os.environ["OMP_NUM_THREADS"]="1"
#os.environ["ASE_AIMS_COMMAND"]=f"mpiexec -n 4 {os.environ['HOME']}/prg/aims/build-exe/aims.220309.scalapack.mpi.x | tee -a aims.total.log.3"
#os.environ["AIMS_SPECIES_DIR"]=f"{os.environ['HOME']}/prg/aims/species_defaults/defaults_2010/tight"
#aims_lib_path = f"{os.environ['HOME']}/prg/aims/build-so-3/libaims.220309.scalapack.mpi.so"
ASI_LIB_PATH = os.environ['ASI_LIB_PATH']

sl = ScaLAPACK4py(CDLL(ASI_LIB_PATH))

def make_aims_calc(iPI=False):
  from ase.calculators.aims import Aims
  calc = Aims(xc='pbe', 
    relativistic="atomic_zora scalar",
    occupation_type="gaussian 0.010",
    sc_accuracy_eev=1E-3,
    sc_accuracy_rho=1e-05,
    sc_accuracy_etot=1e-06,

    sc_accuracy_forces=1e-1, # just to enable force calculation
    species_dir=os.environ["AIMS_SPECIES_DIR"],
    tier = [1, 2],
    #density_update_method="density_matrix"
  )
  if iPI:
    calc.parameters["use_pimd_wrapper"]=("localhost", 12346)
  return calc

def make_socket_calc():
  from ase.calculators.socketio import SocketIOCalculator
  aims_calc = make_aims_calc(True)
  return SocketIOCalculator(aims_calc, log="socket_calc.log", port=12346)

def init_aims(asi):
  make_aims_calc().write_input(asi.atoms)


def predict_dm(water_clusters):
  assert len(water_clusters) % 3 == 0
  n_mols = len(water_clusters) // 3
  total_dm = np.ndarray(shape=(n_mols*44, n_mols*44), order='F')
  total_dm[:,:] = 0
  for i in range(n_mols):
    mol = water_clusters[i*3:(i+1)*3]
    mol_dm = predict_water_dm(mol)
    total_dm[i*44:(i+1)*44, i*44:(i+1)*44] = mol_dm
  return total_dm
  #return np.loadtxt(f"{cwd}/dm_1_1.init.full").T
  
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
    np.savetxt(f"dm_{asi.scf_cnt}.txt", data)
    #print (f"S*D = {np.sum(data * asi.overlap)}")
  except Exception as eee:
    print ("Something happened in dm_calc", eee)

def dyn_step(asi):
  atoms.calc.asi.register_DM_init(dm_init, atoms.calc.asi) # reset dm_init callback
  fmax = max(np.linalg.norm(atoms.get_forces(), axis=-1))
  parprint(f'E = {atoms.get_potential_energy():.6f} fmax = {fmax:.6f}' )


atoms = read(sys.argv[1])


##
## Select calculator here:
##
if True:
  atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_aims, None, atoms)
  atoms.calc.asi.register_DM_init(dm_init, atoms.calc.asi)
  atoms.calc.asi.scf_cnt = 0
  #atoms.calc.asi.register_dm_callback(dm_calc, atoms.calc.asi)
elif False:
  atoms.calc = make_socket_calc()
else:
  atoms.calc = make_aims_calc()

parprint(f'E = {atoms.get_potential_energy():.6f}')
np.set_printoptions(precision=6, suppress=True)
parprint(atoms.get_forces())

dyn = LBFGS(atoms, logfile='asi.temp/dyn.log')
dyn.attach(lambda :dyn_step(atoms.calc.asi), interval=1)
dyn.run(steps=3, fmax=1e-6)



