import sys, os
import numpy as np
from mpi4py import MPI
from mpiprint import parprint, ordprint
from ctypes import cast, py_object, CDLL

from ase.build import molecule
from ase.io import read, write
from asi4py.asecalc import ASI_ASE_calculator
from ase.optimize.lbfgs import LBFGS
from ase import units

cwd = os.getcwd()

try:
  import dmpredictXXX # provoke
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


def make_aims_calc(iPI=False):
  from ase.calculators.aims import Aims
  calc = Aims(xc='pbe', 
    relativistic="atomic_zora scalar",
    occupation_type="gaussian 0.010",
    sc_accuracy_eev=1E-3,
    sc_accuracy_rho=1e-05,
    sc_accuracy_etot=1e-06,

    sc_accuracy_forces=1e-2, # just to enable force calculation
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


def dyn_step(asi):
  #atoms.calc.asi.register_DM_init(dm_init, atoms.calc.asi) # reset dm_init callback
  atoms.calc.asi.init_density_matrix = {(1,1):predict_dm(atoms)}
  fmax = max(np.linalg.norm(atoms.get_forces(), axis=-1))
  parprint(f'E = {atoms.get_potential_energy():.6f} fmax = {fmax:.4f}' )


atoms = read(sys.argv[1])


##
## Select calculator here:
##
if True:
  atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_aims, None, atoms)
  atoms.calc.asi.init_density_matrix = {(1,1):predict_dm(atoms)}
elif False:
  atoms.calc = make_socket_calc()
else:
  atoms.calc = make_aims_calc()

parprint(f'E = {atoms.get_potential_energy():.6f}')
with np.printoptions(formatter={"float_kind":lambda x: f"{(0.0 if abs(x) < 1e-6 else x):.4f}"}):
  parprint(atoms.get_forces())

dyn = LBFGS(atoms, logfile='asi.temp/dyn.log')
dyn.attach(lambda :dyn_step(atoms.calc.asi), interval=1)
dyn.run(steps=3, fmax=1e-6)



