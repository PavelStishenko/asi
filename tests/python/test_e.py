import sys, os
import numpy as np
from ctypes import CDLL, RTLD_LOCAL, RTLD_GLOBAL
ASI_LIB_PATH = os.environ['ASI_LIB_PATH']
#CDLL(ASI_LIB_PATH, mode=RTLD_LOCAL)

#from mpi4py import rc as mpi4py_rc
#mpi4py_rc.initialize = False
#print ("mpi4py_rc", mpi4py_rc.initialize)
#mpi4py.rc(initialize=False, finalize=False)
from mpi4py import MPI
#print ("MPI.Is_initialized = ", MPI.Is_initialized())
#MPI.Init()
#print ("MPI.Is_initialized2 = ", MPI.Is_initialized())

from ase.build import molecule
from asi4py.asecalc import ASI_ASE_calculator
from mpiprint import parprint, ordprint

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
  )
  calc.write_input(asi.atoms)

atoms = molecule('H2O')
atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_via_ase, MPI.COMM_WORLD, atoms)
parprint(f'E = {atoms.get_potential_energy():.6f}')

