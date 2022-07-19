from ctypes import cdll, CDLL, RTLD_LOCAL
from ctypes import POINTER, byref, c_int, c_int64, c_bool, c_char_p, c_double, c_void_p, CFUNCTYPE, py_object, cast, byref
import ctypes


import numpy as np
from numpy.ctypeslib import ndpointer

import sys
from pathlib import Path
import os, shutil
from ase import units

libdl = cdll.LoadLibrary('libdl.so')
dmhs_callback = CFUNCTYPE(None, c_void_p, c_int, c_int, POINTER(c_int), POINTER(c_double))  # void(*)(void *aux_ptr, int iK, int iS, int *blacs_descr, void *blacs_data)

class DFT_C_API:
  def __init__(self, lib_file, initializer, mpi_comm=None, atoms=None, work_dir='asi.temp', logfile='asi.log'):
    self.lib_file = Path(lib_file).resolve()
    self.initializer = initializer
    if mpi_comm is not None:
      self.mpi_comm = mpi_comm
    else:
      from mpi4py import MPI
      self.mpi_comm = MPI.COMM_WORLD
    self.atoms = atoms.copy() if atoms is not None else None
    self.work_dir = Path(work_dir)
    self.work_dir.mkdir(parents=True, exist_ok=True)
    self.logfile = logfile

  def __enter__(self):
    return self.init()

  def __exit__(self, type, value, traceback):
    #Exception handling here
    #print ("__exit__: ", type, value, traceback)
    self.close()  

  def init(self):
    curdir = os.getcwd()
    try:
      os.chdir(self.work_dir)
      
      if self.mpi_comm.Get_rank() == 0:
        self.initializer(self)
    
      # Load the FHI-aims library
      self.lib = CDLL(self.lib_file, mode=RTLD_LOCAL)

      self.lib.ASI_n_atoms.restype = c_int
      self.lib.ASI_energy.restype = c_double
      self.lib.ASI_forces.restype = POINTER(c_double)
      self.lib.ASI_atomic_charges.restype = POINTER(c_double)
      self.lib.ASI_atomic_charges.argtypes  = [c_int,]
      self.lib.ASI_calc_esp.argtypes = [c_int, ndpointer(dtype=np.float64), ndpointer(dtype=np.float64), ndpointer(dtype=np.float64)]
      self.lib.ASI_register_dm_callback.argtypes = [dmhs_callback, c_void_p]
      self.lib.ASI_register_overlap_callback.argtypes = [dmhs_callback, c_void_p]
      self.lib.ASI_register_hamiltonian_callback.argtypes = [dmhs_callback, c_void_p]
      if hasattr(self.lib, "ASI_register_dm_init_callback"):
        self.lib.ASI_register_dm_init_callback.argtypes = [dmhs_callback, c_void_p]
      self.lib.ASI_is_hamiltonian_real.restype = c_bool
      self.lib.ASI_get_basis_size.restype = c_int
      
      input_filename = {1:"dummy", 2:"dftb_in.hsd"}[self.lib.ASI_flavour()]
      self.lib.ASI_init(input_filename.encode('UTF-8'), self.logfile.encode('UTF-8'), c_int(self.mpi_comm.py2f()))
      if (self.lib.ASI_flavour() == 2):
        self.set_coords() # FIXME
      return self
    finally:
      os.chdir(curdir)
  
  def close(self):
    curdir = os.getcwd()
    try:
      os.chdir(self.work_dir)
      self.lib.ASI_finalize()
      handle = self.lib._handle
      del self.lib
      res = libdl.dlclose(handle)
      assert res == 0, "dlclose = {res}"
      if self.mpi_comm.Get_rank() == 0:
        os.system(f"cat {self.logfile} >> total.log")
    finally:
      os.chdir(curdir)
    
  def run(self):
    curdir = os.getcwd()
    try:
      os.chdir(self.work_dir)
      self.lib.ASI_run()
    finally:
      os.chdir(curdir)

  def register_DM_init(self, DM_init_callback, DM_init_aux):
    self.DM_init_callback = dmhs_callback(DM_init_callback)
    self.DM_init_aux = DM_init_aux
    self.lib.ASI_register_dm_init_callback(self.DM_init_callback, c_void_p.from_buffer(py_object(self.DM_init_aux)))

  def register_overlap_callback(self, overlap_callback, overlap_aux):
    self.overlap_callback = dmhs_callback(overlap_callback)
    self.overlap_aux = overlap_aux
    self.lib.ASI_register_overlap_callback(self.overlap_callback, c_void_p.from_buffer(py_object(self.overlap_aux)))

  def register_hamiltonian_callback(self, hamiltonian_callback, hamiltonian_aux):
    self.hamiltonian_callback = dmhs_callback(hamiltonian_callback)
    self.hamiltonian_aux = hamiltonian_aux
    self.lib.ASI_register_hamiltonian_callback(self.hamiltonian_callback, c_void_p.from_buffer(py_object(self.hamiltonian_aux)))

  def register_dm_callback(self, dm_callback, dm_aux):
    self.dm_callback = dmhs_callback(dm_callback)
    self.dm_aux = dm_aux
    self.lib.ASI_register_dm_callback(self.dm_callback, c_void_p.from_buffer(py_object(self.dm_aux)))

  def register_external_potential(self, ext_pot_func, ext_pot_aux_obj):
    '''
      self.ext_pot_func returns potential for positive charges
    '''
    self.ext_pot_func = ext_pot_func
    self.ext_pot_aux_obj = ext_pot_aux_obj
    self.lib.ASI_register_external_potential(self.ext_pot_func, c_void_p.from_buffer(py_object(self.ext_pot_aux_obj)))

  def calc_esp(self, coords):
    n = len(coords)
    esp = np.zeros((n,), dtype=c_double)
    esp_grad = np.zeros((n,3), dtype=c_double)
    self.lib.ASI_calc_esp(c_int(n), coords.ravel(), esp, esp_grad) 
    return esp, esp_grad

  @property
  def n_atoms(self):
    return self.lib.ASI_n_atoms()

  @property
  def n_basis(self):
    return self.lib.ASI_get_basis_size()
  
  @property
  def total_forces(self):
    forces_ptr = self.lib.ASI_forces()
    if forces_ptr:
      return np.ctypeslib.as_array(forces_ptr, shape=(self.n_atoms, 3))
    else:
      return None

  @property
  def atomic_charges(self):
    chg_ptr = self.lib.ASI_atomic_charges(-1)
    if chg_ptr:
      return np.ctypeslib.as_array(chg_ptr, shape=(self.n_atoms,)).copy()
    else:
      return None

  @property
  def total_energy(self):
    return self.lib.ASI_energy()
  
  def set_coords(self, coords=None):
    if coords is not None:
      assert False
      self.atoms.positions[:] = coords
    self.lib.ASI_set_atom_coords((self.atoms.positions / units.Bohr).ctypes.data_as(c_void_p), len(self.atoms))
    
  @property
  def is_eigen_real(self):
    return self.lib.ASI_is_hamiltonian_real()



