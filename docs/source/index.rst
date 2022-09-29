.. asi4py documentation master file, created by
   sphinx-quickstart on Mon Sep 19 21:13:35 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

asi4py is a Python wrapper for ASI API
======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   reference

Usage example
==================================

.. code-block:: python

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
      calc.write_input(asi.atoms, properties=['forces'])


  atoms = molecule('H2O')

  atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_via_ase, None, atoms)
  atoms.calc.asi.keep_density_matrix = True
  atoms.calc.asi.keep_hamiltonian = True
  atoms.calc.asi.keep_overlap = True

  parprint(f'E = {atoms.get_potential_energy():.6f}')

  S = atoms.calc.asi.overlap_storage[(1,1)]
  H = atoms.calc.asi.hamiltonian_storage[(1,1)]
  DM = atoms.calc.asi.dm_storage.get((1,1), None)
  S_cnt = atoms.calc.asi.overlap_calc_cnt[(1,1)]
  H_cnt = atoms.calc.asi.hamiltonian_calc_cnt[(1,1)]
  DM_cnt = atoms.calc.asi.dm_calc_cnt[(1,1)]

  if DM is not None:
    print(f'Nel = {np.sum(S*DM):.6f}')
    print(f'EigSum = {np.sum(H*DM):.6f}')

  parprint(f'S_cnt = {S_cnt}')
  parprint(f'H_cnt = {H_cnt}')
  parprint(f'DM_cnt = {DM_cnt}')



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
