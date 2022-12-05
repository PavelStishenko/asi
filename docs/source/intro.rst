==================================================================================================
Introduction
==================================================================================================


ASI API specification
=======================

ASI API is specified as a C header file `asi.h`_. Codes implementing ASI API must provide linkable library with definitions of functions from `asi.h`_. Depending on particular usage of the implementaions, some functions can be ommited or implemented as stubs, if they are not going to used. To use Python ASI wrapper it is necessary to have all functions from `asi.h` defined, but of course stub definitions can be used.

.. _`asi.h`: https://pvst.gitlab.io/asi/asi_8h.html

Implementations
=====================

* `DFTB+ <https://dftbplus.org/>`_: in the separate branch `api-H-import <https://github.com/PavelStishenko/dftbplus/tree/api-H-import>`_: `in separate branch <https://github.com/PavelStishenko/dftbplus/tree/api-dm-3>`_. Milestone release: `23.1`.

* `FHI-aims <https://fhi-aims.org/>`_: in the main branch.


Building
==========

FHI-aims
-----------

FHI-aims has embedded support of ASI API. Just build latest version of FHI-aims as a shared library and use with your code.


DFTB+
-----------

1. Download and build DFTB+ from `the branch with ASI API <https://github.com/PavelStishenko/dftbplus/tree/api-dm-3>`_ with shared library support.

2. Set environment variables `DFTBP*INCLUDE` and `DFTBP*LIB_DIR` to folders with DFTB+ C-headers and libraries.

3. Optionally export environment variables `INSTALL*PREFIX` and `BUILD*PATH` to set installation and building locations.

4. Run `make && make install` from the root of the working copy of this repository. 

5. The shared library implementing ASI API for DFTB+ will be in `${INSTALL_PREFIX}/lib`.

Python wrapper `asi4py`
========================

A Python wrapper `asi4py` is available for ASI-enabled codes.

To install asi4py use `pip`

```
pip install asi4py
```

or install from `git sources <https://gitlab.com/pvst/asi/-/tree/master/pyasi>`_


Testing
=========

Use `Makefile` in `tests` folder to build native tests. Set environment variables in the header of `tests/Makefile` to link with proper ASI API implementaions.

To run tests go to `tests/testcases` and run `run*dftbp*tests.sh` or `run*aims*tests.sh` to run test.

Usage
=======

See `tests/src` for examples of usage in native code.

See `tests/python` for examples of usage in Python scripts.



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
        #sc_accuracy_forces=1e-1, # enable forces calculation
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


