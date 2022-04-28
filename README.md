# Atomic Similation Interface (ASI) API

Atomic Simulation Interface (ASI) is a native C-style API that includes functions for export and import of data structures that are used in electronic structure calculations and for classical molecular dynamics simulations. ASI aims to be a uniform, generic and efficient interface for connecting various computational chemistry and materials science codes in multiscale simulation workflows, such as QM/MM, QM/ML, QM/QM. ASI specifies functions, data types and calling conventions for export and import of density matrices, overlap and Hamiltonian matrices, electrostatic potential, atomic coordinates, charges, total energy and forces. 

ASI API is specified as a C header file `asi.h`.

[ASI API specification](https://pvst.gitlab.io/asi/asi_8h.html).

## Supported in:

* [DFTB+](https://dftbplus.org/): [in separate branch](https://github.com/PavelStishenko/dftbplus/tree/api-dm-3)
* [FHI-aims](https://fhi-aims.org/): in main branch.


## Building

* FHI-aims has embedded support of ASI API. Just build latest version of FHI-aims as a shared library and use with your code.
* DFTB+: run `make && make install` from the project root folder. See header of the `Makefile` for environment variables controlling the build process.
