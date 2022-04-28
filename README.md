# Atomic Similation Interface (ASI) API

[ASI API specification](https://pvst.gitlab.io/asi/asi_8h.html).

##Supported in:

* [DFTB+](https://dftbplus.org/): [in separate branch](https://github.com/PavelStishenko/dftbplus/tree/api-dm-3)
* [FHI-aims](https://fhi-aims.org/): in main branch.


##Building

* FHI-aims has embedded support of ASI API. Just build latest version of FHI-aims as a shared library and use with your code.
* DFTB+: run `make && make install` from the project root folder. See header of the `Makefile` for environment variables controlling the build process.
