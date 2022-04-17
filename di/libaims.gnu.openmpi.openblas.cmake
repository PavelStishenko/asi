set(TARGET_NAME aims CACHE STRING "")
set(CMAKE_INSTALL_PREFIX "$ENV{AIMS_HOME}" CACHE STRING "")

set(BUILD_SHARED_LIBS ON CACHE STRING "")
set(USE_MPI ON CACHE BOOL "" FORCE)
set(USE_SCALAPACK ON CACHE BOOL "" FORCE)


###############
# Basic Flags #
###############
set(CMAKE_Fortran_COMPILER mpifort CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-O2 -ffree-line-length-none" CACHE STRING "")
set(Fortran_MIN_FLAGS "-O0 -ffree-line-length-none" CACHE STRING "")
#set(LIB_PATHS "$ENV{MKLROOT}/lib/intel64" CACHE STRING "")
set(LIBS "scalapack-openmpi openblas" CACHE STRING "")

###############
# C,C++ Flags #
###############
set(CMAKE_C_COMPILER mpicc CACHE STRING "")
set(CMAKE_C_FLAGS "-O2 -DNDEBUG -funroll-loops" CACHE STRING "")
