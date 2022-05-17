#!/bin/sh

$MPIEXEC -n 2 $TESTING_PYTHON -u -c 'from mpi4py import MPI; print (MPI);'

$MPIEXEC -n 2 $TESTING_PYTHON -u -c 'from mpi4py import MPI; print (MPI.COMM_WORLD.Get_rank(), " of ", MPI.COMM_WORLD.Get_size());'
