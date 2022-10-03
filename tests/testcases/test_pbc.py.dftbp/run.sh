#!/bin/sh

set -e

mkdir -p asi.temp

export DFTB_PREFIX=$DFTB_PREFIX/../pbc-0-3/

$MPIEXEC -n 1 $TESTING_PYTHON -u $PYTESTS/test_pbc.py | tee asi.temp/test1
diff --color -s asi.temp/test1 test1


