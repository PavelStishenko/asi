#!/bin/sh

set -e

mkdir -p asi.temp
export OMP_NUM_THREADS=4
$MPIEXEC -n 1 $TESTING_PYTHON -u $PYTESTS/test_sc.py | tee asi.temp/test1.log
diff --color -s asi.temp/test1.log test1.log

