#!/bin/sh

set -e
ulimit -s unlimited

mkdir -p asi.temp

$MPIEXEC -n 1 $TESTING_PYTHON -u $PYTESTS/test_expdmhs_default.py | tee asi.temp/test1.log
diff --color -s asi.temp/test1.log test1.log


