#!/bin/sh

set -e
ulimit -s unlimited

mkdir -p asi.temp

$MPIEXEC -n 2 $TESTING_PYTHON -u $PYTESTS/test_e.py | tee asi.temp/test1.log

