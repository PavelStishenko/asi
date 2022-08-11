#!/bin/sh

set -e

mkdir -p asi.temp

$MPIEXEC -n 1 $TESTING_PYTHON -u $PYTESTS/test_expdmhs_esp.py OH | tee asi.temp/test1.log
grep -e 'End self' asi.temp/asi.log | wc -l | tee -a asi.temp/test1.log

