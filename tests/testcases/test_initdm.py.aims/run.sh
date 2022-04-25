#!/bin/sh

set -e

mkdir -p asi.temp

$MPIEXEC -n 4 $TESTING_PYTHON -u $PYTESTS/test_dm_init.py h2o-1.xyz | tee asi.temp/test1.log
grep -e 'End self' asi.temp/asi.log | wc -l | tee -a asi.temp/test1.log
diff --color -s asi.temp/test1.log test1.log

$MPIEXEC -n 4 $TESTING_PYTHON -u $PYTESTS/test_dm_init.py TIP4P-2.xyz | tee asi.temp/test2.log
grep -e 'End self' asi.temp/asi.log | wc -l | tee -a asi.temp/test2.log
diff --color -s asi.temp/test2.log test2.log

