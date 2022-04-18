#!/bin/sh
set -e

rm -rf work
mkdir work
cp *.hsd *.skf geo_end.* work/
cd work


$MPIEXEC -n 1 $TESTS/test_expdmhs-dftbp.x dftb_in.2h2o.hsd > test.2h2o.log
diff --color -s test.2h2o.log ../test.2h2o.log

$MPIEXEC -n 1 $TESTS/test_expdmhs-dftbp.x dftb_in.h2o.hsd > test.h2o.log
diff --color -s test.h2o.log ../test.h2o.log

$MPIEXEC -n 1 $TESTS/test_expdmhs-dftbp.x dftb_in.Si2.hsd > test1.log
diff --color -s test1.log ../test1.log

$MPIEXEC -n 1 $TESTS/test_expdmhs-dftbp.x dftb_in.Si2.2_kpts.hsd > test2.log
diff --color -s test2.log ../test2.log

