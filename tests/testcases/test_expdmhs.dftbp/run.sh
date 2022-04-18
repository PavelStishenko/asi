#!/bin/sh

rm -r work
mkdir work
cp *.hsd *.skf geo_end.* work/
cd work

mpiexec -n 1 ../../test_expdmhs-dftbp dftb_in.2h2o.hsd > test.2h2o.log
diff --color -s test.2h2o.log ../test.2h2o.log

mpiexec -n 1 ../../test_expdmhs-dftbp dftb_in.h2o.hsd > test.h2o.log
diff --color -s test.h2o.log ../test.h2o.log

mpiexec -n 1 ../../test_expdmhs-dftbp dftb_in.Si2.hsd > test1.log
diff --color -s test1.log ../test1.log

mpiexec -n 1 ../../test_expdmhs-dftbp dftb_in.Si2.2_kpts.hsd > test2.log
diff --color -s test2.log ../test2.log

