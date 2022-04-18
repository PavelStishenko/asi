#!/bin/sh

ulimit -s unlimited

rm -rf work
mkdir work
cd work

# H2O
cp ../geometry.in.h2o ./geometry.in


TEST_LOG=test.h2o.mixed
echo $TEST_LOG
cp -v ../control.in.h2o.mixed ./control.in
cp -v ../dm_1_1.init ./
$TESTS/test_initdm-aims.x | tee  $TEST_LOG
grep -e 'End self-consistency iteration #' asi.log  | tail -n 1 | tee -a $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG


TEST_LOG=test.h2o.mixed.forces
echo $TEST_LOG
cp -v ../control.in.h2o.mixed.forces ./control.in
cp -v ../dm_1_1.init ./
mpiexec -n 2 $TESTS/test_initdm-aims.x | tee  $TEST_LOG
grep -e 'End self-consistency iteration #' asi.log  | tail -n 1 | tee -a $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG




