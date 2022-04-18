#!/bin/sh

ulimit -s unlimited

rm -rf work
mkdir work
cd work


# O2
cp ../geometry.in.o2 ./geometry.in


cp -v ../control.in.o2 ./control.in
TEST_LOG=test10
echo $TEST_LOG
$MPIEXEC -n 1 $TESTS/test_expdmhs-aims.x dummy-arg | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG

# H2O
cp ../geometry.in.h2o ./geometry.in

cp -v ../control.in.h2o ./control.in
TEST_LOG=test9
echo $TEST_LOG
$MPIEXEC -n 1 $TESTS/test_expdmhs-aims.x dummy-arg | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG

cp -v ../control.in.h2o.mixed ./control.in
TEST_LOG=test.h2o.mixed
echo $TEST_LOG
$MPIEXEC -n 4 $TESTS/test_expdmhs-aims.x dummy-arg | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG


# Si2
cp ../geometry.in ./geometry.in

cp -v ../control.in.1_1_r ./control.in
TEST_LOG=test1
echo $TEST_LOG
$MPIEXEC -n 1 $TESTS/test_expdmhs-aims.x dummy-arg | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG

cp -v ../control.in.1_1_r ./control.in
TEST_LOG=test2
echo $TEST_LOG
$MPIEXEC -n 4 $TESTS/test_expdmhs-aims.x dummy-arg | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG

cp -v ../control.in.2_2_c ./control.in
TEST_LOG=test6
echo $TEST_LOG
$MPIEXEC -n 4 $TESTS/test_expdmhs-aims.x dummy-arg | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG

cp -v ../control.in.1_1_c ./control.in
TEST_LOG=test3
echo $TEST_LOG
$MPIEXEC -n 4 $TESTS/test_expdmhs-aims.x dummy-arg | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG

cp -v ../control.in.4_1_c ./control.in
TEST_LOG=test7
echo $TEST_LOG
$MPIEXEC -n 4 $TESTS/test_expdmhs-aims.x dummy-arg | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG

cp -v ../control.in.2_1_c ./control.in
TEST_LOG=test4
echo $TEST_LOG
$MPIEXEC -n 2 $TESTS/test_expdmhs-aims.x dummy-arg | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG

cp -v ../control.in.2_1_c ./control.in
TEST_LOG=test5
echo $TEST_LOG
$MPIEXEC -n 4 $TESTS/test_expdmhs-aims.x dummy-arg | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG

