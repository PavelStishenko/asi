#!/bin/sh

ulimit -s unlimited

rm -rf work
mkdir work
cp *.skf work/
cd work

TEST_LOG=test1
$MPIEXEC -n 2 ../../test_esp-dftbp | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG


