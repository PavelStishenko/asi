#!/bin/sh
set -e
ulimit -s unlimited

rm -rf work
mkdir work
cd work

TEST_LOG=test1
$MPIEXEC -n 1 $TESTS/test_esp-aims.x | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG


TEST_LOG=test2
$MPIEXEC -n 2 $TESTS/test_esp-aims.x | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG


