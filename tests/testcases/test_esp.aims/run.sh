#!/bin/sh
set -e

rm -rf work
mkdir work
cd work

TEST_LOG=test1
mpiexec -n 1 $TESTS/test_esp-aims.x | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG


TEST_LOG=test2
mpiexec -n 2 $TESTS/test_esp-aims.x | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG


