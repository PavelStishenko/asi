#!/bin/sh
set -e

rm -rf work
mkdir work
cp *.skf work/
cd work

TEST_LOG=test1
$MPIEXEC -n 2 $TESTS/test_esp-dftbp.x | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG


