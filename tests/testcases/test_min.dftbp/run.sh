#!/bin/sh
set -e

rm -rf work
mkdir work
cp *.skf work/
cp dftb_in.hsd work/
cp geo_end.gen work/
cd work

TEST_LOG=test1
$MPIEXEC -n 1 $TESTS/test_min-dftbp.x | tee  $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG


