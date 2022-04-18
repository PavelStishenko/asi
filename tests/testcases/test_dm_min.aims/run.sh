#!/bin/sh
set -e

rm -fr work
mkdir work
cd work

cp ../*.in ./
cp ../dm.txt ./
TEST_LOG=test1
$TESTS/test_dm_min-aims.x | tee $TEST_LOG
diff --color -s $TEST_LOG ../$TEST_LOG
