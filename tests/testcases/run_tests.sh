#!/bin/sh
set -e

export MPIEXEC="${MPIEXEC:-mpiexec}"
export TESTING_PYTHON="${TESTING_PYTHON:-python3}"

# set default path suggesting that PWD is tests/testcases
export TESTS="${TESTS:-${PWD}/../build}"
export PYTESTS="${PYTESTS:-${PWD}/../python/}"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PWD}/../../build/

env

echo "MPIEXEC:" $MPIEXEC
echo "TESTING_PYTHON:" $TESTING_PYTHON
echo "TESTS dir:" $TESTS
echo "PYTESTS dir:" $PYTESTS

for t in $@
do
  echo "Run test " $t
  cd $t
  ./run.sh
  cd ..
done
