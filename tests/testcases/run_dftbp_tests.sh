#!/bin/sh
set -e

export TESTS="${TESTS:-${PWD}/../build}"
export MPIEXEC="${MPIEXEC:-mpiexec}"

echo "Tests dir:" $TESTS
ls -lAhtr $TESTS

env


for t in test_*.dftbp
do
  echo "Run test " $t
  cd $t
  ./run.sh
  cd ..
done
