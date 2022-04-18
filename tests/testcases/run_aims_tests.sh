#!/bin/sh
set -e

ulimit -s unlimited

export TESTS="${TESTS:-${PWD}/../build/}"
export MPIEXEC="${MPIEXEC:-mpiexec}"

echo "Tests dir:" $TESTS

env

for t in test_*.aims 
do
  echo "Run test " $t
  cd $t
  ./run.sh
  cd ..
done
