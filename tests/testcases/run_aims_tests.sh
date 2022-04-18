#!/bin/sh
set -e

ulimit -s unlimited

export TESTS="${TESTS:-${PWD}/../build/}"

echo "Tests dir:" $TESTS

for t in test_*.aims 
do
  echo "Run test " $t
  cd $t
  ./run.sh
  cd ..
done
