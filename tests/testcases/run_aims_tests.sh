#!/bin/sh
set -e

ulimit -s unlimited

export TESTS="${TESTS:-${PWD}/../build/}"
export MPIEXEC="${MPIEXEC:-mpiexec}"
export TESTING_PYTHON="${TESTING_PYTHON:-python3}"
export PYTESTS="${PYTESTS:-${PWD}/../python/}"
#export PYTHONPATH=$PYTHONPATH:${PWD}/../../python
#export DMPREDICT_PATH=${PWD}/../../../dmpredict/
#export PYTHONPATH=$PYTHONPATH:${DMPREDICT_PATH}


echo "Tests dir:" $TESTS

env

for t in test_*.aims
do
  echo "Run test " $t
  cd $t
  ./run.sh
  cd ..
done
