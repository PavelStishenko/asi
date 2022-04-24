#!/bin/sh
set -e

ulimit -s unlimited

export TESTS="${TESTS:-${PWD}/../build/}"
export MPIEXEC="${MPIEXEC:-mpiexec}"
export PYTHONPATH=$PYTHONPATH:${PWD}/../../python
export TESTING_PYTHON="${TESTING_PYTHON:-${HOME}/ve/bin/python3}"
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
