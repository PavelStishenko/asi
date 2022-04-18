#!/bin/sh

ulimit -s unlimited

for t in test_*.dftbp
do
  cd $t
  ./run.sh
  cd ..
done
