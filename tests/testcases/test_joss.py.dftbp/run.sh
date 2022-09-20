#!/bin/sh

set -e

mkdir -p asi.temp

$TESTING_PYTHON -u $PYTESTS/test_joss.py | tee asi.temp/test1.log
diff --color -s asi.temp/test1.log test1.log

