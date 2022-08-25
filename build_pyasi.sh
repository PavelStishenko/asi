#!/bin/sh

PYASI_VERSION=`python3 -c "from build import toml_loads; print(toml_loads(open('pyasi/pyproject.toml').read())['project']['version'])"`

python3 -m build pyasi

pip install --force-reinstall --no-dependencies  `ls -tr pyasi/dist/asi4py-0.0.*.whl | tail -n 1`
