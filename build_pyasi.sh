#!/bin/sh

PYASI_VERSION=`python3 -c "from build import toml_loads; print(toml_loads(open('pyasi/pyproject.toml').read())['project']['version'])"`

python3 -m build pyasi
