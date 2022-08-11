#!/bin/sh

rm dist/*.whl
python3 -m build
pip uninstall -y asi4py
pip install -U  dist/asi4py-0.0.*-py3-none-any.whl

