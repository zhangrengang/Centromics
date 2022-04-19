#!/bin/bash

# install dependency
(cd REPcluster && python3 setup.py install)

# install
python3 setup.py install

