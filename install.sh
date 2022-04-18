#!/bin/bash
# install 
python setup.py install

# install dependency
(cd REPcluster && python setup.py install)
