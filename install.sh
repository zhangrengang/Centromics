#!/bin/bash

# install dependency
[ ! -s REPcluster/setup.py ] && \
	echo 'Submodule REPcluster do not exit. Please use `git clone --recurse-submodules`' && \
	exit
(cd REPcluster && python3 setup.py install)

# install
python3 setup.py install

