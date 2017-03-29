#!/bin/bash

# Setup the environment:
# Browse to src/wrapper directory:
export MODELS_HOME=/home/jovyan/models/models-fork/models
cd $MODELS_HOME/src/wrapper

# Run automake:
libtoolize --force && aclocal && autoheader && automake --force-missing --add-missing

# Run autoconf and generate configure script:
autoconf
./configure CPPFLAGS="-I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL \
-I/opt/conda/pkgs/python-2.7.12-0/include/python2.7" --with-adolc=/usr/local/lib64
make clean && make

# Set up wrapper-specific variables:
export PYTHONPATH=$PYTHONPATH:$MODELS_HOME/src/wrapper
# export WRAP_MOD="sampledir.phys_prop"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib64:$MODELS_HOME/src/wrapper/.libs

# Build the properties library and attempt to run the solver:
cd $MODELS_HOME/src/model_example
make clean && make

# Activate virtual env of the solver:
source activate python2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/pkgs/python-2.7.12-0/include/python2.7
python pyomo_minimal.py