#!/bin/bash

current=$PWD
### change to the path where you want to install hipfort
instdir=$current/install
###
source ./env_elmer.sh

rm -rf $instdir
mkdir -p $instdir

export FC=$(which gfortran-13)
export CC=$(which gcc-13)
export CXX=$(which g++-13)
export GCCAR=$(which gcc-ar-13)
export GCCRANLIB=$(which gcc-ranlib-13)

rm -rf hipfort
git clone https://github.com/ROCmSoftwarePlatform/hipfort
cd hipfort
mkdir build ; cd build
cmake -DHIPFORT_INSTALL_DIR=$instdir -DHIPFORT_COMPILER=$FC \
      -DHIPFORT_AR=$GCCAR -DHIPFORT_RANLIB=$GCCRANLIB ..
make -j install
cd $current
