#!/bin/bash

source ./env_elmer.sh

### change to the path where you want to install hipfort
export ROOT_HIPFORT=$PWD/install
###

rm -rf *.o vecadd
#compile kernel file (.ccp)
hipcc -fno-gpu-rdc -fPIC --offload-arch=gfx90a -c hipfort/test/f2008/vecadd/hip_implementation.cpp

#compile fortran file (.f90)
gfortran-13 -cpp -I$ROOT_HIPFORT/include/hipfort/amdgcn -c hipfort/test/f2008/vecadd/main.f03 

#link object files 
gfortran-13 hip_implementation.o main.o -L$ROOT_HIPFORT/lib -lhipfort-amdgcn -L$ROCM_PATH/lib -lamdhip64 -Wl,-rpath=$ROCM_PATH/lib  -lstdc++ -o vecadd
