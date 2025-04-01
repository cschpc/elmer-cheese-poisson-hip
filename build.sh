ROOT_HIPFORT=./install

hipcc -O3 -fno-gpu-rdc -fPIC --offload-arch=gfx90a -c kernelPoisson.cpp
gfortran-13 -cpp -I$ROOT_HIPFORT/include/hipfort/amdgcn -c fortran_hip_interface.F90

elmerf90 Poisson.F90 fortran_hip_interface.o  kernelPoisson.o -L$ROOT_HIPFORT/lib -lhipfort-amdgcn  -I$ROOT_HIPFORT/include/hipfort/amdgcn/ -L$ROCM_PATH/lib -lamdhip64 -Wl,-rpath=$ROCM_PATH/lib  -I./ -lstdc++  -o poisson.so
# elmerf90 -L$ROOT_HIPFORT/lib -lhipfort-amdgcn  -I$ROOT_HIPFORT/include/hipfort/amdgcn/ -L$ROCM_PATH/lib -lamdhip64 -Wl,-rpath=$ROCM_PATH/lib  -I./ -lstdc++ Poisson.F90 kernelPoisson.o -o Poisson-x.so
