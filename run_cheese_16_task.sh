#!/bin/bash
################################
## SLURM batchjob script for
## Elmer on LUMI
##
## copyleft 2023-06-21
##    CSC-IT Center for Sciencce
##
################################
 
#SBATCH --time=00:15:00
#SBATCH --job-name=ChEESE_test
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=standard-g

####### change to your project #######
#SBATCH --account=project_462000007

####### change to numbers of nodes and MPI tasks ###
####### NB: we provide meshes for 128,256,512 and 1024 partitions #####
#######     do the math by matching the product of next entries   #####
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --mem=400G
#SBATCH --exclusive

##SBATCH --gpus=8
##SBATCH --gpus-per-task=1

################## OpenMP Stuff ##########
## use only if you undersubscribe
## the MPI tasks
##########################################
#SBATCH --cpus-per-task=7
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "running OpenMP on $SLURM_CPUS_PER_TASK"
#export KMP_AFFINITY=compact
#export KMP_DETERMINISTIC_REDUCTION=yes

#export CRAY_ACC_NO_ASYNC=1

## These control USM behaviour
export CRAY_ACC_USE_UNIFIED_MEM=0
export HSA_XNACK=0
#export CRAY_ACC_DEBUG=2
#export LIBOMPTARGET_KERNEL_TRACE=2
#export LIBOMPTARGET_INFO=0
export MPICH_GPU_SUPPORT_ENABLED=1
export FI_CXI_RX_MATCH_MODE=hybrid

. ./8mpi-7omp-8gpu.sh

#ml myelmer/cheese
ml LUMI/24.03 partition/G PrgEnv-gnu/8.5.0 rocm  buildtools cray-hdf5/1.12.2.11 cray-netcdf/4.9.0.11  cray-libsci/24.03.0 
ml myelmer/cheese

# this loads the spack-PrgEnv-gnu cray-libsci (BLAS, LAPACK) version
###### best to use this for audits! #########
#ml use /appl/lumi/spack/23.03/0.19.2/share/spack/modules/linux-sles15-zen2
#ml use /appl/local/csc/soft/eng/elmer/spack/23.03/0.19.2/modules/tcl/linux-sles15-zen2/
###### this loads the container version of Elmer
# module load elmer/latest
###### this loads the PregEnv/gnu version ###
###### (using cray-libsci (BLAS, LAPACK)
#module load elmer/gcc-cray
###### make it so! ######### 


export CASE_SOLVER_MODULE=poisson.so

export ELMER_MESH_NAME=./winkel_audit
echo %%%% MESH_NAME=$ELMER_MESH_NAME
#srun rocprof --basenames on --timestamp on -i metrics.txt -o ${ELMER_MESH_NAME}.csv ElmerSolver case.sif
#srun rocprof --basenames on --timestamp on -i metrics.txt --hip-trace ElmerSolver case.sif
#srun bash -c '. select_gpu && env > env_${SLURM_LOCALID}.txt'
srun bash -c '. select_gpu; ElmerSolver case.sif'
