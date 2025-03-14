#!/bin/bash -l
#SBATCH --job-name=hipfort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --partition=standard-g
#SBATCH --hint=nomultithread
#SBATCH --account=project_462000031
#SBATCH --exclusive
#SBATCH --gres=gpu:1

source ./env_elmer.sh

srun ./vecadd
