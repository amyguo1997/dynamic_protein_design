#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q gpu.q
#$ -pe smp 1
#$ -N colabfold
#$ -l mem_free=32G
#$ -l h_rt=1:00:00
#$ -l compute_cap=61,gpu_mem=5000M
#$ -r y
#$ -j y
#$ -o job_outputs
#$ -e job_outputs
#$ -t 1:10000
#$ -tc 500

conda activate colabfold

export CUDA_VISIBLE_DEVICES=$SGE_GPU
module load Sali cuda

pdb_id=$1

alias python=python3
python run_colabfold.py $pdb_id

