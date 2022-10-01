#!/bin/bash
#SBATCH --job-name=kinfin
#SBATCH --nodes=1
#SBATCH --time=200:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G


###need python 2.7 for kinfin###
source activate py27
###get these files from Orthofinder results###
./kinfin/kinfin --cluster_file Orthogroups.txt --config_file config.txt --sequence_ids_file SequenceIDs.txt --functional_annotation functional_annotation.txt
