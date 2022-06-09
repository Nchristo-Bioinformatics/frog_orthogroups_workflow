#!/bin/bash
#SBATCH --job-name=transdecoder
#SBATCH --nodes=1
#SBATCH --time=100:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G

for file in *fasta;
do
TransDecoder.LongOrfs -t $file
TransDecoder.Predict -t $file --single_best_only
done
