#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --nodes=1
#SBATCH --time=200:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G

#conda create -n buscoing -c conda-forge -c bioconda busco=5.4.2
conda activate buscoing
ln -s ~/polar_workflow/Transdecoder/new_names/cdhit_cds/*cds .
for file in *.transdecoder.cds; do f2=${file%%.transdecoder.cds}"_busco"; busco -i $file -l metazoa_odb10 -o $f2 -m trans -c 8; done
