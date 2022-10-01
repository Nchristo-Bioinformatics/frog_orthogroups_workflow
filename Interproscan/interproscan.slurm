#!/bin/bash
#SBATCH --job-name=interpro
#SBATCH --nodes=1
#SBATCH --time=200:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=250G


cd ~/polar_workflow/Interproscan/interproscan-5.52-86.0
for proteinName in *fasta; do ./interproscan.sh -i $proteinName -d out/ -t p --goterms -appl Pfam -f TSV; done



