mkdir Interproscan
cd interproscan

#######Let us install interproscan######
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.52-86.0/interproscan-5.52-86.0-64-bit.tar.gz
tar -xzvf interproscan-5.52-86.0-64-bit.tar.gz
cd ./interproscan-5.52-86.0
mkdir ./out
####need peptide files for this####
 cp ../../Orthofinder/01-AA_cdhit/*fasta .
####script for annotating all your protein fasta files####
#!/bin/bash
#SBATCH --job-name=interpro
#SBATCH --nodes=1
#SBATCH --time=200:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G


for proteinName in *fasta; do ./interproscan.sh -i $proteinName -d out/ -t p --goterms -appl Pfam -f TSV; done

# Now concatenate all of the annotated protein files into a single file which can be passed to KinFin for orthogroup-level functional annotation.
cat ./interproscan-5.52-86.0/out/*.tsv > all_proteins.tsv
