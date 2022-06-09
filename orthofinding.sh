#!/bin/bash
#SBATCH --job-name=Orthofrogs
#SBATCH --nodes=1
#SBATCH --time=200:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --cpus-per-task=20
#SBATCH --mem=350G

##MAY 2022, REDOING TO CHECK####
#./orthofinder -t 20 -a 16 -f NO_CDHIT_NEW_TREE_CHECK -s iqtree_checking_pub.contree_rooted_mzuatopology -T fasttree -M msa

###using supermatrix from phylopypruner###
./orthofinder -t 20 -a 16 -f NO_CDHIT_NEW_TREE_CHECK -s supermatrix_mzuaGTR20.contree_rooted -T fasttree -M msa
