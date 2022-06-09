#######Let us install interproscan######
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.52-86.0/interproscan-5.52-86.0-64-bit.tar.gz
tar -xzvf interproscan-5.52-86.0-64-bit.tar.gz
cd ./interproscan-5.52-86.0
mkdir ./out

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

      # Install kinfin:
git clone https://github.com/DRL/kinfin.git
cd kinfin
./install
./kinfin

      # This will convert the concatenated InterProScan files into a format readable by KinFin.
./kinfin/scripts/iprs2table.py -i all_proteins.tsv --domain_sources Pfam

      # Copy Orthofinder files that are needed to the KinFin directory:
cp ~/OrthoFinder/NO_CDHIT_NEW_TREE_CHECK/OrthoFinder/Results_Jun01/WorkingDirectory/SequenceIDs.txt .
cp ~/OrthoFinder/NO_CDHIT_NEW_TREE_CHECK/OrthoFinder/Results_Jun01/WorkingDirectory/SpeciesIDs.txt .
cp ~/OrthoFinder/NO_CDHIT_NEW_TREE_CHECK/OrthoFinder/Results_Jun01/Orthogroups/Orthogroups.txt .

      # Create the KinFin configuration file:
echo '#IDX,TAXON' > config.txt
sed 's/: /,/g' SpeciesIDs.txt | \
cut -f 1 -d"." >> config.txt

      # Run the KinFin functional annotation script:
      # This failed for me until I commented out 'ax.set_facecolor('white')' on lines 681 and 1754 of ./kinfin/src/kinfin.py. Once I commented them out with vim it worked fine.
./kinfin/kinfin --cluster_file Orthogroups.txt --config_file config.txt --sequence_ids_file SequenceIDs.txt --functional_annotation functional_annotation.txt
