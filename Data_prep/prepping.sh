conda create --prefix /home/$$$$/orthofinder_workflow
conda activate /home/$$$$/orthofinder_workflow
conda install -y -c bioconda perl-db-file   
conda install -y -c bioconda orthofinder
conda install -y -c conda-forge r-ape
conda install -y -c bioconda hyphy
conda install -c bioconda perl-uri
cpan URI::Escape   

PAML
There's an issue with PAML from conda. If you want to run PAML on your own machine download a binary from here: http://abacus.gene.ucl.ac.uk/software/
It is already installed and in your PATH on nzinga. 


 
conda activate /home/$$$$/orthofinder_workflow

###put the species name in front of the sequences, will help for later###
perl -pi.orig -e 's/^>/>SPECIES_NAME|/' SPECIES.fasta
for file in *.fasta; do
  f1=${file%%.fasta}
  echo "$f1"
  perl -pi.orig -e "s/^>/>${f1}|/" "$file"
done
