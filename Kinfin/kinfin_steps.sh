     # Install kinfin:
####you will need python 2.7, and some packages, probably easiest way to do this is as follows:###
conda create -y -n py27 python=2.7
conda activate py27

pip install SciPy
pip install ete3
pip install Matplotlib
pip install Docopt

git clone https://github.com/DRL/kinfin.git
cd kinfin
./install
./kinfin

      # This will convert the concatenated InterProScan files into a format readable by KinFin.
./kinfin/scripts/iprs2table.py -i all_proteins.tsv --domain_sources Pfam

      # Copy Orthofinder files that are needed to the KinFin directory:
cp ~/polar_workflow/Orthofinder/01-AA_cdhit/OrthoFinder/Results_Jul29/WorkingDirectory/SequenceIDs.txt .
cp ~/polar_workflow/Orthofinder/01-AA_cdhit/OrthoFinder/Results_Jul29/WorkingDirectory/SpeciesIDs.txt .
cp ~/polar_workflow/Orthofinder/01-AA_cdhit/OrthoFinder/Results_Jul29/Orthogroups/Orthogroups.txt .

      # Create the KinFin configuration file:
echo '#IDX,TAXON' > config.txt
sed 's/: /,/g' SpeciesIDs.txt | \
cut -f 1 -d"." >> config.txt

      # Run the KinFin functional annotation script:
      # This failed for me until I commented out 'ax.set_facecolor('white')' on lines 681 and 1754 of ./kinfin/src/kinfin.py. Once I commented them out with vim it worked fine.
./kinfin/kinfin --cluster_file Orthogroups.txt --config_file config.txt --sequence_ids_file SequenceIDs.txt --functional_annotation functional_annotation.txt
