###Run the predict function of transdecoder to get the open reading frames from each of your transcriptomes###
ls -1 | cut -f1 -d. | while read -r LINE; do TransDecoder.LongOrfs -t $LINE.fasta > ../01-TRANSDECODER_NC/td$LINE.out 2>../01-TRANSDECODER_NC/td$LINE.err; done

####This will create directories  “Species_A.fasta.transdecoder_dir” which contain:
 
###base_freqs.dat
###base_freqs.dat.ok
###longest_orfs.cds
###longest_orfs.gff3
###longest_orfs.pep

####Use BLAST to identify translated sequences that have similarity to the swissprot database. 
#####The resulting BLAST reports can be used by TransDecoder.predict to more accurately predict ORFS 

###make sure you are in the directory with all your *transdecoder_dir subdirectorires, also you will need a swissprot diamond database for this###
ls *transdecoder_dir | grep dir | cut -f1 -d. | while read -r LINE; do diamond blastp -p 18 -e 1e-5 -d /usr/local/uniprot/swissprot \
-q $LINE.fasta.transdecoder_dir/longest_orfs.pep > $LINE.diamond.out 2> $LINE.diamond.err; done

####use your blast results to predict the likely coding regions ####
ls *transdecoder_dir | grep dir | cut -f1 -d. | while read -r LINE; do TransDecoder.Predict -t $LINE.fasta \
--retain_blastp_hits $LINE.diamond.out --cpu 18 > $LINE.td.p.out 2> $LINE.td.p.err; done

#####This will generate several files in your $SPECIES.transdecoder_dir/ directories 
#####and the following files in in the folder that you ran the command
 
####$SPECIES.fasta.transdecoder.bed
####$SPECIES.fasta.transdecoder.gff3
####$SPECIES.fasta.transdecoder.cds (used later in selection tests)
####$SPECIES.fasta.transdecoder.pep (used in next step)

###Now use cd-hit for clustering and further reducing transcript/gene redundancy###
for file in *transdecoder.pep; do cd-hit -i $file -o $file"_cdhit90.pep" -c 0.9; done &
###rename the sequences and files for your orthofinder run####
mkdir ../new_names/
mkdir ../new_names/cdhit_cds
mkdir ../new_names/cdhit_peps
###first peps###
for file in *.transdecoder.pep_cdhit90.pep; do f2=${file%%.transdecoder.pep_cdhit90.pep}; sed 's/::/$/' $file | sed 's/Gene.*\$//g' | sed 's/lcl|/lcl_/g' > ../new_names/cdhit_peps/$f2; done
###then cds###
for file in *.transdecoder.cds; do sed 's/::/$/' $file | sed 's/Gene.*\$//g' | sed 's/lcl|/lcl_/g' > ../new_names/cdhit_cds/$file; done
###now use the peps for your Orthofinder directory###
cd ../new_names/cdhit_peps
cp *fasta ../../../Orthofinder/01-AA_cdhit/

