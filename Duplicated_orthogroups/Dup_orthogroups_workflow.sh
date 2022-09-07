grep -e "N21 -e "N24" -e "N27"-e "N30"-e "N33 -e "N36" -e "N38" Duplications_per_Species_Tree_Node.tsv
#N36     36      16
#N30     36      3
#N21     70      1
#N24     89      4
#N33     90      42
#N27     87      10
#N38     11      11



awk '$2 == "N21"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n1 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N24"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n4 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N27"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n10 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N30"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n3 | sort -u -k1,1  >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N33"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n42 | sort -u -k1,1  >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N36"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n16 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N38"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n11 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt


#############How many duplicated orthogroups?############################
cut -f1 ALL_dupgroups_nodes_support.txt | sort | uniq | wc -l

#####Let's gather the protein fastas for duplicated orthogroups###
mkdir 01-SEQS
cut -f1 ALL_dupgroups_nodes_support.txt | sort | uniq | while read -r LINE; do cp ../Orthofinder/01-AA_cdhit/OrthoFinder/Results_Jul29/Orthogroup_Sequences/$LINE".fa" 01-SEQS/; done

####how many have an annotation?########
cut -f1 ALL_dupgroups_nodes_support.txt | sort | uniq | while read -r LINE; do grep $LINE all_annotations.tsv >> Pristi_oreo_dupgroup_kinfin_annotations.txt; done


###PRUNING PARALOGUES TO MAKE ORTHOGROUPS SMALLER###
mkdir DUPLICATED_PHYLOPYPRUNER
##cut -f1 ALL_dupgroups_nodes_support.txt | sort | uniq | while read -r LINE; do cp ../Orthofinder/01-AA_cdhit/OrthoFinder/Results_Jul29/MultipleSequenceAlignments/$LINE".fa" Duplicated_Phylotreepruner/; done
#cut -f1 ALL_dupgroups_nodes_support.txt | sort | uniq | while read -r LINE; do cp ../Orthofinder/01-AA_cdhit/OrthoFinder/Results_Jul29/Gene_Trees/$LINE* Duplicated_Phylotreepruner/; done
#cd DUPLICATED_PHYLOPYPRUNER
#for file in *_tree.txt; do f2=${file%%_tree.txt}".tree"; mv $file $f2; done
#cd ..
###LETS SEE IF THIS WORKS##########
#phylopypruner --dir DUPLICATED_PHYLOPYPRUNER --overwrite --threads 8 --include Oreobates_cruralis Pristimantis_andinognomus Pristimantis_orestes Pristimantis_sp1 Pristimantis_sp2 Pristimantis_sp3 Pristimantis_sp4 Pristimantis_sp5 --mask longest --min-support 0.5 --min-taxa 20 --prune MI
####OKAY IT DOESNT WORK, ONLY ONE ORTHOGROUP ALIGNMENT REMAINS########
mkdir 02-INDEX_FASTAS
cd 01-SEQS/
###GET TOP 2 ISOFORMS AND MAKE CSV###
for file in *fa; 
do samtools faidx $file; done
for file in *.fa.fai; do f2=${file%%.fa.fai}".csv"; cut -f1 -d"|" $file | sort -u |  parallel -k "grep {} -m 1 <(sort -k2,2 -r --sort g $file ) | cut -f1" > $f2; done
cp *csv ../01-SEQS
###SUBSET PROTEIN ORTHOFINDER FILES####
cd ..
mkdir 03-PRUNED_SEQS
cd 01-SEQS
for file in *.csv; do f1=${file%%.csv}".fa" ;f2=${file%%.csv}"_pruned.fa"; samtools faidx $f1 -r $file > ../03-PRUNED_SEQS/$f2; done
####Now lets get the CDS NUCLEOTIDE SEQUENCES###
cd ..
perl ../perl_scripts/get_corresponding_cds.pl ../Transdecoder/new_names/cdhit_cds/ 03-PRUNED_SEQS/ 04-CDS
