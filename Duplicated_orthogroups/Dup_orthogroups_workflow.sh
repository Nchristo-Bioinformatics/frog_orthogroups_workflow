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

mkdir OG_alignments
cut -f1 dupgroup_kinfin_annotations.txt | cut -f2 -d: | sort | uniq| while read -r LINE; do cp ~/OrthoFinder/NO_CDHIT_NEW_TREE_CHECK/OrthoFinder/Results_Jun01/MultipleSequenceAlignments/$LINE* OG_alignments/; done
cut -f1 dupgroup_kinfin_annotations.txt | cut -f2 -d: | sort | uniq| while read -r LINE; do cp ~/OrthoFinder/NO_CDHIT_NEW_TREE_CHECK/OrthoFinder/Results_Jun01/Resolved_Gene_Trees/$LINE* OG_alignments/; done
