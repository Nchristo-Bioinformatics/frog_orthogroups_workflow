###Getting duplicated orthogroups###

grep -e "N18" -e "N23" -e "N27" -e "N30" -e "N33" -e "N36" -e "N38" Duplications_per_Species_Tree_Node.tsv
# N36     111     64
# N27     169     29
# N23    143     7
# N18     57      4
# N38     47      47
# N30     121     13
# N33     245     145

awk '$2 == "N18"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n4 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N23"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n7 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N27"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n29 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N30"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n13 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N33"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n145 | sort -u -k1,1  >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N36"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n64 | sort -u -k1,1  >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N38"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n47 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt

####Getting orthogroups with at least two pristimantis sequences, for filtering with phylopypruner####
####getting the trees for phylopypruner###
awk '{nzeros=0; for(col=36; col<=NF; col++) {if($col == 0) {nzeros++}} {if(nzeros < 6) {print $0}}}' Orthogroups.GeneCount.tsv | cut -f1 | while read -r LINE; do cp ../OrthoFinder/Results_May17/Resolved_Gene_Trees/$LINE* OG_trees/; done
###getting the orthofinder alignnments for phylopypruner####
awk '{nzeros=0; for(col=36; col<=NF; col++) {if($col == 0) {nzeros++}} {if(nzeros < 6) {print $0}}}' Orthogroups.GeneCount.tsv | cut -f1 | while read -r LINE; do cp ../OrthoFinder/Results_May17/MultipleSequenceAlignments/$LINE* OG_pep_alignments/; done
###getting the original fastas for aligning with pal2nal####
awk '{nzeros=0; for(col=36; col<=NF; col++) {if($col == 0) {nzeros++}} {if(nzeros < 6) {print $0}}}' Orthogroups.GeneCount.tsv | cut -f1 | while read -r LINE; do cp ../OrthoFinder/Results_May17/Orthogroup_Sequences/$LINE* OG_fastas/; done

###Put all the alignments and trees into a directory, run phylopypruner on the directory as so####
phylopypruner --overwrite --dir filtering_seqs_orthofinder_alignments/ --output pruning_output_2/ --threads 16 --include Oreobates_cruralis mzua2401_Trinity mzua2448_Trinity mzua2499_Trinity mzua2511_Trinity mzua2538_Trinity mzua2591_Trinity mzua2610_Trinity --subclades subclades.txt --trim-divergent .80 --trim-lb 5 --min-taxa 20

####remove all phylopypruner alignments without mzua in them####
grep -r -L -Z 'mzua' . | xargs --null rm
sed -i 's/|/_/g' *fa
sed -i 's/_lcl_/_lcl|/g' *fa
