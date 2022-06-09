

grep -e "N23" -e "N27" -e "N30" -e "N32" -e "N35" -e "N37" -e "N39" Duplications_per_Species_Tree_Node.tsv
#N32     121     13
#N27     143     7
#N23     58      4
#N39     45      45
#N30     169     29
#N35     245     145
#N37     113     55



awk '$2 == "N23"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n4 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N27"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n7 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N30"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n29 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N32"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n13 | sort -u -k1,1  >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N35"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n145 | sort -u -k1,1  >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N37"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n55 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
awk '$2 == "N39"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n45 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt

#############How many duplicated orthogroups?############################
cut -f1 ALL_dupgroups_nodes_support.txt | sort | uniq | wc -l
#236

####how many have an annotation?########
cut -f1 ALL_dupgroups_nodes_support.txt | sort | uniq | while read -r LINE; do grep $LINE cluster_domain_annotation.* >> dupgroup_kinfin_annotations.txt; done
cut -f1 dupgroup_kinfin_annotations.txt | cut -f2 -d: | sort | uniq | wc -l
##228###
