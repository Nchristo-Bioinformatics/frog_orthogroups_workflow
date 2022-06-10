####Run first orthofinder on the directory with all your translated proteins fastas from Transdecoder to make orthogroups###
./orthofinder -t 20 -a 16 -f NO_CDHIT_ORTHOGROUPS_CHECKING -T fasttree -M msa

####Getting orthogroups with at least two pristimantis sequences, for filtering with phylopypruner####
####getting the trees for phylopypruner###
awk '{nzeros=0; for(col=36; col<=NF; col++) {if($col == 0) {nzeros++}} {if(nzeros < 6) {print $0}}}' Orthogroups.GeneCount.tsv | cut -f1 | while read -r LINE; do cp ../OrthoFinder/Results_May17/Resolved_Gene_Trees/$LINE* OG_trees/; done
###getting the orthofinder alignnments for phylopypruner####
awk '{nzeros=0; for(col=36; col<=NF; col++) {if($col == 0) {nzeros++}} {if(nzeros < 6) {print $0}}}' Orthogroups.GeneCount.tsv | cut -f1 | while read -r LINE; do cp ../OrthoFinder/Results_May17/MultipleSequenceAlignments/$LINE* OG_pep_alignments/; done
###getting the original fastas for aligning with pal2nal####
awk '{nzeros=0; for(col=36; col<=NF; col++) {if($col == 0) {nzeros++}} {if(nzeros < 6) {print $0}}}' Orthogroups.GeneCount.tsv | cut -f1 | while read -r LINE; do cp ../OrthoFinder/Results_May17/Orthogroup_Sequences/$LINE* OG_fastas/; done

###Put all the alignments and trees into a directory, run phylopypruner on the directory as so####
phylopypruner --overwrite --dir filtering_seqs_orthofinder_alignments/ --output pruning_output_2/ --threads 16 --include Oreobates_cruralis mzua2401_Trinity mzua2448_Trinity mzua2499_Trinity mzua2511_Trinity mzua2538_Trinity mzua2591_Trinity mzua2610_Trinity --subclades subclades.txt --trim-divergent .80 --trim-lb 5 --min-taxa 20


###make a new species tree from concatenated phylopypruner alignment, will be using this to re-run orthofinder###

iqtree -s supermatrix.fas -pre supermatrix_mzuaGTR20 -nt 10 -m GTR20 -bb 5000 -o Bombina_orientalis,Bombina_bombina,Bombina_variegata

###iqtree gives you an unrooted tree, and Orthofinder needs it rooted, so let's root it at our outgroup in R###
R
library(ape)
tree <- read.tree('supermatrix_mzuaGTR20.contree')
outgroop <- c("Bombina_orientalis", "Bombina_variegata", "Bombina_bombina")
treetry <- root(tree2, resolve.root = TRUE, outgroup = outgroop)
treetry

#Phylogenetic tree with 41 tips and 40 internal nodes.

#Tip labels:
	#Xenopus_laevis_lcl, Xenopus_tropicalis_lcl, Odorrana_tormota, Odorrana_margaretae, Rana_clamitans_Rana_clamitans, Rana_catesbeiana, ...
#Node labels:
	#Root, 100, 100, 100, , 100, ...

#Rooted; includes branch lengths.
write.tree(treetry, 'supermatrix_mzuaGTR20.contree_rooted')


####re-run orthofinder with your new species tree, to make sure that duplications are correctly put on the tree####
./orthofinder -t 20 -a 16 -f NO_CDHIT_NEW_TREE_CHECK -s supermatrix_mzuaGTR20.contree_rooted -T fasttree -M msa
