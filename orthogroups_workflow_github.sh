####Run first orthofinder to make orthogroups###
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







####remove all phylopypruner alignments without mzua in them####
grep -r -L -Z 'mzua' . | xargs --null rm
sed -i 's/|/_/g' *fa
sed -i 's/_lcl_/_lcl|/g' *fa


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



#########Get protein sequences for your selected proteins from the original, unaligned protein fasta orthogroup files from Orthofinder###
for file in *.fa; do f2=${file%%.fa}"_pruned.csv"; f3=${file%%.fa}"_pruned.pep.fa" ; seqtk subseq $file ../pruned_csvs/$f2 > $f3; done
####align the protein sequences####
for file in *pep.fa; do f2=${file%%pep.fa}"pep.aln"; mafft --anysymbol --auto --quiet --thread 8 $file > ../pruned_pep_alignments/$f2; done 

####Now make nucleotide alignment files from your pruned nucleotide fastas, and your aligned, pruned protein fastas you just made####
for file in *_pruned.aln; do f2=${file%%_pruned.aln}"_pruned.nuc.fa"; f3=${file%%_pruned.aln}".nuc.aln"; pal2nal.pl $file ../pruned_nucs/$f2 -output fasta > ../nuc_alignments_from_pal2nal/$f3; done
###trim your nucleotide alignments###
for file in *aln; do f2=${file%%.aln}"trim.aln"; trimal -in $file -out ../trimmed_nuc_alignments_pal2nal/$f2 -fasta -gappyout; done
####now remove stop codons####
hyphy ~/miniconda3/envs/dammit/lib/hyphy/TemplateBatchFiles/CleanStopCodons.bf Universal OG0001347.alntrim No/Yes OG0001347.alntrim_nostopls
