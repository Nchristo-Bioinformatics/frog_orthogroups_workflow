####Go to your phylopypruner results directory####

cd output_alignments
mkdir ../filtered_phylopypruner_output_alignments
cp *pruned.fa ../filtered_phylopypruner_output_alignments

####remove all phylopypruner alignments without mzua in them####
cd ../filtered_phylopypruner_output_alignments
grep -r -L -Z 'mzua' . | xargs --null rm
sed -i 's/|/_/g' *fa
sed -i 's/_lcl_/_lcl|/g' *fa


###let's make some csv files to help get the pruned nucleotide and peptide files for pal2nal###
mkdir ../pruned_csvs/
for file in OG*fa; do f2=${file%%fa}"csv"; grep ">" $file | sed 's/>//g' > ../pruned_csvs/$f2; done

###trim all these pruned alignment files, make new protein trees for Hyphy to use####
for file in *pruned.fa; do f2=${file%%pruned.fa}"trim.aln"; trimal -in $file -out $f2 -fasta -gappyout; done
for file in *trim.aln; do iqtree -s $file -st AA -fast -mset LG,WAG,JTT,DAYHOFF  -lbp 1000 -nt 16; done
###label the foregrounds for these trees, for pristimantis AND direct developers####
for file in *treefile; do f2=${file%%treefile}"dd_tree"; hyphy ~/Cleaning_up_mylife/NO_CDHIT/hyphy-analyses/LabelTrees/label-tree.bf --tree $file  --regexp 'mzua|Oreo' --output ../$f2; done
for file in *treefile; do f2=${file%%treefile}"pristi_tree"; hyphy ~/Cleaning_up_mylife/NO_CDHIT/hyphy-analyses/LabelTrees/label-tree.bf --tree $file  --regexp 'mzua' --output ../$f2; done

###because we need the original, unaltered coding sequences, I find it just easier to concatenate all our *cds files from Transdecoder into one file, aka. bigtranscripts_1.fa###
cd ..
mkdir CDs
cd CDs
mv ~/OrthoFinder/transdecoding/*cds .
for file in *.fasta.transdecoder.cds; do f2=${file%%.fasta.transdecoder.cds}".fasta"; mv $file $f2; done
sed -i "s|\.|_|g" *fasta
for file in *.fasta; do sed -i "s/>/>${file}/g" $file; done
sed -i 's/\.fasta/|/g' *fasta
cat *fasta >> ../pruned_csvs/bigtranscripts_1.fasta

#########Get nucleotide sequences for your selected proteins from bigtranscripts_1.fasta, which you just created###
cd pruned_csvs
mkdir ../pruned_nucs
for file in *.csv; do f2=${file%%.csv}".nuc.fa"; seqtk subseq bigtranscripts_1.fasta $file > ../pruned_nucs/$f2; done
cd ..
#########Get protein sequences for your selected proteins from the original, unaligned protein fasta orthogroup files from Orthofinder###
cd OG_fastas
mkdir ../pruned_peps
for file in *.fa; do f2=${file%%.fa}"_pruned.csv"; f3=${file%%.fa}"_pruned.pep.fa" ; seqtk subseq $file ../pruned_csvs/$f2 > ../pruned_peps/$f3; done
cd ..
####align the protein sequences####
cd pruned_peps
mkdir ../pruned_pep_alignments
for file in *pep.fa; do f2=${file%%pep.fa}"pep.aln"; mafft --anysymbol --auto --quiet --thread 8 $file > ../pruned_pep_alignments/$f2; done 





####Now make nucleotide alignment files from your pruned nucleotide fastas, and your aligned, pruned protein fastas you just made####
for file in *_pruned.aln; do f2=${file%%_pruned.aln}"_pruned.nuc.fa"; f3=${file%%_pruned.aln}".nuc.aln"; pal2nal.pl $file ../pruned_nucs/$f2 -output fasta > ../nuc_alignments_from_pal2nal/$f3; done
###trim your new pal2nal nucleotide alignments###
for file in *aln; do f2=${file%%.aln}"trim.aln"; trimal -in $file -out ../trimmed_nuc_alignments_pal2nal/$f2 -fasta -gappyout; done
####now remove stop codons####
for file in *aln; do f2=${file%%aln}"nostop.fa"; HYPHYMP ~/miniconda3/pkgs/hyphy-2.3.14-h5466e78_0/lib/hyphy/TemplateBatchFiles/CleanStopCodons.bf Universal $file No/Yes ../nostop_alignments/$f2; done
