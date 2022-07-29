diamond blastp -p 18 -e 1e-5 -d /usr/local/uniprot/swissprot \
  -q Species_A.fasta.transdecoder_dir/longest_orfs.pep \
  > Sp_A.diamond.out 2> Sp_A.diamond.err
