###From genecount file, get only orthogroups where values for everything EXCEPT pristimantis species is 0####

awk '{ for (i=2; i<=21; ++i) if ($i != 0) next; print }' Orthogroups.GeneCount.tsv > halfway_pristi_only_orthos.tsv

awk '{ for (i=29; i<=42; ++i) if ($i != 0) next; print }' halfway_pristi_only_orthos.tsv > only_pristi_orthos.tsv
sed 's/\t/,/g' only_pristi_orthos.tsv > only_pristi_orthos.csv


###Now let's check what these are###
cut -f1 only_pristi_orthos.csv -d, | while read -r LINE; do grep $LINE annotation_file_bilu_enrichment.csv; done
OG0017548,PF00042|Globin,GO:0020037|heme binding ,IPR000971|Domain Globin
OG0029708,PF14214|Helitron helicase-like domain at N-terminus,NA,IPR025476|Domain Helitron helicase-like domain
OG0029709,PF01200|Ribosomal protein S28e,GO:0003735|structural constituent of ribosome ;GO:0005840|ribosome ;GO:0006412|translation ,IPR000289|Family Ribosomal protein S28e
OG0029719,PF01157|Ribosomal protein L21e,GO:0003735|structural constituent of ribosome ;GO:0005840|ribosome ;GO:0006412|translation ,IPR001147|Family Ribosomal protein L21e
OG0042432,PF00042|Globin,GO:0020037|heme binding ,IPR000971|Domain Globin
OG0042435,PF07742|BTG family,NA,IPR002087|Domain Anti-proliferative protein
OG0042457,PF13602|Zinc-binding dehydrogenase,NA,NA
OG0042487,PF07710|P53 tetramerisation motif,GO:0051262|protein tetramerization ,"IPR010991|Domain p53, tetramerisation 
