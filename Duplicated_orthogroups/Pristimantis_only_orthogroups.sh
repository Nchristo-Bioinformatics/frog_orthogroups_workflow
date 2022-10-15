awk '{ for (i=2; i<=21; ++i) if ($i != 0) next; print }' Orthogroups.GeneCount.tsv > halfway_pristi_only_orthos.tsv

awk '{ for (i=29; i<=42; ++i) if ($i != 0) next; print }' halfway_pristi_only_orthos.tsv > only_pristi_orthos.tsv
sed 's/\t/,/g' only_pristi_orthos.tsv > only_pristi_orthos.csv
