 2002  cd polar_workflow/Duplicated_orthogroup_comparisons/
 2003  ls
 2004  head Duplications.tsv
 2005  awk '$2 == "Scuttiger_sikimmensis"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n |  sort -u -k1,1
 2006* grep Scut Duplicationsawk '$2 == "N21"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n | tail -n1 | sort -u -k1,1 >> ALL_dupgroups_nodes_support.txt
 2007  grep Scut Duplications.tsv
 2008  awk '$2 == "Scutiger_sikimmensis"' Duplications.tsv | cut -f1,2,3,4 | sort -k4 -n
 2009  awk '$2 == "Scutiger_sikimmensis"' Duplications.tsv | cut -f1,2,3,4 | wc -l
 2010  awk '$2 == "Scutiger_sikimmensis"' Duplications.tsv | cut -f1,2,3,4 > S_sikimmensis_dups.txt
 2011  awk '$2 == "Nanorana_parkeri"' Duplications.tsv | cut -f1,2,3,4 > N_parkeri_dups.txt
 2012  ls
 2013  head ALL_dupgroups_nodes_support.txt
 2014  head N_parkeri_dups.txt
 2015  cut -f1 N_parkeri_dups.txt | sort | uniq | wc -l
 2016  join <(sort S_sikimmensis_dups.txt) <(sort N_parkeri_dups.txt)
 2017  join <(sort S_sikimmensis_dups.txt) <(sort N_parkeri_dups.txt) | cut -f1 | sort | uniq | wc -l
 2018  join <(sort S_sikimmensis_dups.txt) <(sort N_parkeri_dups.txt) | tr ' ' '\t' > common_Nparkeri_Sscutiger_dup_orthos.txt
 2019  cat common_Nparkeri_Sscutiger_dup_orthos.txt
 2020  cut -f1 common_Nparkeri_Sscutiger_dup_orthos.txt  | sort | uniq | wc -l
 2021  cut -f1 common_Nparkeri_Sscutiger_dup_orthos.txt  | sort | uniq
 2022  join <(sort common_Nparkeri_Sscutiger_dup_orthos.txt) <(sort ALL_dupgroups_nodes_support.txt) | tr ' ' '\t'
 2023  join <(sort common_Nparkeri_Sscutiger_dup_orthos.txt) <(sort ALL_dupgroups_nodes_support.txt) | tr ' ' '\t' > possible_high_elevation_dups.txt
 2024  cut -f1 possible_high_elevation_dups.txt | sort | uniq | wc -l
 2025  cat possible_high_elevation_dups.txt
 2026  sort -u k1,1 possible_high_elevation_dups.txt
 2027  sort -u -k1,1 possible_high_elevation_dups.txt
