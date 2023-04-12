###The purpose of this script is to get the messy output from Kinfin, aka a species having or not having a certain Pfam or IPR domain,
###and put it into a count matrix for regression analysis. This seems more appropriate for fragmented sequences rather than just lumping them all into 
###an orthogroup
import csv, re
import sys
import pandas as pd

#inny = open(sys.argv[1])
inny = open('cluster_domain_annotation.IPR.pristimantis_groups.txt')
#outy = open(sys.argv[2], 'w')

reader = csv.DictReader(inny, delimiter='\t')
domains = {}
orthogroups = []
for r in reader:
        Orthogroup = r['#cluster_id']
        desc = r['domain_description']
        combined_ortho_domainss = '%s_%s' % (Orthogroup, desc)
        combined_ortho_domains = combined_ortho_domainss.replace(',', '')
        fracts = r['TAXA_with_domain']
        sp_fracts = fracts.split(',')
        #print(sp_fracts)
        sp_counts = {}
        for i in sp_fracts:
                #print(i)
                #fract_total = sp_fracts.split('/')[2]
                sp_fracts_str = str(i)
                fract_nototal = sp_fracts_str.split('/')[0]
                species = fract_nototal.split(':')[0]
                count_nototal = fract_nototal.split(':')[1]
                sp_counts[species] = count_nototal
                #count_list.append(sp_counts)
        domains[combined_ortho_domains] = sp_counts
transposed_df = pd.DataFrame(domains).transpose()
transposed_df_2 = pd.DataFrame(transposed_df).fillna(0)
transposed_df_3 = transposed_df_2.reindex(sorted(transposed_df_2.columns), axis=1)
pd.DataFrame(transposed_df_3).to_csv('pristi_domain_counts_final.csv')
