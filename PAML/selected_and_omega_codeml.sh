###get background and foreground omegas and proportions in different files
for file in *codeml; do grep -A1 "proportion" $file  | cut -f1,11,12,13,14 -d" " | paste -s -d, - | sed 's/proportion  //g' | sed 's/background /background,/g'| sed 's/ /,/g' | sed 's/,,/,/g' >> background_props_omegas.csv; done
for file in *codeml; do grep -A2 "proportion" $file  | sed '/background/d'| cut -f1,11,12,13,14 -d" " | paste -s -d, - | sed 's/proportion  //g' | sed 's/foreground /foreground,/g'| sed 's/ /,/g' | sed 's/,,/,/g' >> foreground_props_omegas.csv; done
cat *props_omegas.csv > fore_back_props_omegas.csv

####make this the header for fore_back_props_omegas.csv#### ---> prop_2a,prop_2b,foreground_background,omega_2a,omega_2b


###for t-tests on proportions of selected sites, omegas for each codeml hypothesis###
for file in *codeml; do echo $file; grep -e proportion -e "background w" $file | sed 's/background w/background.w/g' | \
sed 's/proportion /proportion/g' | sed 's/  / /g' | cut -f1,6,7,8 -d" " | cut -f2,3 -d" " | paste -d " " - - | while read -r pine; 
do echo "$file $pine background" | sed 's/.nuc.nostop.trim.fa.alt.codeml//g' | sed 's/ /,/g'  >> pristi_fore_back_props_omegas.csv; done; done

for file in *codeml; do echo $file; grep -e proportion -e "foreground w" $file | sed 's/foreground w/foreground.w/g' | \
sed 's/proportion /proportion/g' | sed 's/  / /g' | cut -f1,6,7,8 -d" " | cut -f2,3 -d" " | paste -d " " - - | while read -r pine; 
do echo "$file $pine foreground" | sed 's/.nuc.nostop.trim.fa.alt.codeml//g' | sed 's/ /,/g'  >> pristi_fore_back_props_omegas.csv; done; done

##add headers to your files##
