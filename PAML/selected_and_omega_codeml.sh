###get background and foreground omegas and proportions in different files
for file in *codeml; do grep -A1 "proportion" $file  | cut -f1,11,12,13,14 -d" " | paste -s -d, - | sed 's/proportion  //g' | sed 's/background /background,/g'| sed 's/ /,/g' | sed 's/,,/,/g' >> background_props_omegas.csv; done
for file in *codeml; do grep -A2 "proportion" $file  | sed '/background/d'| cut -f1,11,12,13,14 -d" " | paste -s -d, - | sed 's/proportion  //g' | sed 's/foreground /foreground,/g'| sed 's/ /,/g' | sed 's/,,/,/g' >> foreground_props_omegas.csv; done
cat *props_omegas.csv > fore_back_props_omegas.csv

####make this the header for fore_back_props_omegas.csv#### ---> prop_2a,prop_2b,foreground_background,omega_2a,omega_2b
