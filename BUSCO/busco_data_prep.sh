####pristimantis sequences would not be taken by BUSCO for some reason, perhaps clustering with cdhit and transdecoder made them not compatible?###

###so we renamea them###
for file in Pristimantis*cds; do awk -v awkvar="$file" '/^>/{print ">awkvar" ++i; next}{print}' < $file | sed "s/awkvar/$file/g" > $file.busco_names.fasta; done
