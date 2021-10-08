#!/bin/bash
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Trims paired end reads in a directory. Make sure read names end in _R1.fastq.gz and _R2.fastq.gz"
   echo
   echo "Syntax: ./trimming_from_reads.sh [-d|n|q|a|m]"
   echo "options:"
   echo "h     Print this Help."
   echo "d     Directory of untrimmed paired end reads."
   echo "n     Name of directory for trimmed read output"
   echo "q     Minimum quality score for trimmomatic (Default:20)"
   echo "a     Adaptor file to use (Default:NexteraPE-PE.fa)"
   echo "m     Minimum read length after trimming (Default:60)"
   echo
}
####set variables###
untrim_reads_dir="untrim_reads_dir"
paired_trim_reads_dir="trim_reads_dir"
quality="20"
adaptorfile="NexteraPE-PE.fa"
minlen="60"
############################################################
# Process the input options. Add options as needed.        #
############################################################
############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
while getopts ":hd:n:q:a:m:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      d) # Enter a name
         untrim_reads_dir=$OPTARG;;
      n) # Enter a name
         paired_trim_reads_dir=$OPTARG;;
      q) # Enter a number
         quality=${OPTARG};;
      a) # Enter fasta file
         adaptorfile=$OPTARG;;
      m) # Enter seq length
         minlen=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

                                                                 ####TRIMMING######
#export untrim_reads_dir="$1"
echo "$untrim_reads_dir is the directory of untrimmed paired reads"
echo "$paired_trim_reads_dir is where the trimmed reads will be"
echo "$quality is the minimum quality for trimming reads"
echo "$adaptorfile will be the adaptor file used"
echo "Minumum sequence length after trimming is $minlen"
for f1 in $untrim_reads_dir/*R1.fastq.gz;
do
f2=${f1%%R1.fastq.gz}"R2.fastq.gz"
f3=${f1%%R1.fastq.gz}"R1_p.fastq.gz"
f4=${f1%%R1.fastq.gz}"R1_u.fastq.gz"
f5=${f1%%R1.fastq.gz}"R2_p.fastq.gz"
f6=${f1%%R1.fastq.gz}"R2_u.fastq.gz"
#echo $f1
#echo $f2
#echo $f3
#echo $f4
#echo $f5
#echo $f6
java -jar /share/apps/trinityrnaseq-Trinity-v2.6.6/trinity-plugins/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 $f1 $f2 $f3 $f4 $f5 $f6 \
ILLUMINACLIP:/share/apps/trinityrnaseq-Trinity-v2.6.6/trinity-plugins/Trimmomatic-0.36/adapters/$adaptorfile:2:$quality:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$minlen 
done

mkdir $paired_trim_reads_dir
mv $untrim_reads_dir/*R*_p.fastq.gz $paired_trim_reads_dir
rm $untrim_reads_dir/*R*_u.fastq.gz
