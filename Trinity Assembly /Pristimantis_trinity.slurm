#!/bin/bash
#SBATCH --job-name=T-2628
#SBATCH --nodes=1
#SBATCH --time=125:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G

module load trinity
module load samtools/1.9
orig_directory=`pwd`
job_directory="trinity$$"
mkdir -p /scratch/$job_directory ##move your normalized read pairs to the job directory##
cd /scratch/$job_directory
ln -s $orig_directory/2628frog_S43_R1_p.fastq.gz $orig_directory/2628frog_S43_R2_p.fastq.gz .
echo "job directory is $job_directory"
/share/apps/trinityrnaseq-Trinity-v2.6.6/Trinity --seqType fq --SS_lib_type RF --left 2628frog_S43_R1_p.fastq.gz --right 2628frog_S43_R2_p.fastq.gz --full_cleanup --CPU 10 --inchworm_cpu 10 --max_memory 100G --verbose --output /scratch/$job_directory
# put the results in a tar file
cd ..
# copy results back to your home directory
cp *.fasta $orig_directory
