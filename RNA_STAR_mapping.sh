#!/bin/sh

# Grid Engine options
#$ -N map
#$ -cwd
#$ -l h_vmem=8G
#$ -m baes
#$ -pe sharedmem 4

# If you plan to load any software modules, then you must first initialise the modules framework.
. /etc/profile.d/modules.sh
# Choose the staging environment
export OMP_NUM_THREADS=$NSLOTS

# Then, you must load the modules themselves
module load igmm/apps/STAR/2.7.8a

#STAR Index
ref='LALO.fa'
anno='LALO.gtf'
#basic options to generate genome indices are as follows:
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./Mask_ref --genomeSAindexNbases 10 --genomeFastaFiles ${ref} --sjdbGTFfile ${anno} --sjdbOverhang 149

#STAR mapping
for i in `less Sample_list.txt`;do \
prefix=$1 #Sample_list.txt $i

listR1=./raw_data/${prefix}/${prefix}_1.fq.gz
listR2=./raw_data/${prefix}/${prefix}_2.fq.gz
output=./mapping/${prefix}
ref=./Mask_ref
ls  ${listR1} ${listR2}
echo ${output}
#Prepare the STAR index first: ref

STAR --genomeDir ${ref} --runThreadN 4 --readFilesCommand zcat --readFilesIn ${listR1} ${listR2} --outFilterType BySJout --outSAMunmapped None --outReadsUnmapped Fastx --outFileNamePrefix ${output} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 7900000000

done
