#!/bin/bash
# Grid Engine options
#$ -N feature
#$ -cwd
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -l h_rt=16:00:00
#$ -m baes


date
featurecounts=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/Install/subread-2.0.0-Linux-x86_64/bin/featureCounts


wkdir=./mapping/
anno=./Mask_LALO.fasta_rm_MT.gtf
#anno=./Mask_LALO.fasta_rm_MT_mRNA_CDS_rmMT.gff3

$featurecounts -p -P -C -t exon -B -D 5000 -g gene_id -T 8 -a ${anno} -o Output_samples.txt ${wkdir}/*Aligned.sortedByCoord.out.bam  >& Output_samples.featurecounts_sortedByCoord.log

#END#
