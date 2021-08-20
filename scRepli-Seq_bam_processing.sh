### Script for processing bam files from 10X Genomics Cell Ranger DNA (v1.1.0) for scRepli-Seq
# author: Qian Du
# date: July 2021

## Splitting bam files from 10X Genomics Cell Ranger DNA (v1.1.0) into individual cells/barcodes
# example of HCT116 G1 sample

# extract barcodes
awk -F "\"*,\"*" '{print $1}' ./outputs/hg19/HC1-G1/per_cell_summary_metrics.csv | tail -n +2 | sed 's/-1//' > ./hct_G1_barcodes.txt

# split bam, sort, index and flagstat. This was done in a cluster environment
FILE=./GE/HCT116-G1/outs/possorted_bam.bam

IDXS=$( sed -n ${SGE_TASK_ID}p < ./hct_G1_barcodes.txt )

samtools view -b -h -u -d CB:${IDXS}-1 -o ./sc_bams_hct_G1/hct_G1_${IDXS}.bam $FILE
samtools sort -T ./sc_bams_hct_G1/tmp/hct_G1_${IDXS}_Aligned.out.sorted -o ./sc_bams_hct_G1/hct_G1_${IDXS}.asd.bam ./sc_bams_hct_G1/hct_G1_${IDXS}.bam
samtools index ./sc_bams_hct_G1/hct_G1_${IDXS}.asd.bam
samtools flagstat ./sc_bams_hct_G1/hct_G1_${IDXS}.asd.bam > ./sc_bams_hct_G1/hct_G1_${IDXS}.asd.bam.stats

## Filter bams - remove duplicates, filter for MAPQ > 10 and remove reads mapping to blacklist regions
# example of HCT116 G1 sample, "cell_1"

# hg19.DACblacklist.bed blacklist file
# removing 'chr' from bed to match bam header
sed 's/chr//' hg19.DACblacklist_orig.bed > hg19.DACblacklist.bed

# filter bam, sort, index and flagstat
samtools view -b -h -u -f 2 hct_G1_cell1.asd.bam | samtools view -@ 4 -b -h -F 1024 -q 10 - | bedtools intersect -abam stdin -b ./hg19.DACblacklist.bed -v > ./sc_bams_hct_G1_filt/hct_G1_cell1.bam
samtools sort -T ./sc_bams_hct_G1_filt/tmp/hct_G1_cell1_Aligned.out.sorted -o ./sc_bams_hct_G1_filt/hct_G1_cell1.asd.bam ./sc_bams_hct_G1_filt/hct_G1_cell1.bam
samtools index ./sc_bams_hct_G1_filt/hct_G1_cell1.asd.bam
samtools flagstat ./sc_bams_hct_G1_filt/hct_G1_cell1.asd.bam > ./sc_bams_hct_G1_filt/hct_G1_cell1.asd.bam.stats

## merging control G1 cells
cd ./sc_bams_hct_G1_filt/
mkdir bams_to_pool

# softlinks of bams to pool
cd ./sc_bams_hct_G1_filt/bams_to_pool/
while read IDXS
do
    echo "$IDXS"
    samplename=$( ls .. | grep "$IDXS" | grep ".bam$" )
    echo "$samplename"
    ln -s ../"$samplename" ./"$samplename"
done < idx_merge_hctG1.txt

# merge bams
samtools merge ./sc_bams_hct_G1_filt/pooled_hct_G1.asd.bam ./sc_bams_hct_G1_filt/bams_to_pool/*.bam
samtools index ./sc_bams_hct_G1_filt/pooled_hct_G1.asd.bam
