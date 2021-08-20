### Script for allelically mapping Repli-Seq data using phased SNPs
# author: Qian Du
# date: July 2021

# Using SNPs called from Nanopore data, phased with medaka
# Extracting shared heterozygous SNPs between HCT116 and DKO1 for filtering in R

bgzip hct_phased.vcf
bcftools index hct_phased.vcf.gz

bcftools isec -Oz -p intersection hct_phased.vcf.gz dko_phased.vcf.gz
#intersection/0000.vcf.gz	for records private to	hct_phased.vcf.gz
#intersection/0001.vcf.gz	for records private to	dko_phased.vcf.gz
#intersection/0002.vcf.gz	for records from hct_phased.vcf.gz shared by both	hct_phased.vcf.gz dko_phased.vcf.gz
#intersection/0003.vcf.gz	for records from dko_phased.vcf.gz shared by both	hct_phased.vcf.gz dko_phased.vcf.gz

cd intersection
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%TYPE[\t%SAMPLE=%GT][\t%SAMPLE=%GQ][\t%SAMPLE=%PS]\t%INFO/pos1\t%INFO/q1\t%INFO/pos2\t%INFO/q2\n' 0002.vcf.gz > 0002_extracted.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%TYPE[\t%SAMPLE=%GT][\t%SAMPLE=%GQ][\t%SAMPLE=%PS]\t%INFO/pos1\t%INFO/q1\t%INFO/pos2\t%INFO/q2\n' 0003.vcf.gz > 0003_extracted.vcf

# creating reference genomes for each haplotype
bcftools view -R hct_phased_filtered.bed -Oz 0002.vcf.gz -o hct_hap_filt.vcf.gz
bcftools view -R dko_phased_filtered.bed -Oz 0003.vcf.gz -o dko_hap_filt.vcf.gz

bcftools index hct_hap_filt.vcf.gz
bcftools index dko_hap_filt.vcf.gz

bcftools consensus -H 1pIu -f ./genomes/hg19/genome.fa ./intersection/hct_hap_filt.vcf.gz > ./phased_fa/hg19_hct_hap1.fa
bcftools consensus -H 2pIu -f ./genomes/hg19/genome.fa ./intersection/hct_hap_filt.vcf.gz > ./phased_fa/hg19_hct_hap2.fa
bcftools consensus -H 1pIu -f ./genomes/hg19/genome.fa ./intersection/dko_hap_filt.vcf.gz > ./phased_fa/hg19_dko_hap1.fa
bcftools consensus -H 2pIu -f ./genomes/hg19/genome.fa ./intersection/dko_hap_filt.vcf.gz > ./phased_fa/hg19_dko_hap2.fa

samtools faidx hg19_hct_hap1.fa
samtools faidx hg19_hct_hap2.fa
samtools faidx hg19_dko_hap1.fa
samtools faidx hg19_dko_hap2.fa

bowtie-build ./phased_fa/hg19_hct_hap1.fa ./bowtieindex/hg19_hct_hap1
bowtie-build ./phased_fa/hg19_hct_hap2.fa ./bowtieindex/hg19_hct_hap2
bowtie-build ./phased_fa/hg19_dko_hap1.fa ./bowtieindex/hg19_dko_hap1
bowtie-build ./phased_fa/hg19_dko_hap2.fa ./bowtieindex/hg19_dko_hap2

# HCT116 and DKO1 Repli-Seq was mapped to allelic haplotypes, e.g.
zcat HCT116_G1_A.fastq.gz | bowtie -v 0 -m 1 --tryhard --best --strata --time --trim5 6 --phred33-quals --threads 8 --sam ./bowtieindex/hg19_hct_hap2 - HCT116_G1_A_hap2.aln.sam
samtools view -Sbt ./bowtieindex/hg19_hct_hap2.fai HCT116_G1_A_hap2.aln.sam > HCT116_G1_A_hap2.aln.bam
samtools sort HCT116_G1_A_hap2.aln.bam HCT116_G1_A_hap2.ash.bam
picard.sam.CleanSam INPUT=HCT116_G1_A_hap2.ash.bam OUTPUT=HCT116_G1_A_hap2.cleaned.bam TMP_DIR=[./tmp] VALIDATION_STRINGENCY=LENIENT
picard.sam.MarkDuplicates INPUT=HCT116_G1_A_hap2.cleaned.bam OUTPUT=HCT116_G1_A_hap2.asd.bam METRICS_FILE=HCT116_G1_A_hap2.asd.bam.dupl ASSUME_SORTED=true TMP_DIR=[./tmp] VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=9 CREATE_MD5_FILE=true

# extracting only reads that map to haplotype SNPs
samtools view -b -L phased_filtered_snps.bed HCT116_G1_A_hap2.asd.bam > ./bowtie_filt/HCT116_G1_A_hap2.asd.bam

# next step is Du_et_al_Repli-Seq_processing_bam_to_PNDV.R
