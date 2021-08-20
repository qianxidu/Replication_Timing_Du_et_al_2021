### Script for mapping Nanopore data with NGMLR and calling SVs with sniffles
# author: Qian Du
# date: July 2021

# DNA was sequenced on one PromethION flow cell
# reads were base called with Guppy (v3.3.0)

# reads were mapped to hg19 with NGMLR
ngmlr -t 8 -r ./genomes/hg19/genome.fa -q ./HCT116.fastq.gz -o ./HCT116.sam -x ont --bam-fix

samtools view -b ./HCT116.sam | samtools sort -o ./HCT116.sorted.bam

# calculate read depth to get subsampling ratio for DKO1
samtools view -c HCT116.sorted.bam
samtools view -c DKO1.sorted.bam

# subsample DKO1 reads
samtools view -s 0.63405435 -b DKO1.sorted.bam > DKO1_subsamp.sorted.bam

# checking median depth is equivalent between HCT116 and DKO1 - do below for both HCT116 and DKO1 bam files
samtools depth DKO1_subsamp.sorted.bam > DKO1_subsamp_bam.depth

awk '{print $3}' DKO1_subsamp_bam.depth | sort -n  > DKO1_subsamp_count_sorted.txt

wc -l DKO1_subsamp_count_sorted.txt
  #calculate row number that is middle of order by dividing row length by 2
sed '[middle row number]q;d' samtools_count1.sh

# calling SVs with sniffles
sniffles -m HCT116.sorted.bam -v sniffles_HCT116_hg19.vcf --genotype --cluster
