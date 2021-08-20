#---
# title: processing Repli-Seq from bam to bigwig
# authors: Nicola Armstrong, Phuc Loi Luu, Qian Du
#---

### Compute PNDV
library(GenomicRanges)
library(GenomicAlignments)
library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)

#--- 150pb
thres.nread.high <- 20

#--- 50kb
thres.nread.low <- 5 # cutoff of 4 was used for allelic replication timing

samples <- c("HCT116", "DKO1")
reps <- c("A", "B")

#--- output dir
outpath <- ".../pndv/"
inpath <- "..." #location of bam files

#--- big window
wb <- 50000

#--- small window
ws <- 150

#--- spacing for big window
s <- 1000

#--- normalize to 1 million
n <- 1000000

#--- collect datasets in each sample
names(samples) <- samples
bamfiles <- lapply(samples, function(s){
			temp <- grep(s, grep(".bam$", list.files(inpath, full.names=T), value=T), value=T)
	})
bamfiles <- do.call(c, bamfiles)
names(bamfiles) <- basename(bamfiles)

#--- read the BAM file by readGAlignments in GenomicRanges
gr_list <- lapply(bamfiles, function(bam_file) as(readGAlignments(bam_file),"GRanges"))
names(gr_list) <- names(bamfiles)

#--- get all the chrom names of human
chromNames <- names(Hsapiens)[1:23]
chromLengths <- seqlengths(Hsapiens)[chromNames]

#--- Creates a compact GRanges representation of bins across specified chromosomes of a human genome
blocks <- genomeBlocks(Hsapiens, chromNames, width=wb, spacing=s)
blocks <- blocks[width(blocks)==wb]
# For allelic Repli-Seq processing, the 50kb sliding windows/blocks were modified so that they only covered regions within each haplotype block and did not extend between haplotype blocks.

#--- Genomic regions were excluded from further analysis if they either contained greater than x reads in a 150 bp window
#--- NB definition of "bad" changes depending on the dataset (total number of reads)
#--- See Du et al. (2019) Nat. Commun. for further information

badblocks <- genomeBlocks(Hsapiens, chromNames, width=ws, spacing=ws)
badblocks <- badblocks[width(badblocks)==ws]

#--- Counts reads inside blocks
badcounts <- annotationBlocksCounts(bamfiles, badblocks)
summary(badcounts)
save(badcounts, file=paste0(outpath, "badcounts.RData"))

#--- get index of bad region (e.g. contained greater than x reads in a 150 bp window)
badind <- c()
for (i in 1:length(bamfiles)){ badind <- c(badind, which(badcounts[,i] > thres.nread.high)) }
badind <- unique(badind)

#--- get coordinates of bad regions
removebad <- badblocks[badind,]

#--- keep reads if there is no overlap with bad regions
bamfiles <- names(bamfiles)
names(bamfiles) <- bamfiles

gr_good <- lapply(bamfiles, function(y){
	t <- endoapply(gr_list[y], function(x) x[!x %over% removebad])
	return(t[[y]])
})
save(gr_good, file=paste0(outpath, "gr_good.RData"))

#--- count number of reads in each blocks after filtering the bad regions
gr_counts <- lapply(gr_good, function(x){
	count <- blocks
	score(count) <- annotationBlocksCounts(x, blocks)
	return(count)
})

#--- total number of alignments for each fraction after removing reads overlapping 150bp bad bins
gr_total <- lapply(gr_good, function(x){
	length(x)-length(x[which(seqnames(x)=="chrM")])-length(x[which(seqnames(x)=="chrY")])
})

#--- normalise to 1 million
gr_norm <- lapply(bamfiles, function(x){
		(score(gr_counts[[x]])/gr_total[[x]])*n
})
save(gr_norm, file=paste0(outpath, "gr_norm.RData"))

#--- form PNDV dataframe for each sample of 6 fractions
pndvs <- list()
inds <- list()
for(s in samples){
	for(r in reps){
		p <- paste0(r, "_", s)
		print(p)
		rep_name <- paste0(r, "_", s)
		rt <- bamfiles[grep(p, bamfiles)]
		names(rt) <- rt
		count <- do.call(cbind, lapply(rt, function(x) gr_norm[[x]]) )
		colnames(count) <- names(rt)
		rs <- rowSums(count)
		pndvs[[rep_name]] <- count/rs*100
		for(i in rt){ pndvs[[rep_name]][which(rs==0), i] <- 0  }
		# Genomic regions were excluded from further analysis if all 6 fracs contained less than x reads in a 50 kb window
		inds[[rep_name]] <- apply(count, 1, function(x){all(x <= thres.nread.low)})
	}
}

#--- merging the bad regions in both replicates in both samples
inds_dt <- do.call(cbind, inds)
inds_merged <- apply(inds_dt, 1, function(x) any(x==T))

#--- remove the bad regions in both replicates in both samples
pndvs_good <- lapply(pndvs, function(x){ x[!inds_merged,] })
blocks <- blocks[!inds_merged,]
seqlengths(blocks) <- seqlengths(Hsapiens)[1:23]
blocks <- resize(blocks, width=1000, fix="center")
save(pndvs, pndvs_good, blocks, file=paste0(outpath, "PNDVs.RData"))

for(rep_name in names(pndvs_good)){
	count <- pndvs_good[[rep_name]]
	for(i in colnames(count)){
		T <- blocks
		score(T) <- count[,i]
		export(T, paste0(outpath, gsub("_filt.bam", "", i), ".bw"), format="bw")
	}
}

#--- WA = 0.917*G1 + 0.75*S1 + 0.583*S2 + 0.417*S3 + 0.25*S4
WA <- list()
for(rep_name in names(pndvs_good)){
  pndv_good <- pndvs_good[[rep_name]]
  G1 <- grep("_G1", colnames(pndv_good), value=T)
  S1 <- grep("_S1", colnames(pndv_good), value=T)
  S2 <- grep("_S2", colnames(pndv_good), value=T)
  S3 <- grep("_S3", colnames(pndv_good), value=T)
  S4 <- grep("_S4", colnames(pndv_good), value=T)
  WA[[rep_name]] <- 0.917*pndv_good[,G1] + 0.75*pndv_good[,S1] + 0.583*pndv_good[,S2] + 0.417*pndv_good[,S3] + 0.25*pndv_good[,S4]
}
names(WA)

#--- make bigWig file
for(x in names(WA)){
  wa <- blocks
  score(wa) <- WA[[x]]
  export(wa, paste0(outpath, x, ".bw"), format="bw")
}
