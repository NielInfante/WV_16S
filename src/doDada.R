library(dada2); packageVersion("dada2")
library(tidyverse)


# Set working directory
baseDir <- "~/"
setwd(baseDir)

# List of where the files are

fnFs <- sort(list.files('~/', pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files('~/', pattern="_R2.fastq.gz", full.names = TRUE))


# Read the meta
meta <- read.table('meta', header=T, sep="\t", stringsAsFactors = F)
sample.names <- meta$Sample

# Take a look at some quality
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Create filtered files
filt_path <- paste0(baseDir,"/filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,230), 
										 trimLeft = 10,     
                     maxN=0,  truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,
#                     maxEE=c(3,4)
                     maxEE=2
                     )
head(out)

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# Learn Error Rates   - Might consider pooling for small sample sets (pool=T)
errF <- learnErrors(filtFs, multithread=T)
png('pics/error_F.png')
	plotErrors(errF, nominalQ = T)
dev.off()

errR <- learnErrors(filtRs, multithread=T)
png('pics/error_R.png')
	plotErrors(errR, nominalQ = T)
dev.off()

# Sample Inference
# Maybe run again with pool=T when I've got a lot of extra time
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=F)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=F)
dadaFs[[1]]
dadaRs[[1]]

# Merge Paired Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct Sequence Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

saveRDS(seqtab, 'Data/seqtab_before_sizing.rds')

# Remove sequences that are inappropriatly sized
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(401,438)]
dim(seqtab2)
#table(nchar(getSequences(seqtab2)))
seqtab <- seqtab2

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
100 * sum(seqtab.nochim)/sum(seqtab)

# 95% of the reads are still there, even though 96% of identified sequences were removed as chimeas

# Overview of where things were filtered - see if you need to go back and tweek some parameters
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track)
track$Sample <- sample.names
#track$Survived <- round(100 * track$nonchim / track$input)
head(track)

tidyTrack <- track %>% gather(Step, Count, input:nonchim)
tidyTrack$Step <- factor(tidyTrack$Step, levels = c('input', 'filtered','denoised','merged','tabled','nonchim')) 
head(tidyTrack)

png('pics/filtering.png', width=1400, height=1000)
	ggplot(tidyTrack) + geom_line(aes(x=Step, y=Count, group=Sample, color=Sample))
dev.off()

seqtab <- seqtab.nochim
saveRDS(seqtab, file='Data/seqtab.rds')

# Assign Taxonomy   - If you want species, need to use RDP or Silva, not green genes
taxa <- assignTaxonomy(seqtab.nochim, "~/genome/gg_13.8/gg_13_8_train_set_97.fa.gz", multithread = T)

taxa.print <- taxa
row.names(taxa.print) <- NULL
head(taxa.print)

#write.table(seqtab.nochim, file="test", sep="\t", col.names = NA, quote=F)

# Make a phylogenetic tree Using phangorn
library(phangorn)
library(DECIPHER)
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# Here we first construct a neighbor-joining tree, and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree using the neighbor-joining tree as a starting point.

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)



# Some Extra
library(phyloseq)
#packageVersion("phyloseq")
library(ggplot2)

meta <- read.table('meta', header=T, sep="\t", stringsAsFactors = F)
rownames(meta) <- meta$Sample

# Issue making tree, try it without for the moment
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), sample_data(meta), tax_table(taxa),phy_tree(fitGTR$tree))

#ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), sample_data(meta), tax_table(taxa))



saveRDS(ps, file="data/PhyloseqObject.rds")
save(ps, file="data/PS.Rdata")

plot_richness(ps, x="Treatment", measures=c("Shannon", "Simpson"), color="Cage") + theme_bw()


ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
plot_ordination(ps, ord.nmds.bray, color="Cage", title="Bray NMDS")


# Barplot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="SampleID", fill="Family") + facet_wrap(~Treatment, scales="free_x")+geom_bar(stat="identity")
plot_bar(ps.top20, x="SampleID", fill="Phylum") + facet_wrap(~Condition, scales="free_x")+geom_bar(stat="identity")
