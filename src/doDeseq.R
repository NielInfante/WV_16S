# DESeq in Phyloseq

library(phyloseq)
library(DESeq2)
library(RColorBrewer)
library(tidyverse)
library(gplots)
library(ggrepel)

setwd('~/projects/Anderson/Oct_2018')
ps <- readRDS('data/Phyloseq_filtered.rds')

outDir <- "~/projects/Anderson/Oct_2018/deseq"
outPrefix <- 'Naive_Vehicle'

#sample_data(ps)$Age_Day <- as.factor(sample_data(ps)$Age_Day)
#names(sample_data(ps))[7] <- 'SurgeryType'

psn <- subset_samples(ps, Treatment=='Naive' | Treatment=='Vehicle')


#head(sample_data(ps)$Experiment)

psn <- prune_samples(sample_sums(psn) > 100, psn)

intGroup <- 'Treatment'
contrast <- c('Treatment', 'Naive', 'Vehicle')
dds <- phyloseq_to_deseq2(psn, ~  Cage + Treatment)
dds <- DESeq(dds, test="Wald", fitType = "parametric")
#intGroup <- 'Experiment'
#dds <- phyloseq_to_deseq2(psn, ~  Experiment)
#dds <- DESeq(dds, test="Wald", fitType = "parametric")

# If that throws an error (every gene contains at least one zero...) The do this:
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
	exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
dds <- phyloseq_to_deseq2(psn, ~ Cage + Treatment)
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = DESeq(dds, fitType="local")

# Moving on

res <- results(dds, contrast=contrast)
sigtab <- res[which(res$padj < 0.05),]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psn)[rownames(sigtab), ], "matrix"))
sigtab <- sigtab[order(sigtab$padj, decreasing = F),]
head(sigtab)
names(sigtab)


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
	scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
#x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
#x = sort(x, TRUE)
#sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))


# family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))



png(paste0(outDir, '/', outPrefix, '_overview.png'), width=1600, height=1000)
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
	geom_abline(slope=0, intercept = 0) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=15))
dev.off()

#ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
#	theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=15))



name <- paste(outDir, '/', outPrefix, '_significant.txt', sep="") 
write.table(sigtab, file=name, sep="\t", quote=F, row.names=F)

# Now do the standard stuff
#vsd <- vst(dds, blind=F)
vsd <- varianceStabilizingTransformation(dds, blind=F)

##########  Sanity Check
# Plot counts of most significant, to check if fold change is right
png(paste0(outDir, '/',outPrefix,'_sanity.check.png'))
plotCounts(dds, gene=rownames(sigtab)[1], intgroup = intGroup, main=sigtab[1,]$Genus, pch=19)
dev.off()

#########  MA Plot   #########

name <- paste(outDir, '/', outPrefix, '_MAplot.png', sep="") 
png(name)
plotMA(dds, main=outPrefix)
dev.off()

name <- paste(outDir, '/', outPrefix, '_MA_sig.png', sep="") 
png(name)
sigtab[1:20,] %>% ggplot(aes(x=baseMean,y=log2FoldChange)) + geom_point() + 
			geom_text_repel(aes(label=Genus)) + 
			scale_x_log10() + geom_hline(yintercept=0, color='black') +  theme_bw()
dev.off()


#########  Heatmap   #########

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(SampleID))

name <- paste(outDir, '/', outPrefix, '_heatmap.png', sep="") 
png(name)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(6, 6), density.info = 'none')
dev.off()


###########  Cluster   ###########

name <- paste(outDir, '/', outPrefix, '_cluster.png', sep="") 
png(name)
plot(hclust(dist(t(assay(vsd)))), label=with(colData(dds), paste0(SampleID,':',Treatment)), main=outPrefix, xlab='', sub='')
dev.off()

############  PCA    ###########s

name <- paste(outDir, '/', outPrefix, '_PCA.png', sep="") 
png(name)
print(plotPCA(vsd, intgroup=intGroup))
dev.off()

#tiff(file=name, width=1800, height=1200, units='px', res=300)

name <- paste(outDir, '/', outPrefix, '_PCA_names.png', sep="") 
png(name)
p <- plotPCA(vsd, intgroup=intGroup)
#p <- p + geom_text_repel(aes_string(x="PC1", y="PC2", label=colData(dds)$Sex), point.padding = unit(2,"points"))
p <- p + geom_text_repel(aes(x=PC1, y=PC2, label=colData(dds)$SampleID), point.padding = unit(2,"points"))
print(p)
dev.off()

####  Volcano

name <- paste(outDir, '/', outPrefix, '_volcano.png', sep="") 
png(name)

par(pch = 16)
with(res, plot(log2FoldChange, -log10(pvalue), main = outPrefix))
with(subset(res, padj < 0.05), points(log2FoldChange, -log10(pvalue), col = "red"))
with(subset(res, abs(log2FoldChange) > 2), points(log2FoldChange, -log10(pvalue),  col = "orange"))

with(subset(res, padj < 0.05 & abs(log2FoldChange) > 2), points(log2FoldChange,  -log10(pvalue), col = "green"))

# Add legend
legend("topleft", legend = c("FDR<0.05", "|LFC|>2", "both"), pch = 16, col = c("red", "orange", "green"))

# Label Extra significant points
#with(subset(res, padj < 0.05 & abs(log2FoldChange) > 2), textxy(log2FoldChange, -log10(pvalue), labs = Gene, cex = 1))

# Label all significant
#with(subset(res, padj < 0.05), textxy(log2FoldChange, -log10(pvalue), labs = Gene, cex = 1))

dev.off()





