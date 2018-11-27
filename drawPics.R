# Drawing several plots

library(phyloseq)
library(tidyverse)

ps <- readRDS('Data/PhyloseqObject.rds')


plot_richness(ps, x="Treatment", measures=c("Shannon", "Simpson"), color="Cage") + theme_bw()


ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
plot_ordination(ps, ord.nmds.bray, color="Cage", title="Bray NMDS")


# Barplot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="SampleID", fill="Family") + facet_wrap(~Treatment, scales="free_x")+geom_bar(stat="identity")
plot_bar(ps.top20, x="SampleID", fill="Phylum") + facet_wrap(~Condition, scales="free_x")+geom_bar(stat="identity")






#qplot(sample_data(ps)$Experiment, geom = "histogram") + xlab("Experiment")
#qplot(log10(rowSums(otu_table(ps)))) + xlab("Logged counts-per-sample")
qplot((rowSums(otu_table(ps)))) + xlab("Logged counts-per-sample")

pslog <- transform_sample_counts(ps, function(x) log(1 + x))
#sample_data(pslog)$age_binned <- cut(sample_data(pslog)$Time), breaks = c(0, 100, 200, 400))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "bray")



evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "Treatment", shape='Cage') +
	labs(col = c('Sex',"Genotype")) +
	coord_fixed(sqrt(evals[2] / evals[1])) + geom_point(size=3)

evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "Treatment") +
	labs(col = "Treatment") +
	coord_fixed(sqrt(evals[2] / evals[1])) + geom_point(size=3)


# PCoA

out.bc.log <- ordinate(pslog, method = "MDS", distance = "bray")

evals <- out.bc.log$values$Eigenvalues
plot_ordination(pslog, out.bc.log, color = "Treatment", shape="Cage") +
	coord_fixed(sqrt(evals[2] / evals[1])) +
	labs(col = "Treatment", shape = "Cage") + geom_point(size=3)


#DPCoA

out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$values$Eigenvalues
plot_ordination(pslog, out.dpcoa.log, color = "Treatment", shape="Cage") +
	coord_fixed(sqrt(evals[2] / evals[1])) +
	labs(col = "Genotype", shape = "Tab")

#? Not working


out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "Treatment", shape="Cage") +
	coord_fixed(sqrt(evals[2] / evals[1])) +
	labs(col = "Treatment", shape = "Cage") + geom_point(size=3)


# PCA on Ranks

plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
	coord_fixed(sqrt(evals[2] / evals[1]))

evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "Experiment",
								shape = "Time") +
	coord_fixed(sqrt(evals[2] / evals[1])) +
	labs(col = "Experiment", shape = "Time")

plot_ordination(pslog, out.wuf.log, type = "species", color = "Phylum") +
	coord_fixed(sqrt(evals[2] / evals[1]))

abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))




# Some quickies for my presentation


psn <- subset_samples(ps,  Genotype=='WT' & (SurguryType != '.') & (SurguryType != 'N/A') & 
												(Age_Day == 0 | Age_Day == 3) )




pslog <- transform_sample_counts(ps, function(x) log(1 + x))

out.bc.log <- ordinate(pslog, method = "MDS", distance = "bray")

evals <- out.bc.log$values$Eigenvalues
plot_ordination(pslog, out.bc.log, color = "Genotype", shape="Sex") +
	coord_fixed(sqrt(evals[2] / evals[1])) +
	labs(col = "Genotype", shape = "Sex")

png('pics/all_PCA.png')
dev.off()



png('pics/WT_Alpha.png')
#plot_richness(psn, x="SurguryType", measures=c("Shannon", "Simpson"), color="Age_Day") + theme_bw()
plot_richness(psn, x="SurguryType", measures=c("Shannon", "Simpson")) + theme_bw() + labs(title="WT at Zero or Three Days")
dev.off()

plot_richness(ps, x="interaction(SurguryType,Age_Day)", measures=c("Shannon", "Simpson"), color="Genotype") + 
	theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))



psra = transform_sample_counts(psn, function(x){x / sum(x)})

png('pics/T3_bar.png', width=1200, height=1000)
plot_bar(psra, fill="Phylum") + geom_bar(stat="identity")
dev.off()








