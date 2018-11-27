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



