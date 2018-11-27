library(phyloseq)
library(ggplot2)
library("gridExtra")
library(dplyr)

setwd('~/projects/Anderson')

# Play with shiny, if you like
#install.packages("shiny")
#shiny::runGitHub("shiny-phyloseq","joey711")

# Load Data
ps <- readRDS('data/PhyloseqObject.rds')






# Refactor data - I want to put in more metadata  - Old, metadata in phyloseq should be good
#tree <- phy_tree(ps)
#tax <- tax_table(ps)
#otu <- otu_table(ps)
#sam <- sample_data(ps)

#m2 <- read.table("~/projects/Cuff/Rashel/Both/meta", header=T, sep="\t", stringsAsFactors = F)
#rownames(m2) <- m2$SampleID
#m2$Experiment <- substr(m2$SampleID,1,1)
#m2$Age_cat <- cut(m2$Age, breaks=c(0,18,22,26,30))
#m2$Group <- substr(m2$Individual, 1, 2)

#ps <- phyloseq(otu, sample_data(m2), tax,tree)

#####  Filtering

rank_names(ps)

# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Get rid of NA's and uncharacterized phyla - not always the best thing to do.
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
table(tax_table(ps0)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

# Are any phyla minimally represented?
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# or

myPrev <- prevdf %>% group_by(Phylum) %>% summarize(n=n(),sum=sum(Prevalence),mean=mean(Prevalence),max=max(Prevalence))
myPrev %>% print(n=100)

# If so, filter them out
# Define phyla to filter
filterPhyla <- myPrev %>% filter(max <= 3) %>% select(Phylum)

#filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla$Phylum)

# Don't do any filtering
#ps1 <- ps0


## Prevalence Filtering

# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)


# Just check what this does
# Compute prevalence of each feature, store as data.frame
prevdf2 = apply(X = otu_table(ps2),
							 MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),
							 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf2 = data.frame(Prevalence = prevdf2,
										TotalAbundance = taxa_sums(ps2),
										tax_table(ps2))



# skipping prevelance filtering
ps2 <- ps1



# Agglomerate taxa

# I don't think I will do this at this time, but I'm including it for later reference

# How many genera would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))  # combine all features that descend from the same genus
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)

# Or can do it using tree height, if you don't trust the taxonomy
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)

# See what the agglomerated trees look like:
multiPlotTitleTextSize = 8
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))

# group plots together
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)


###
#ps2 <- ps1

###  Abundance value transformation
##     Normalize by read counts, etc

plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("p__Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Genotype",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# Transform to relative abundance. Save as new object.
psra = transform_sample_counts(ps, function(x){x / sum(x)})

plotBefore = plot_abundance(ps2,"")
plotAfter = plot_abundance(ps2ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2, plotBefore, plotAfter)

# Subset by Taxonomy
psOrd = subset_taxa(ps2ra, Order == "o__Erysipelotrichales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)


### Moving On


ps <-ps2

saveRDS(ps, file='data/Phyloseq_filtered.rds')


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

# Flatten the low end
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1

ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))

tax <- tax_table(ps)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)

main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))

row_scores <- row_scores %>%
  left_join(sample_data(pslog))
col_scores <- col_scores %>%
  left_join(tax)






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








