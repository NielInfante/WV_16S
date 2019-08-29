library(phyloseq)
library(ggplot2)
library("gridExtra")
library(dplyr)

setwd('~/')

# Play with shiny, if you like
#install.packages("shiny")
#shiny::runGitHub("shiny-phyloseq","joey711")

# Load Data
ps <- readRDS('Data/PhyloseqObject.rds')


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

myPrev <- prevdf %>% group_by(Phylum) %>% summarize(n=n(),sum=sum(Prevalence),
																										mean=mean(Prevalence),
																										max=max(Prevalence),
																										totAbund=sum(TotalAbundance))
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

saveRDS(ps, file='Data/Phyloseq_filtered.rds')






