library(tidyverse)
library(phyloseq)
library(vegan)

setwd('')
ps <- readRDS('Data/PhyloseqObject_filtered.rds')
ps <- transform_sample_counts(ps, function(otu) {otu/sum(otu)})

setseed(3)


# Subset the data, then call
psn <- subset_samples(ps, Run=='2' & Gonad=='intact')
multiBeta(psn, "Run2 Intact XXF vs XYF", "vegan/Run2_Intact_XXF_XYF", "Chromosome", "Chromosome", "NMDS")



# Call this function, it will do all the beta measures you add.
multiBeta <- function(ps, title, outfile_base, label="", variable, ordMethod){

	pv_bray     <- beta_with_plot(ps, title, paste0(outfile_base, '_bray.png'), variable, label, ordMethod, 'bray', "Bray Curtis")
	pv_jaccard  <- beta_with_plot(ps, title, paste0(outfile_base, '_jaccard.png'), variable, label, ordMethod, 'jaccard', "Jaccard")
	pv_unifrac  <- beta_with_plot(ps, title, paste0(outfile_base, '_unifrac.png'), variable, label, ordMethod, 'unifrac', "Unifrac")
	pv_wunifrac <- beta_with_plot(ps, title, paste0(outfile_base, '_wunifrac.png'), variable, label, ordMethod, 'wunifrac', "Weighted Unifrac")

	# If you want to put the results in a md table
	print(paste(pv_bray, pv_jaccard, pv_unifrac, pv_wunifrac, sep=" | "))
	
}

beta_with_plot <- function(ps, title, outfile, label="", variable, ordMethod, distance='bray', distanceName="Bray Curtis"){
	
	# Calculate p value
	distance_matrix <- phyloseq::distance(ps, method=distance)
	sample_df <- data.frame(sample_data(ps))

	formula <- as.formula(paste0('distance_matrix ~ ', variable, collapse = ''))
	ad <- adonis(formula, data=sample_df)
	pval <- ad$aov.tab$`Pr(>F)`[1]

	# Only consider the 50 most numerous taxa
	ps = prune_taxa(names(sort(taxa_sums(ps), TRUE)[1:50]), ps)

	formula <- as.formula(paste0(' ~ ', variable, collapse = ''))
	ord = ordinate(ps, formula = formula, ordMethod, distance)
	ordplot <- plot_ordination(ps, ord, "samples", color=variable)

	p <- ordplot + 
		stat_ellipse(type = "norm", linetype = 2) +
		stat_ellipse(type = "t", geom='polygon', alpha=0.1, aes_string(fill=variable), show.legend = F) +
		labs(title = title, subtitle = paste('P-value:',pval), caption=distanceName, color=label) +
	theme_bw()
 
	ggsave(outfile, plot=p, device='png', width=3, height=3)
	
	return(pval)
}



