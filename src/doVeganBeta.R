# From: http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html

library(vegan)
library(tidyverse)
library(ggplot2)
library(phyloseq)

setwd('~/projects/Anderson/Oct_2018')
ps <- readRDS('data/Phyloseq_filtered.rds')

psn <- subset_samples(ps, Treatment=='Naive' | Treatment=='Vehicle')

dat <- psmelt(ps)
dat <- psmelt(psn)


dat$Taxa <- paste(dat$Kingdom,dat$Phylum,dat$Class,dat$Order,dat$Family,dat$Genus,sep = ";")

dat$Kingdom <- NULL
dat$Phylum <- NULL
dat$Class <- NULL
dat$Order <- NULL
dat$Family <- NULL
dat$Genus <- NULL
dat$Species <- NULL



#abund_table <- dat %>% select(OTU, Taxa, Sample, Abundance) %>% spread(Sample, Abundance)
								
abund_table_Base <- dat %>% select(Taxa, Sample, Abundance) %>% 
							group_by(Taxa, Sample) %>% summarize(Ab=sum(Abundance)) %>% 
							spread(Sample, Ab) %>% ungroup()


meta<-read.table("meta", sep="\t", header=T, stringsAsFactors = F)


row.names(meta) <- meta$SampleID


# Filter to test some subsets

meta_table <- meta %>% filter(Treatment=='Naive' | Treatment=='Vehicle')

#meta_table <- meta %>% filter(Age_Day==270 & Genotype=='CVN')

filteredNames <- meta_table$SampleID

abund_table <- abund_table_Base %>% select(Taxa, filteredNames)
abund_table <- as.data.frame(abund_table)
row.names(abund_table) <- abund_table$Taxa
abund_table$Taxa <- NULL

abund_table <- t(abund_table)

at_sums <- rowSums(abund_table)
abund_table <- abund_table[at_sums>0,]

meta_table$Variable <- meta_table$Cage

doBeta(title='Naive vs Vehicle', outfile='vegan/naive_vehicle_bray.png', label='Treatment', distance = 'bray', distanceName = 'Bray Curtis')
doBeta(title='Naive vs Vehicle', outfile='vegan/naive_vehicle_jaccard.png', label='Treatment', distance = 'jaccard', distanceName = 'Jaccard')
doBeta(title='Cage', outfile='vegan/cage_bray.png', label='Cage', distance = 'bray', distanceName = 'Bray Curtis')




#doBeta(title='CVN MvF 270 Days', outfile='vegan/CVN_MvF_270_bray.png', label='Sex', distance = 'bray', distanceName = 'Bray Curtis')
#doBeta(title='CVN vs WT 270 Days', outfile='vegan/CVNvWT_270_bray.png', label='Genotype', distance = 'bray', distanceName = 'Bray Curtis')
#doBeta(title='CVN vs WT 270 Days', outfile='vegan/CVNvWT_270_jaccard.png', label='Genotype', distance = 'jaccard', distanceName = 'Jaccard')

#doBeta(title='HN tMCAO vs Sham 28 Days', outfile="vegan/HN_tMCAOvSham_28_bray.png", label='Surgery', distance='bray', distanceName="Bray Curtis")
#doBeta(title='HN tMCAO vs Sham 28 Days', outfile="vegan/HN_tMCAOvSham_28_jaccard.png", label='Surgery', distance='jaccard', distanceName="Jaccard")
#doBeta(title='WT vs HN tMCAO 14 Days', outfile="vegan/WTvHN_tMCAO_14_chao.png", label='Genotype', distance='chao', distanceName="Chao")

#doBeta(title='WT vs HN tMCAO 28 Days', outfile="vegan/WTvHN_tMCAO_28_bray.png", label='Genotype', distance='bray', distanceName="Bray Curtis")
#doBeta(title='WT vs HN tMCAO 28 Days', outfile="vegan/WTvHN_tMCAO_28_jaccard.png", label='Genotype', distance='jaccard', distanceName="Jaccard")



# Need to set up meta_table and abund_table first
doBeta <- function(title="", outfile, label="",  distance='bray', distanceName="Bray Curtis"){

	print(paste("outfile is", outfile))

#Get MDS stats
sol<-metaMDS(abund_table, distance = distance, k = 2, trymax = 50)

NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],Variable=as.factor(meta_table$Variable))

#Get spread of points based on variable
plot.new()
ord<-ordiellipse(sol, as.factor(meta_table$Variable) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()

# Define Function
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#Generate ellipse points
df_ell <- data.frame()
for(g in levels(NMDS$Variable)){
  if(g!="" && (g %in% names(ord))){
    
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Variable==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                  ,Variable=g))
  }
}

#Generate mean values from NMDS plot grouped on Genotype
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Variable),mean)

# Beta Stat

bc.d <- vegdist(abund_table, method=distance)
ad=adonis(bc.d~Variable, data=meta_table, permutations=999)
pval <- ad$aov.tab$`Pr(>F)`[1]

# Plot

#shape_values<-seq(1,11)

p<-ggplot(data=NMDS,aes(x,y,colour=Variable))
#p<-p+ annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)
p<-p+ geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=1)
p <- p + geom_point()
#p <- p + labs(title = title, subtitle = paste('P-value:',pval), caption=distanceName, color=legend)
p <-  p + labs(title = title, subtitle = paste('P-value:',pval), caption=distanceName, color=label)
#p<-p+geom_point(aes(shape=Depth))+scale_shape_manual(values=shape_values)+theme_bw() 
png(outfile)
print(p)
dev.off()

}



#braycurtis.mds <- metaMDS()
#plot(braycurtis.mds, type="t")
#plot(sol)

