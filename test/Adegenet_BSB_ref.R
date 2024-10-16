###################
## ADGENET 2.1.1 ##
###################
rm(list = ls())
#install.packages("adegenet", dep=TRUE)

library("adegenet")
library("ape") 
library("RColorBrewer")
library("poppr") # create trees and treat missing data

#show Adegenet version
packageDescription("adegenet", fields = "Version") 

#set working directory
#setwd("/Users/remygatins/Documents/Drive/Projects/2021_Black Sea Bass/RADs/popmap_ref/BSB_all")
#setwd("/Users/remygatins/Documents/Drive/Projects/2021_Black Sea Bass/RADs/popmap_ref/BSB_all_r0.8")
#setwd("/Users/remygatins/Documents/Drive/Projects/2021_Black Sea Bass/RADs/popmap_ref/BSB_all_r0.8_maf0.05")
#setwd("/Users/remygatins/GoogleDrive_gmail/Work/Projects/2021_Black\ Sea\ Bass/RADs/popmap_ref/BSB_all_r0.8_maf0.01")
setwd("/Users/remygatins/GoogleDrive_gmail/Work/Projects/2021_Black\ Sea\ Bass/RADs/popmap_ref/BSB_all_r0.8_maf0.01_wsnp")
setwd("/Users/remygatins/GoogleDrive_gmail/Work/Projects/2021_Black\ Sea\ Bass/RADs/popmap_ref/BSB_all_r0.8_maf0.01_wsnp_noNC")


#Load .vcf data output from populations (STACKS)
#install.packages("vcfR")
library(vcfR)
vcf <- read.vcfR(file = "populations.snps.vcf")
#vcf <- read.vcfR(file = "filtered_ind_2.recode.vcf")
#vcf <- read.vcfR(file = "filtered.recode.vcf")
data <- vcfR2genlight(vcf)
data


#__________________
#Label populations 
#------------------

#popmap1
pop(data)<- c(rep("NM",20), rep("MD",20), rep("ME",18), rep("NC",13), rep("NJ",17), rep("SN",29)) 

#BSB_noNC
pop(data)<- c(rep("NM",20), rep("MD",20), rep("ME",18), rep("NJ",17), rep("SN",29)) 


Dgen <- dist.gene(data)
Dgeo <- dist(data)
ploidy(data) <- 2

#--------------------------------
# Convert vcf file to structure
#--------------------------------
#if (!require("devtools")) install.packages("devtools")
#devtools::install_github("thierrygosselin/radiator")
library(radiator)

#To verify your file is detected by radiator as the correct format:
radiator::detect_genomic_format(
  data="filtered_ind_2.recode.vcf")

strata <- radiator::read_strata(strata = "BSB.strata.filt.txt")

BSB_data <- genomic_converter(
  data = "filtered_ind_2.recode.vcf", 
  strata = "BSB.strata.filt.txt",
  output = c("structure","genepop"))

structure <- BSB_data

write_structure(BSB_data$tidy.data, filename = "structure")

tidy.data <-  tidy_genlight(data)

tidy.data$GT <- tidy.data$GT_BIN

tidy.vcf.data <-  tidy_vcf(data="filtered_ind_2.recode.vcf",
                       strata = "BSB.strata.filt.txt")

#rename colummn GT_BIN to GT
library(dplyr)
tidy.vcf.data_2<- tidy.vcf.data %>% 
                      rename(GT= GT_BIN)

write_structure(tidy.vcf.data_2, filename = "str")





##############
# Check data #
##############
#plot missing data
glPlot(data, posi="topleft")
glPlot(data, col=bluepal(6))

#plot allele frequencies
myFreq <- glMean(data)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies", main="Distribution of (second) allele frequencies")
temp <- density(myFreq) 
lines(temp$x, temp$y*1.8,lwd=3)


#map missing data across loci
temp <- density(glNA(data)) 
plot(temp, type="n", xlab="Position in the alignment", main="Location of the missing values (NAs)")
polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3)) 
points(glNA(SRE_light), rep(0, nLoc(SRE_light)), pch="|", col="blue")



#########################
#Create a Distance tree #
#########################

library(poppr)
tree <- aboot(data, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

#add colors to differentiate between populations
cols <- brewer.pal(n = nPop(data), name = "Dark2")
cols <- c(brewer.pal(n = 8, name="Dark2"), brewer.pal(name="Paired", n=12),brewer.pal(name="Pastel1", n=5) )
plot.phylo(tree, cex = 0.4, font = 2, adj = 0, tip.color =  cols[pop(data)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend('topleft', legend = c("H. clarionensis","H. limbaughi","H. passer", "Hybrid"), fill = cols, border = FALSE, bty = "n", cex = 1)
#legend('topleft', legend = c("Clarion","Los Cabos","Roca Partida", "San Benedicto", "Socorro"), fill = cols, border = FALSE, bty = "n", cex = 1)
legend('topleft', legend = c("Clarion","San Benedicto", "Socorro"), fill = cols, border = FALSE, bty = "n", cex = 1)
legend('topleft', legend = c("Mazatlán","B. Magdalena","El Faro","Clipperton","Galapagos","Panama","Cabo Pulmo","Angel Guarda","B. Muertos","Santa Cruz","La Paz","San Pedro Martir","Zihuatanejo","Los Cabos"), fill = cols, border = FALSE, bty = "n", cex = 1)


axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")


###########################
#Minimum Spanning networks#
###########################

library(igraph)

data.dist <- bitwise.dist(data)
data.msn <- poppr.msn(data, data.dist, showplot = FALSE, include.ties = T)

node.size <- rep(2, times = nInd(data))
names(node.size) <- indNames(data)
vertex.attributes(data.msn$graph)$size <- node.size

set.seed(9)
plot_poppr_msn(data, data.msn , palette = brewer.pal(n = nPop(data), name = "Dark2"), gadj = 70)


######################
#Principal Components#
######################

data.pca <- glPca(data, nf = 4) #number of PC axes to retain (here, 4)
data.pca # prints summary

barplot(100*data.pca$eig/sum(data.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

scatter(data.pca, posi="topright", cex=0.5, alpha= 0.2)
title("PCA of data\n axes 1-2")

#Assess PCA using colors
myCol <- colorplot(data.pca$scores,data.pca$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey") 
add.scatter.eig(data.pca$eig[1:40],2,1,2, posi="bottomright", inset=.1, ratio=.2)
legend(10,-4, col=c("orange","red", "blue", "green"))

#View PCA by population
#To view the results of the PCA we can use the package ggplot2. 
#We need to convert the data frame that contains the principal components (data.pca$scores) 
#into the new object data.pca.scores. In addition, we will add the population values as a new 
#column in our data.pca.scores object, in order to be able to color samples by population.

data.pca.scores <- as.data.frame(data.pca$scores)
data.pca.scores$pop <- pop(data)
#write.csv(data.pca.scores, "data.pca_scores.csv")

#add the rest of your metadata to your pcs table, then reload containing all metadata


cols <- brewer.pal(n = nPop(data), name = "Dark2")
cols <- c(brewer.pal(n = 8, name="Dark2"), brewer.pal(name="Paired", n=12), brewer.pal(name="Set2", n=5))

library(ggplot2)
set.seed(9)
ggplot(data.pca.scores, aes(x=PC1, y=PC2, colour=pop)) +
            geom_point(size=2, alpha=0.6) +
            #stat_ellipse(level = 0.95, size = 1) +
            scale_color_manual(values = cols) +
            geom_hline(yintercept = 0) +
            geom_vline(xintercept = 0) +
            theme_bw()

#Variance explained of PC1 and PC2
vPC1 <- data.pca$eig[1]/sum(data.pca$eig)*100
vPC2 <- data.pca$eig[2]/sum(data.pca$eig)*100
vPC1
vPC2

ggplot(data.pca.scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=2, alpha=0.6) +
  #stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("PC1 (1.28% exp var)")+
  ylab("PC2 (1.27% exp var)")+
  labs(colour = "")+
  theme_bw()

data.pca.scores$pop <- factor(data.pca.scores$pop, levels = c("ME","NM","SN","NJ","MD","NC"))

ggplot(data.pca.scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=2, alpha=0.6) +
  #stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("PC1 (1.28% exp var)")+
  ylab("PC2 (1.27% exp var)")+
  labs(colour = "")+
  theme_bw()

write.csv(data.pca.scores, "data.pca.scores.csv")

data_noNC <- data %>%
  mutate(pop != NC)


########
# DAPC #
########

data.dapc1 <- dapc(data, n.pca = 60, n.da = 2)

data.dapc1$pop <- factor(data.pca.scores$pop, levels = c("ME","NM","SN","NJ","MD","NC"))


dscatter(data.dapc1, col = cols, cex = 2, legend = TRUE, clabel = F, cstar=1, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topright", cleg = 0.75, ratio.pca = 0.2, ratio.da = 0.2)

# find the number of clusters in your dataset
grp <- find.clusters(data, max.n.clust=7)
table(pop(data), grp$grp)
table.value(table(pop(data), grp$grp), col.lab=paste("inf", 1:6), row.lab=paste("ori", 1:6))


#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(data.dapc,col = cols, posi = 'top')

#transform data to comoplot in ggplot
dapc.results <- as.data.frame(data.dapc$posterior)
dapc.results$pop <- pop(data)
dapc.results$indNames <- rownames(dapc.results)

library(reshape2)
dapc.results <- melt(dapc.results)
colnames(results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = cols) +
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
        