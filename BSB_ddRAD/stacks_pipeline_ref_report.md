Run reference pipeline using BSB_all_r0.8_maf0.01 and `--write-single-snp`

```
#write-single-snp
#populations -P $DIR/stacks/ -M $DIR/popmap/BSB_all --write-single-snp -r 0.80 --min-maf 0.01 --vcf --genepop --structure --fstats --hwe -t 10
```


```r
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


#Load .vcf data output from populations (STACKS)
#install.packages("vcfR")
library(vcfR)
vcf <- read.vcfR(file = "populations.snps.vcf")
#vcf <- read.vcfR(file = "filtered_ind_2.recode.vcf")
#vcf <- read.vcfR(file = "filtered.recode.vcf")
data <- vcfR2genlight(vcf)
data

```
<img width="538" alt="image" src="https://user-images.githubusercontent.com/26288352/228410995-9ddb102c-a61a-4526-8952-fbb49b879264.png">


```r
#__________________
#Label populations 
#------------------

#popmap1
pop(data)<- c(rep("NM",20), rep("MD",20), rep("ME",18), rep("NC",13), rep("NJ",17), rep("SN",29)) 


Dgen <- dist.gene(data)
Dgeo <- dist(data)
ploidy(data) <- 2

######################
#Principal Components#
######################

data.pca <- glPca(data, nf = 4) #number of PC axes to retain (here, 4)
data.pca # prints summary

barplot(100*data.pca$eig/sum(data.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

```
<img width="563" alt="image" src="https://user-images.githubusercontent.com/26288352/228411350-1ba145e9-ba60-42a1-8547-595d0e963853.png">


```r
#View PCA by population
#To view the results of the PCA we can use the package ggplot2. 
#We need to convert the data frame that contains the principal components (data.pca$scores) 
#into the new object data.pca.scores. In addition, we will add the population values as a new 
#column in our data.pca.scores object, in order to be able to color samples by population.

data.pca.scores <- as.data.frame(data.pca$scores)
data.pca.scores$pop <- pop(data)
write.csv(data.pca.scores, "data.pca_scores.csv")

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
```
<img width="565" alt="image" src="https://user-images.githubusercontent.com/26288352/228411896-b4beb7d2-d72c-411e-b2da-6dd7f7cbf7fb.png">


```r
#rearrange populations from North to South
data.pca.scores$pop <- factor(data.pca.scores$pop, levels = c("ME","NM","SN","NJ","MD","NC"))

ggplot(data.pca.scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=2, alpha=0.6) +
  #stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("PC1 (1.81% exp var)")+
  ylab("PC2 (1.50% exp var)")+
  labs(colour = "")+
  theme_bw()

```

<img width="685" alt="image" src="https://user-images.githubusercontent.com/26288352/228410934-7df51af5-dd1e-49ad-9746-9f4a5ec2529b.png">



Now filter out individuals from NC and re-run the analyses

```bash
#filter out NC samples
populations -P $DIR/stacks/ -M $DIR/popmap/BSB_noNC --write-single-snp -r 0.80 --min-maf 0.01 --vcf --genepop --structure --fstats --hwe -t 10

```
<img width="546" alt="image" src="https://user-images.githubusercontent.com/26288352/228420542-72654c3c-fe97-4df5-9a85-f6f86fcac641.png">


<img width="687" alt="image" src="https://user-images.githubusercontent.com/26288352/228420092-94a4e9ba-44cc-4f2c-8cef-5f864eb725f4.png">


# SNMF on LEA to get structure plot

run interactive mode

            srun -p lotterhos -N 1 --pty /bin/bash

1. Create chromosome map

```bash
module load miniconda3
conda activate plink
conda activate bcftools

bcftools view -H populations.snps.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > populations.snps.vcf.chrom-map.txt
```

2. Make ped file using this chromosome map

```bash
module load vcftools

vcftools --vcf populations.snps.vcf --out BSB_r0.8_maf0.01 --plink --chrom-map populations.snps.vcf.chrom-map.txt
```

```bash
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf populations.snps.vcf
	--chrom-map populations.snps.vcf.chrom-map.txt
	--out BSB_r0.8_maf0.01
	--plink

After filtering, kept 117 out of 117 Individuals
Writing PLINK PED and MAP files ...
	Read 39 chromosome mapping file entries.
Done.
After filtering, kept 24896 out of a possible 24896 Sites
Run Time = 1.00 seconds
```
This should have created a .ped and .map file. Now download this to your local computer or wherever you will run R from

In R, install LEA and convert from a .ped file to a geno file to run SNMF (We will first need to convert to a lfmm file and then to a geno file)

http://membres-timc.imag.fr/Olivier.Francois/LEA/tutorial.htm

```r
############
## SNMF ####
############
rm(list = ls())

#===== BSB =========
# 6 source populations
# 117 individuals
# 24,896 diploid loci
# ? env variables (none for now)

setwd("/Users/remygatins/GoogleDrive_gmail/Work/Projects/2021_Black\ Sea\ Bass/RADs/popmap_ref/BSB_all_r0.8_maf0.01_wsnp")

###################
##### LEA #########
###################

#Install LEA
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LEA")
library(LEA)


#Convert ped to lfmm and lfmm to geno
output = ped2lfmm("BSB_r0.8_maf0.01.ped")
#output = BSB_r0.8_maf0.01.lfmm (writen to the working directory)
# Create file:	"genotypes.geno".
output = lfmm2geno("BSB_r0.8_maf0.01.lfmm", "genotypes.geno")
```

Now choose the most likely number of clusters

```r
# Choose the number of clusters
obj.snmf = snmf("genotypes.geno", K = 1:6, ploidy = 2, entropy = T, alpha = 100, project = "new")
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)
```
<img width="941" alt="image" src="https://user-images.githubusercontent.com/26288352/228557667-ebffd7b5-10fa-45be-b5f6-90313a663098.png">

K=1 seems to be the most likely number of cluster. 

I will run K=3 to see the plot 

```r
obj.snmf = snmf("genotypes.geno", K = 3, alpha = 100, project = "new") 
qmatrix = Q(obj.snmf, K = 3)
barplot(t(qmatrix), col = c("orange","violet","lightgreen"), border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")
```
![image](https://user-images.githubusercontent.com/26288352/228557976-bc1b753a-2310-4a57-bdfa-d130e56677f0.png)






