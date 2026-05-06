## Outlier detection

```r
################
### pcadapt ###
################

library(pcadapt)

vcf_path <- "minDP10_maxmiss0.7_filtInd.recode.vcf"
input_BSB <- read.pcadapt(vcf_path, type = "vcf")

#Choose number of K
x <- pcadapt(input = input_BSB, K = 20) 
plot(x, option = "screeplot")
#Choose the #PC to the left of the elbow or straight line... so here K=1

#run pcadapt with best K
obj <- pcadapt(input_BSB, K = 1)
summary(obj)

# manhattan plot
plot(obj, option = "manhattan") # ok nice, we have a base plot here

#check the expected uniform distribution of the p-values using a Q-Q plot
plot(obj,option="qqplot")

#histogram of p-values
hist(obj$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
#The excess of small p-values indicates the presence of outliers.

#apply bonferroni correction 
padj <- p.adjust(obj$pvalues,method="bonferroni")
alpha <- 0.01
outliers <- which(padj < alpha)
n_outliers <- length(outliers) 
cat("\nOutlier SNPs:",n_outliers,"\n") #155 outliers
head(outliers)


outliers_bsb <- vcf@fix[outliers, c("CHROM", "POS", "ID")]
outliers_bsb <- as.data.frame(outliers_bsb, stringsAsFactors = FALSE)
outliers_bsb$LocusID <- paste0(outliers_bsb$CHROM, ":", outliers_bsb$POS)

write.csv(outliers_bsb, "pcadapt_outlier_snps_meta.csv", row.names = FALSE)
writeLines(outliers_bsb$LocusID, "outlier_locus_ids.txt")





vcf_out <- vcf
vcf_out@fix <- vcf@fix[outliers,,drop=FALSE]
vcf_out@gt <- vcf@gt[outliers,,drop=FALSE]
write.vcf(vcf_out,"BSB_outliers.vcf.gz")
system(paste0("gunzip -f ","/Users/madelineeppley/Desktop/manta26pub/mobular_outliers.vcf.gz"))




# Improve Manhattan Plot Visualization
# extract p-values
pvals <- obj$pvalues
neglog10_pvals <- -log10(pvals)

# Extract chromosome info from VCF
chrom_info <- vcf@fix[, "CHROM"]

# calculate Bonferroni threshold
bonferroni_threshold <- -log10(0.01 / length(pvals))
manhattan_data <- data.frame(
  SNP_index = 1:length(pvals),
  neglog10_pval = neglog10_pvals,
  significant = neglog10_pvals > bonferroni_threshold,
  chrom = chrom_info
)

# Create alternating color index by chromosome
chrom_ordered <- unique(chrom_info)  # keeps original order
manhattan_data$chrom_num <- as.integer(factor(manhattan_data$chrom, 
                                              levels = chrom_ordered))
manhattan_data$chrom_color <- factor(manhattan_data$chrom_num %% 2)  # alternates 0 and 1
manhattan_data$chrom_label <- gsub("Scaffold_", "", manhattan_data$chrom)

manhattan_data <- manhattan_data %>%
  filter(as.numeric(chrom_label) < 25) #only keep chromosomes 1-24 as the rest don't show anything and just clump up

# Calculate midpoint of each chromosome for x-axis labels
chrom_midpoints <- manhattan_data %>%
  group_by(chrom_label) %>%
  summarise(mid = mean(SNP_index))

scale_x_continuous(breaks = chrom_midpoints$mid, 
                   labels = chrom_midpoints$chrom_label)



#remove NA's from dataset
manhattan_data <- manhattan_data[!is.na(manhattan_data$neglog10_pval), ]

# Plot
manhattan_plot <- ggplot(manhattan_data, aes(x = SNP_index, y = neglog10_pval)) +
  geom_point(aes(color = chrom_color), size = 1, alpha = 0.7) +
  geom_hline(yintercept = bonferroni_threshold, linewidth = 0.8) +
  geom_point(data = subset(manhattan_data, significant == TRUE), 
             shape = 21, size = 2, fill = NA, color = "red") +
  #scale_color_manual(values = c("0" = "#5BA4CF", "1" = "#A6611A")) +
  scale_color_manual(values = c("0" = "#BABABA", "1" = "#0571B0")) +
  scale_x_continuous(breaks = chrom_midpoints$mid, 
                     labels = chrom_midpoints$chrom_label) +
  labs(x = "Pseudo-Chromosome",
       y = expression(-log[10](p-value)),
       #title = "Manhattan Plot of PCAdapt Results",
       ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8)  # rotate for readability
  )

manhattan_plot

ggsave("BSB_pcadapt_outliers.png", manhattan_plot, width = 8, height = 4, dpi = 600)
```

### BSB_filt_relate_r0.7

Outlier SNPs: 155 
<img width="1469" height="722" alt="image" src="https://github.com/user-attachments/assets/fbd0ef5e-058c-4c20-ba19-b9c94a700d64" />



  
  
  
