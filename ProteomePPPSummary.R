# Install packages
library(ggplot2)
library(dplyr)
library(extrafont)
library(ggrepel)
library(tibble)

# Set working directory
setwd("./MR/ProteomePPP/Results")

# Read results csv
SmkInitProteomeRes <- read.csv("ProteomePPPRes.csv", header=TRUE)

# Get valid proteins list
ValidProteins <- read.csv('ValidProteins.csv')
ProteomePPPRes <- ProteomePPPRes %>% filter(exposure %in% ValidProteins[[1]])

# Only retain proteins with threshold measured at p<1e-8 (45)
ProteomePPPRes <- ProteomePPPRes %>% filter(SNP_pval == '1e-08')
nrow(ProteomePPPRes)

# Assign IVW MRE & FE results to one IVW column
ProteomePPPRes$Beta <- ProteomePPPRes$IVW_FE_Beta 
ProteomePPPRes$Pval <- ProteomePPPRes$IVW_FE_pval 
for (i in 1:nrow(ProteomePPPRes)) {
  if (is.na(ProteomePPPRes$Beta[i])) {
    ProteomePPPRes$Beta[i] <- ProteomePPPRes$IVW_MRE_Beta[i]
    ProteomePPPRes$Pval[i] <- ProteomePPPRes$IVW_MRE_pval[i]
  }
}

# Set default shape value to 'Failed'
ProteomePPPRes$shape <- 'Failed'

# Check if all MR tests have the same direction of effect and change shape to 'Passed'
for (i in 1:nrow(ProteomePPPRes)) {
  if ((ProteomePPPRes$Beta[i]>0 & ProteomePPPRes$MR_Egger_Beta[i]>0 & ProteomePPPRes$Simple_Median_Beta[i]>0 & ProteomePPPRes$Weighted_Median_Beta[i]>0) | 
      (ProteomePPPRes$Beta[i]<0 & ProteomePPPRes$MR_Egger_Beta[i]<0 & ProteomePPPRes$Simple_Median_Beta[i]<0 & ProteomePPPRes$Weighted_Median_Beta[i]<0)){
    ProteomePPPRes$shape[i] <- 'Passed'
  }
}

# Correct IVW p-value using BH method
ProteomePPPRes <- ProteomePPPRes %>% arrange(Pval)
ProteomePPPRes <- rowid_to_column(ProteomePPPRes, "rank")
for (i in 1:nrow(ProteomePPPRes)) {
  ProteomePPPRes$Pval_BH[i] <- (ProteomePPPRes$Pval[i] * nrow(ProteomePPPRes))/ProteomePPPRes$rank[i]
}

# Using BH-corrected p-value, identify proteins that are significantly effected (2)
sig <- nrow(ProteomePPPRes %>% filter(Pval_BH < 0.05)) 

# Correct pleiotropy p-value using BH method and set shape to 'Failed' for proteins where significant pleiotropy is detected
ProteomePPPRes <- ProteomePPPRes %>% arrange(Pleio_pval)
ProteomePPPRes <- rowid_to_column(ProteomePPPRes, "rank_pleio")
for (i in 1:nrow(ProteomePPPRes)) {
  ProteomePPPRes$Pleio_pval_BH[i] <- (ProteomePPPRes$Pleio_pval[i] * nrow(ProteomePPPRes))/ProteomePPPRes$rank_pleio[i]
  if (ProteomePPPRes$Pleio_pval_BH[i]<0.05) {
    ProteomePPPRes$shape[i] <- 'Failed'
  }
}

# Identify significant proteins that pass both sensitivity tests (2)
ProteomePPPResPassed <- ProteomePPPRes %>% filter(Pval_BH < 0.05) %>% filter(shape == 'Passed')
nrow(ProteomePPPResPassed)

# Add significance labels for two alpha values
ProteomePPPRes$colour <- ifelse(ProteomePPPRes$Pval_BH<(0.05), "P<0.05", 'NS')
ProteomePPPRes$colour <- ifelse(ProteomePPPRes$Pval_BH<(0.005), "P<0.005", ProteomePPPRes$colour)

# Label 2 proteins with smallest p-value
ProteomePPPRes <- ProteomePPPRes %>% arrange(Pval_BH)
ProteomePPPRes$label <- NA
for (i in (1:2)) {ProteomePPPRes$label[i] <- ProteomePPPRes$exposure[i]}

# Plot
volcanoplot <- ggplot(ProteomePPPRes, aes(x=Beta, y=-log10(Pval_BH), colour=colour, label=label, shape=shape))+
  geom_point(size=5)+
  theme_classic()+
  theme(text=element_text(family='Franklin Gothic Book', size=18))+
  geom_hline(yintercept=-log10(0.05), linewidth=1, colour='#7db2c5')+
  geom_hline(yintercept=-log10(0.005), linewidth=1, colour='#f26324')+
  scale_colour_manual(values=c('#bdbdbd', '#7db2c5', '#f26324'),
                      labels=c('NS', 'P<0.05', 'P<0.005'),
                      limits=c('NS', 'P<0.05', 'P<0.005'),
                      drop=FALSE)+
  scale_shape_manual(values=c(17, 19))+
  labs(y=bquote(-log[10](P)), x='\u03B2', colour='Adjusted P-value (BH)', shape='Sensitivity Tests')+
  theme(legend.title = element_text(family='Franklin Gothic Medium'))+
  geom_text_repel(colour='black', box.padding = 0.55, size=4)

# Save
ggsave(plot=scatter, filename="ProteomePPPRes.png", width=10, height=8, units='in', dpi=500)
