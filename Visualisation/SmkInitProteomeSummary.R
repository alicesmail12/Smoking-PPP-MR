# Install packages
library(ggplot2)
library(dplyr)
library(extrafont)
library(ggrepel)
library(tibble)

# Set working directory
setwd("./MR/SmkInitProteome/Results")

# Read results csv
SmkInitProteomeRes <- read.csv("SmkInitProteomeRes.csv", header=TRUE)

# Identify the total number of proteins tested in this analysis (160)
nrow(SmkInitProteomeRes) 

# Assign IVW MRE & FE results to one IVW column
SmkInitProteomeRes$Beta <- SmkInitProteomeRes$IVW_FE_Beta 
SmkInitProteomeRes$Pval <- SmkInitProteomeRes$IVW_FE_pval 
for (i in 1:nrow(SmkInitProteomeRes)) {
  if (is.na(SmkInitProteomeRes$Beta[i])) {
    SmkInitProteomeRes$Beta[i] <- SmkInitProteomeRes$IVW_MRE_Beta[i]
    SmkInitProteomeRes$Pval[i] <- SmkInitProteomeRes$IVW_MRE_pval[i]
  }
}

# Set default shape value to 'Failed'
SmkInitProteomeRes$shape <- 'Failed'

# Check if all MR tests have the same direction of effect and change shape to 'Passed'
for (i in 1:nrow(SmkInitProteomeRes)) {
  if ((SmkInitProteomeRes$Beta[i]>0 & SmkInitProteomeRes$MR_Egger_Beta[i]>0 & SmkInitProteomeRes$Simple_Median_Beta[i]>0 & SmkInitProteomeRes$Weighted_Median_Beta[i]>0) | 
      (SmkInitProteomeRes$Beta[i]<0 & SmkInitProteomeRes$MR_Egger_Beta[i]<0 & SmkInitProteomeRes$Simple_Median_Beta[i]<0 & SmkInitProteomeRes$Weighted_Median_Beta[i]<0)){
    SmkInitProteomeRes$shape[i] <- 'Passed'
  }
}

# Correct IVW p-value using BH method
SmkInitProteomeRes <- SmkInitProteomeRes %>% arrange(Pval)
SmkInitProteomeRes <- tibble::rowid_to_column(SmkInitProteomeRes, "rank")
for (i in 1:nrow(SmkInitProteomeRes)) {
  SmkInitProteomeRes$Pval_BH[i] <- (SmkInitProteomeRes$Pval[i] * nrow(SmkInitProteomeRes))/SmkInitProteomeRes$rank[i]
}

# Using BH-corrected p-value, identify proteins that are significantly effected (53)
nrow(SmkInitProteomeRes %>% filter(Pval_BH < 0.05)) 

# Correct pleiotropy p-value using BH method and set shape to 'Failed' for proteins where significant pleiotropy is detected
SmkInitProteomeRes <- SmkInitProteomeRes %>% arrange(Pleio_pval)
SmkInitProteomeRes <- tibble::rowid_to_column(SmkInitProteomeRes, "rank_pleio")
for (i in 1:nrow(SmkInitProteomeRes)) {
  SmkInitProteomeRes$Pleio_pval_BH[i] <- (SmkInitProteomeRes$Pleio_pval[i] * nrow(SmkInitProteomeRes))/SmkInitProteomeRes$rank_pleio[i]
  if (SmkInitProteomeRes$Pleio_pval_BH[i]<0.05) {
    SmkInitProteomeRes$shape[i] <- 'Failed'
  }
}

# Identify significant proteins that pass both sensitivity tests (48)
SmkInitProteomeResPassed <- SmkInitProteomeRes %>% filter(Pval_BH < 0.05) %>% filter(shape == 'Passed')
nrow(SmkInitProteomeResPassed)

# Add significance labels for two alpha values
SmkInitProteomeRes$colour <- ifelse(SmkInitProteomeRes$Pval_BH<(0.05), "P<0.05", 'NS')
SmkInitProteomeRes$colour <- ifelse(SmkInitProteomeRes$Pval_BH<(0.005), "P<0.005", SmkInitProteomeRes$colour)

# Label 5 proteins with smallest p-value
SmkInitProteomeRes <- SmkInitProteomeRes %>% arrange(Pval_BH)
SmkInitProteomeRes$label <- NA
for (i in (1:5)) {SmkInitProteomeRes$label[i] <- SmkInitProteomeRes$outcome[i]}

# Plot
volcanoplot <- ggplot(SmkInitProteomeRes, aes(x=Beta, y=-log10(Pval_BH), colour=colour, label=label, shape=shape))+
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

ggsave(plot=volcanoplot, filename="SmkInitProteomeRes.png", width=10, height=8, units='in', dpi=500)

# Get list of cytokines for Step 2
write.csv(SmkInitProteomeResPassed$outcome, "ValidProteins.csv", row.names=FALSE)

