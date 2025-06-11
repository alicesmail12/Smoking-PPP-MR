# Packages
library(tidyverse)
library(TwoSampleMR)
library(ggnewscale)
library(lemon)

# Set wd
setwd('./Proteome-PPP-Harmonised/')
theme <- theme_classic() + theme(text=element_text(size=12, family='Franklin Gothic Book'))
filelist = list.files(pattern = ".*.txt")

# FUNCTIONS ##########################################################################################
# MR Scatter function
mr_scatter_plot <- function(mr_results, dat){
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d){
    d <- plyr::mutate(d)
    
    # Generate blank plot if not enough SNPs
    if(nrow(d) < 2 | sum(d$mr_keep) == 0){return(blank_plot("Insufficient number of SNPs"))}
    
    # Get rows where mr_keep=TRUE
    d <- subset(d, mr_keep)
    
    # Get exposure and outcome info
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    
    # Perform MR Egger regression
    if("MR Egger" %in% mrres$method)
    {temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
    mrres$a[mrres$method == "MR Egger"] <- temp$b_i}
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
    mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i}
    
    # Plot scatterplot
    ggplot(data=d, aes(x=beta.exposure, y=beta.outcome, colour=colour)) +
      theme_bw() + 
      theme(text = element_text(size = 12, family='Franklin Gothic Book'))+
      geom_hline(yintercept=0, colour='grey')+
      geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
      geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      geom_point(size=3, aes(colour=colour)) + guides(colour = "none") + scale_colour_manual(values=c("#8AB42F", "black"))+
      new_scale_color() +
      geom_abline(data=mrres,aes(intercept=a, slope=b, colour=method), show.legend=TRUE, linewidth=1) +
      scale_colour_manual(values=c("#FF7540", "#CCEA8B", "#FAE36F", "#79ACC5", "black"),
                          labels=c("IVW (MR Effects)", "MR-Egger", "Simple Median", "Weighted Median")) +
      labs(colour=NULL, x=paste("SNP \u03B2 effect on", d$exposure[1]), y=paste("SNP \u03B2 effect on", d$outcome[1])) +
      theme(legend.position=c(0.8,0.2),
            legend.direction="vertical",
            legend.title=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text=element_text(size = 20, family='Franklin Gothic Book'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 

  })
  mrres
}

# Modified MR scatter plot
mr_scatter_plot2 <- function(mr_results, dat){
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d){
    d <- plyr::mutate(d)
    
    # Generate blank plot if not enough SNPs
    if(nrow(d) < 2 | sum(d$mr_keep) == 0){return(blank_plot("Insufficient number of SNPs"))}
    
    # Get rows where mr_keep=TRUE
    d <- subset(d, mr_keep)
    
    # Get exposure and outcome info
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    
    # Perform MR Egger regression
    if("MR Egger" %in% mrres$method)
    {temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
    mrres$a[mrres$method == "MR Egger"] <- temp$b_i}
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
    mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i}
    
    # Plot scatterplot
    ggplot(data=d, aes(x=beta.exposure, y=beta.outcome, colour=colour)) +
      theme_bw() + 
      theme(text = element_text(size = 12, family='Franklin Gothic Book'))+
      geom_vline(xintercept=0, colour='grey')+
      geom_hline(yintercept=0, colour='grey')+
      new_scale_color() +
      geom_abline(data=mrres,aes(intercept=a, slope=b, colour=method), show.legend=TRUE, linewidth=1, alpha=0.75) +
      scale_colour_manual(values=c("#FF7540", "#CCEA8B", "#FAE36F", "#79ACC5"),
                          labels=c("IVW (MR Effects)", "MR-Egger", "Simple Median", "Weighted Median")) +
      labs(colour=NULL, x=paste("SNP \u03B2 effect on", d$exposure[1]), y=paste("SNP \u03B2 effect on", d$outcome[1])) +
      new_scale_color() +
      geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="black", width=0, alpha=0.65) +
      geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="black", height=0, alpha=0.65) +
      new_scale_color() + scale_colour_manual(values=c("black", "#79ACC5")) +
      geom_point(size=2.5, aes(colour=colour, fill=colour)) + guides(colour = "none", fill='none') + 
      theme(legend.position=c(0.8,0.2),
            legend.direction="vertical",
            legend.title=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text=element_text(size = 20, family='Franklin Gothic Book'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.box.background = element_rect(colour = "black", linewidth=1)) 
  })
  mrres
}

# GATHER IV FILES ####################################################################################
# Read one file
x <- read.table(filelist[1], header=TRUE, fill=TRUE, sep=',')
x <- x %>% select(SNP, exposure, chr.outcome, pos.outcome, beta.exposure, pval.exposure)
columns <- colnames(x)

# rbind all files
for (file in 2:length(filelist)){
  y<-read.csv(filelist[file], fill=TRUE, sep=',', header=TRUE)
  y <- y %>% select(SNP, exposure, chr.outcome, pos.outcome, beta.exposure, pval.exposure)
  print(y)
  x <- rbind(x,y)
}

# TNFRSF4 RESULTS ####################################################################################
# Get TNFRSF4 SNPs
TNFRSF4Harm <- read.table('Harm_TNFRSF4_PPP_MHC_Removed_Manual.txt', header=TRUE, fill=TRUE, sep=',') %>% 
  filter(mr_keep == TRUE)
x %>% count(SNP, chr.outcome, pos.outcome) %>% filter(SNP %in% TNFRSF4Harm$SNP)

# MR
res <- mr(TNFRSF4Harm, method_list = c('mr_egger_regression', 'mr_simple_median', 'mr_weighted_median', 'mr_ivw_mre'))
TNFRSF4Harm$colour <- ifelse(TNFRSF4Harm$SNP %in% c('rs11066309', 'rs597808', 'rs7968960'), 'Pleiotropic', 'NA')
TNFRSF4Harm$label <- ifelse(TNFRSF4Harm$SNP %in% c('rs11066309', 'rs597808', 'rs7968960'), 'Pleiotropic', NA)

# Get results
beta <- res[res$method == 'Inverse variance weighted (multiplicative random effects)',]$b
se <- res[res$method == 'Inverse variance weighted (multiplicative random effects)',]$se
pval <- res[res$method == 'Inverse variance weighted (multiplicative random effects)',]$pval

# Transform
padj <- pval*45
OR <- exp(beta)

ci <- exp(beta)*se
lower <- exp(beta)-ci*qnorm(0.975)
upper <- exp(beta)+ci*qnorm(0.975)

# Scatter plot
scatter <- mr_scatter_plot2(res, TNFRSF4Harm)
scatter[[1]] 

# Save plot
ggsave(filename='PleiotropicSNPsTNFRSF4.png', plot=scatter[[1]], width=989/3, height=357/3, units='mm', dpi = 300)

# Single SNP
singleSNP <- mr_singlesnp(TNFRSF4Harm)
mr_forest_plot(singleSNP)[[1]] + theme_classic()+ theme(text = element_text(size = 12, family='Franklin Gothic Book'))

# TNFRSF4 IV PLEIOTROPY ##############################################################################
# Write SNPs to file
lapply(c(TNFRSF4Harm$SNP), write, "SNPList.txt", append=TRUE)

# Get all files
setwd('./ProteomeSNPList/')
filelist = list.files(pattern = ".*.txt")

# Read one file
x <- read.table(filelist[1], header=FALSE, fill=TRUE, sep='\t')
colnames(x) <- c('CHR', 'POS', 'NAME', 'SNP', 'A1', 'A2', 'BETA', 'PVAL', 'LOG10PVAL', 'SE', 'N', 'MAF')
x$EXP <- gsub('.SNPs.txt', '', filelist[1])
x <- x %>% select(SNP, EXP, CHR, POS, BETA, PVAL)

# rbind all files
for (file in 2:length(filelist)){
  y<-read.table(filelist[file], fill=TRUE, sep='\t', header=FALSE)
  colnames(y) <- c('CHR', 'POS', 'NAME', 'SNP', 'A1', 'A2', 'BETA', 'PVAL', 'LOG10PVAL', 'SE', 'N', 'MAF')
  y$EXP <- gsub('.SNPs.txt', '', filelist[file])
  y <- y %>% select(SNP, EXP, CHR, POS, BETA, PVAL)
  print(y)
  x <- rbind(x,y)
}

# Add TNFRSF4 values
filt <- TNFRSF4Harm %>% select(SNP, exposure, chr.outcome, pos.outcome, beta.exposure, pval.exposure)
colnames(filt) <- c('SNP', 'EXP', 'CHR', 'POS', 'BETA', 'PVAL')
heat <- rbind(x, filt)

# Add PPP values
filt <- TNFRSF4Harm %>% select(SNP, outcome, chr.outcome, pos.outcome, beta.outcome, pval.outcome)
colnames(filt) <- c('SNP', 'EXP', 'CHR', 'POS', 'BETA', 'PVAL')
heat <- rbind(heat, filt)

# Set exposures as factors
filt <- heat %>% filter(SNP %in% TNFRSF4Harm$SNP) %>% slice(order(EXP != 'TNFRSF4' & EXP != 'PPP'))
filt$EXP <- factor(filt$EXP, levels=unique(filt$EXP))

# Get y-axis labels
filtLevels <- filt %>% filter(EXP == 'TNFRSF4') %>% arrange(-BETA) %>% select(SNP)
z <- filt %>% filter(EXP == 'TNFRSF4') %>% arrange(PVAL)
median(z$PVAL)

# Pval labels
filt$label <- ifelse(filt$PVAL<1e-08, '*', NA)
filt$sig <- ifelse(filt$PVAL<1e-08, 1, 0)

# Shared pleiotropy
filtSum <- filt %>% group_by(SNP) %>% summarise(sigNum=sum(sig)) %>% as.data.frame()
filtSum <- filtSum[match(rev(filtLevels$SNP), filtSum$SNP),]

# Minus 1 for TNFRSF4
filtSum$sigNum <- filtSum$sigNum-1
nums <- filtSum$sigNum

# TNFRSF4 IV PLEIOTROPY PLOT #########################################################################
# Heatmap
p <- ggplot(filt, aes(x=factor(EXP), y=SNP, fill=BETA)) + 
  
  # Mapping
  geom_tile(alpha=1) +
  xlab("") + ylab("")+
  theme +
  
  # Colour
  scale_fill_gradientn(limits=c(-0.15, 0.15), 
                       colours=c('#57818c', 'white', '#ad2424'),
                       na.value = '#57818c',
                       guide=guide_colorbar(frame.colour="black", 
                                            ticks.colour="black", 
                                            alpha=0.85, 
                                            title='Beta Effect'))+
  
  # Axis labels
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  
  # Labels
  geom_text(
    label=filt$label, 
    nudge_x = 0.0, nudge_y = -0.3, 
    check_overlap = F,
    size=6, colour='white') +
  
  # First two cols
  annotate("rect", xmin=c(0.5), xmax=c(1.5), ymin=c(0.5), ymax=c(14.5), 
           colour="black", fill="transparent", size=0.5)+
  annotate("rect", xmin=c(1.5), xmax=c(2.5), ymin=c(0.5), ymax=c(14.5), 
           colour="black", fill="transparent", size=0.5)+

  # Scale
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0), limits=rev(filtLevels$SNP))+
  
  # Bold
  theme(axis.text.y = element_text(family=c('Franklin Gothic Book', 'Franklin Gothic Medium',
                                            rep('Franklin Gothic Book', 2), 'Franklin Gothic Medium',
                                            rep('Franklin Gothic Book', 7), 'Franklin Gothic Medium'))) +
  
  # Table
  annotate("rect", xmin=c(rep(47)), xmax=c(rep(50)), 
           ymin=c(seq(0.5, 13.5, by=1)), ymax=c(seq(1.5, 14.5, by=1)), 
           colour="white", fill="white", size=0.5)+
  annotate("text", x = 47, y = seq(1, 14, by=1), label=nums, 
           family='Franklin Gothic Book', hjust = 0) +
  
  # Legend
  theme(legend.position = "right", axis.ticks.y=element_blank(),
        axis.title.y = element_text(colour = "#8c8c8c"),
        axis.title.x = element_text(colour = "#8c8c8c"),
        legend.box.spacing = unit(0, "pt")) +
  
  # Axes
  coord_flex_cart(ylim=c(0.5, 16), left=capped_vertical(capped='top', gap = 0)) +
  coord_flex_cart(xlim=c(0.5, 53), bottom=capped_horizontal(capped='right', gap = 0)) +
  
  # Border
  annotate("rect", xmin=c(0.5), xmax=c(46.49), ymin=c(0.5), ymax=c(14.49), 
           colour="black", fill="transparent", size=0.75) +
  
  # y axis label
  annotate(geom = "text", x = 49.5, y = 7.5, label = "Number of Proteins Sharing pQTL", color = "#8c8c8c",
           angle = 270, family='Franklin Gothic Book') +
  ylab('SNPs') + xlab('Protein')

# Plot
ggsave(filename='SummaryPlot.png', plot=p, width=989/4, height=357/4, units='mm', dpi = 300)

# END ################################################################################################


