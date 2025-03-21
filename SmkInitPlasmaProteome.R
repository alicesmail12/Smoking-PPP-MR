# Install packages
message('\nInstalling packages...')
source("MR_Functions.R")
library(tidyverse)
library(LDlinkR)
library(TwoSampleMR)
library(ieugwasr)
library(RadialMR)
library(glue)
library(extrafont)
library(gsubfn)
library(gwasrapidd)
library(MRPRESSO)
library(data.table)
library(R.utils)
message('Packages installed\n')

################################ EXPOSURE #######################################
# Get smoking initiation SNPs
message('\nReading smoking file...')
options(datatable.fread.datatable=FALSE)
smoking <- data.table::fread('./GWAS_Smoking/GSCANSmkInit2022EUR.txt', fill=TRUE, header=TRUE)
message('Completed\n')

######################### IDENTIFYING EXPOSURE IVS ##############################
# Get significant SNPs and set phenotype
message('Formatting smoking file...')
smoking$CHR <- as.factor(str_split(smoking$CHR, 'chr', simplify = TRUE)[,2])
SNPpval <- 1e-8
smokingExp <- smoking %>% filter(P <= SNPpval)
smokingExp$Phenotype <- 'Smoking'

# Remove MHC region
# GRCh38.p14 Location: chr6:28,510,120-33,480,577
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC

# Get MHC region
MHC <- smokingExp %>%
  arrange(CHR, POS) %>%
  filter(CHR == 6) %>%
  filter(POS >= 28510120 & POS <= 33480577)

# Create MHC SNP list
MHCSNP <- MHC %>% dplyr::select(RSID)

# Remove MHC SNPs from GWAS
smokingExpNoMHC <- smokingExp %>% dplyr::filter(!(RSID %in% unlist(MHCSNP)))

######################### FORMATTING EXPOSURE DATA ##############################
# Rename columns 
smokingExp <- smokingExpNoMHC %>%
  rename('BP' = 'POS',
         'effect_allele' = 'EFFECT_ALLELE',
         'other_allele' = 'OTHER_ALLELE',
         'eaf' = 'AF_1000G',
         'beta' = 'BETA',
         'se' = 'SE',
         'pval' = 'P',
         'SNP' = 'RSID',
         'samplesize' ='N')

# Format data using TwoSampleMR
smokingExpMR <- TwoSampleMR::format_data(smokingExp, type='exposure')
message('Completed\n')

################################# CLUMPING ######################################
# LD clumping
message('LD Clumping...')
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
smokingExpClumped <- clump_data(smokingExpMR, clump_kb=500, clump_r2=0.01)
message('Completed\n')

################################ OUTCOME ####################################### 
# Increase timeout to allow downloads
options(timeout=10000000)
setwd('./Mendelian_Randomisation/SmkInit_Proteome_PPP/MR_SmkInit_Proteome/')

# Get deCODE GWAS URL
message('\nGetting URLs...')
deCODEurls <- readLines("./deCODE_GWAS/deCODE_URLs.txt")
message('Completed\n')

for (i in 1:10){ 

  # Set start time
  startTime<-Sys.time()
  
  # Get protein name from URL
  if (grepl('OLINK2', deCODEurls[i], fixed = TRUE)){protein <- gsub("_adjAgeSexPC_InvNorm.txt.gz", "", gsub("^.{90}_.{3}_.{6}_.{8}_", "", deCODEurls[i]))}
  if (grepl('OLINK_', deCODEurls[i], fixed = TRUE)){protein <- gsub("_adjAgeSexBatPC_InvNorm.txt.gz", "", gsub("^.{90}_.{3}_.{5}_.{8}_", "", deCODEurls[i]))}
  protein <- sub("_.*", "", protein)

  # Display start time for each protein    
  message(paste0('(', i, ') Performing MR for: ', protein))
  message(paste0('Start time: ', startTime))

  # Download file
  message('\nDownloading deCODE file...')
  print(Sys.time())
  download.file(deCODEurls[i], destfile=paste0('./deCODE_GWAS/', protein, '.txt.gz'), method='libcurl', quiet=FALSE)
  message('Completed\n')
  
  # Read file
  message('Reading deCODE file...')
  options(datatable.fread.datatable=FALSE)
  deCODE <- data.table::fread(paste0('./deCODE_GWAS/', protein, '.txt.gz'), header=TRUE)
  deCODE$Phenotype <- protein
  deCODE$Chrom <- as.factor(str_split(deCODE$Chrom, 'chr', simplify = TRUE)[,2])
  message('Completed\n')

  # Remove downloaded file after use
  message('Removing deCODE download...')
  if (file.exists(paste0('./deCODE_GWAS/', protein, '.txt.gz'))) {
  file.remove(paste0('./deCODE_GWAS/', protein, '.txt.gz'))}
  message('Completed\n')

  ######################### FORMATTING OUTCOME DATA ###############################
  # Rename columns
  message('Formatting deCODE file...')
  deCODE <- deCODE %>% 
      rename('BP' = 'Pos',
             'effect_allele' = 'A1',
             'other_allele' = 'A0',
             'eaf' = 'ImpFreqA1',
             'beta' = 'Beta',
             'se' = 'SE',
             'pval' = 'Pval',
             'SNP' = 'rsids',
             'samplesize'='N',
             'chr'='Chrom')

  # Format
  deCODEMR <- TwoSampleMR::format_data(dat=deCODE, snps=smokingExpClumped$SNP, type='outcome')
  message('Completed\n')
  
  ########################## IDENTIFYING PROXY SNPS ###############################
  # Get exposure SNPs present and absent in outcome
  message('\nIdentifying proxy SNPs...')
  SNPsPresent <- deCODEMR %>% filter(SNP %in% smokingExpClumped$SNP) 
  SNPsAbsent <- smokingExpClumped %>% filter(SNP %!in% deCODEMR$SNP)
  message('Completed')
  message('\nAbsent SNPs:')
  print(SNPsAbsent$SNP)

  # If SNPs are absent, identify proxy SNPs
  if(nrow(SNPsAbsent)>0){

    # Submit request to LDlink using modified function (MR_Functions.R)
    proxy <- LDproxy_batch(SNPsAbsent$SNP, 
                           pop = "EUR",             
                           r2d = "r2", 
                           token = 'e9f2ef6a36ed', 
                           append = TRUE,           
                           genome_build = "grch38",
                           code=protein) 
    
    # Munge proxies with proteome outcome data
    deCODEMR <- munge_proxies(paste0(protein, "_combined_query_snp_list_grch38.txt"), deCODEMR, SNPsPresent)
    message('Completed\n')
    message('Viewing munged file:')
    print(deCODEMR)
    
    # Remove proxy SNP file      
    message('Removing proxy SNPs file...')
    if (file.exists(paste0(protein, "_combined_query_snp_list_grch38.txt"))){
        file.remove(paste0(protein, "_combined_query_snp_list_grch38.txt"))}
    message('Completed')
    
  }

  # Harmonize exposure and outcome data
  harm <- harmonise_data(
    exposure_dat = smokingExpClumped, 
    outcome_dat = deCODEMR) 
  write.table(harm, file = paste0("./Results/Harm_SmkInit_", protein, ".txt"), sep = ",", quote = FALSE)
  message('Completed\n')

  ########################### IDENTIFYING OUTLIERS ################################
  # Format data for RadialMR
  radialdat <- harm %>% filter(mr_keep == TRUE) 
  radialdat <- radialdat %>% dat_to_RadialMR()
  rmr <- radialdat[[1]]
  
  # Run Radial MR for IVW
  message('\nRunning RadialMR...')
  radialIVW <- ivw_radial(rmr, alpha = 0.05/nrow(rmr)) 
  
  # Get outliers detected from RadialMR
  message('\nRadialMR Outliers:')
  outliers <- radialIVW$outliers
  print(outliers)
  message('Completed\n')

  # Create IVW and Egger plots (containing outliers if present)
  radialIVW <- ivw_radial(rmr, alpha = 0.05/nrow(rmr)) 
  radialEgger <- egger_radial(rmr, alpha = 0.05/nrow(rmr))
  radialIVWPlot <- plot_radial(radialIVW, radial_scale=F, show_outliers=T, scale_match=F)
  radialEggerPlot <- plot_radial(radialEgger, radial_scale=F, show_outliers=T, scale_match=F)
  RadialPlots <- cowplot::plot_grid(radialIVWPlot + theme_bw() + 
                                    theme(legend.position = 'bottom', text = ggplot2::element_text(size = 10, family='Franklin Gothic Book')) + 
                                    scale_colour_manual(breaks=c("Variant","Outlier","IVW"),values=c("#EE4800","#009900","#B23600")), 
                                    radialEggerPlot + theme_bw() + 
                                    theme(legend.position = 'bottom', text = ggplot2::element_text(size = 10, family='Franklin Gothic Book'))+ 
                                    scale_colour_manual(breaks=c("variant","Outlier","MR-Egger"), values=c("#A1D5D2","#009900","#42968F"), labels=c("Variant","Outlier","MR-Egger")), rel_widths = c(6,6))

  # Save plots
  ggsave(plot=RadialPlots, filename=paste0('./Results/Radial/', protein,'_PPP_Outliers.pdf'), width=7, height=5, units='in')
  
  # Set outliers removed to NO
  outliers_removed <- 'NO'

  ############################ REMOVING OUTLIERS ##################################
  # If outliers are present, remove from exposure GWAS and re-harmonise data
  if (class(outliers) == 'data.frame'){

    # Remove outliers from clumped data
    message('\nOutliers detected: removing outliers...')
    outliers <- outliers$SNP
    smokingExpClumpednoOutliers <- smokingExpClumped %>% filter(!SNP %in% outliers)
    outliers_removed <- 'YES'
    
    # Harmonize exposure and outcome data with outliers removed
    message('\n')
    harm <- harmonise_data(
      exposure_dat = smokingExpClumpednoOutliers, 
      outcome_dat = deCODEMR) 
    write.table(harm, file = paste0("./Results/Harm_SmkInit_", protein, "_OutliersRemoved.txt"), sep = ",", quote = FALSE)
    message('Completed\n')
  }

  ############################### RUNNING MR ###################################### 
  message('Performing MR...')
  message('Completed\n')
  res <- mr(harm, method_list = c("mr_egger_regression", "mr_simple_median", "mr_weighted_median", "mr_ivw_fe", "mr_ivw_mre"))
  write.table(res, file = paste0("./Results/Res_SmkInit_", protein, ".txt"), sep = ",", quote = FALSE)
  message('Results:')
  print(res)

  ######################### RUNNING ADDITIONAL TESTS ############################## 
  # Estimate heterogeneity
  message('\nPerforming Additional Tests...')
  message('\nHeterogeneity:')
  het <- mr_heterogeneity(harm) 
  print(het)
        
  # Testing for pleiotropy
  message('\nPleiotropy:')
  pleio <- mr_pleiotropy_test(harm)
  print(pleio)
        
  # Remove IVW method that is not applicable
  if(het$Q_pval[2]<0.05) res[4,] <- NA else res[5,] <- NA

  # Get individual F-stats
  message('\nIndividual F-statistics:')
  dat <- harm %>% filter(mr_keep == TRUE) 
  rows <- nrow(dat)
  total_F = 0
  num = 0
  for (row in 1:nrow(dat)){
    SNP = dat[row,]
    F_STAT = (SNP$beta.exposure**2)/(SNP$se.exposure**2)
    total_F <- total_F+F_STAT 
    num <- num+1
    message(num, ': ', F_STAT)}

  # Get mean F-stat
  message('\nMean F-statistic:')
  F_STAT = total_F/rows
  print(F_STAT)
     
  # Run Radial MR for IVW and MR Egger
  message('\nRunning RadialMR...')
  radialdat <- harm %>% filter(mr_keep == TRUE) 
  radialdat <- radialdat %>% dat_to_RadialMR()
  rmr <- radialdat[[1]]
  radialIVW <- ivw_radial(rmr, alpha = 0.05/nrow(rmr)) 
  radialEgger <- egger_radial(rmr, alpha = 0.05/nrow(rmr))
        
  # Get F-stat
  message('RadialMR F-stat:')
  F_STAT <- radialIVW$meanF
  print(F_STAT)
  
  # RadialMR
  message('\nRadialMR Outliers:')
  outliers <- radialIVW$outliers
  print(outliers)
  
  ############################## CREATING PLOTS ################################### 
  # Radial plots (should not contain outliers)
  message('\nCreating plots...')
  radialIVWPlot <- plot_radial(radialIVW, radial_scale=F, show_outliers=T, scale_match=F)
  radialEggerPlot <- plot_radial(radialEgger, radial_scale=F, show_outliers=T, scale_match=F)
    
  # Scatter plot
  Scatter <- mr_scatter_plot(res, harm) 
  ScatterPlot <- Scatter[[1]] + theme_bw() + 
    guides(color=guide_legend(ncol = 1)) + 
    theme(text = element_text(size = 12, family='Franklin Gothic Book'))
    
  # Leave-one-out analysis 
  LOO <- mr_leaveoneout(harm, method = mr_ivw_fe) 
  LOOPlot <- mr_leaveoneout_plot(LOO)

  # Single SNP analysis 
  SingleSNP <- mr_singlesnp(harm, all_method = c("mr_ivw_fe", "mr_egger_regression", "mr_weighted_median")) %>% as_tibble()
  SingleSNPPlot <- mr_forest_plot(SingleSNP) 
    
  # Funnel plot
  Funnel <- mr_funnel_plot(SingleSNP)
  FunnelPlot <- Funnel[[1]] + theme_classic() + 
    guides(color=guide_legend(ncol=1)) + 
  theme(text = ggplot2::element_text(size = 12, family='Franklin Gothic Book'))
    
  message('Completed\n')

  ############################## SAVING RESULTS ################################### 
  message('\nSaving results...')
  summary <- data.frame('outcome'=res$outcome[1], 
                        'exposure'=res$exposure[1], 
                        'nSNP'=res$n[1], 
                        'MR_Egger_Beta'=res$b[1], 
                        'Simple_Median_Beta'=res$b[2], 
                        'Weighted_Median_Beta'=res$b[3], 
                        'IVW_FE_Beta'=res$b[4], 
                        'IVW_MRE_Beta'=res$b[5],
                        'MR_Egger_SE'=res$se[1], 
                        'Simple_Median_SE'=res$se[2], 
                        'Weighted_Median_SE'=res$se[3], 
                        'IVW_FE_SE'=res$se[4], 
                        'IVW_MRE_SE'=res$se[5],
                        'MR_Egger_pval'=res$pval[1], 
                        'Simple_Median_pval'=res$pval[2], 
                        'Weighted_Median_pval'=res$pval[3], 
                        'IVW_FE_pval'=res$pval[4], 
                        'IVW_MRE_pval'=res$pval[5],
                        'Het_pval'=het$Q_pval[2], 
                        'Pleio_pval'=pleio$pval[1], 
                        'Outliers'=paste(unlist(outliers), collapse='|'), 
                        'Mean_F_statistic'=F_STAT[1],
                        'SNP_pval'=SNPpval,
                        'Outliers_Removed?'=outliers_removed)

  write.table(summary, file = "./Results/Results.txt", sep = ",", append = TRUE, quote = FALSE, 
              col.names = !file.exists("./Results/Results.txt"), row.names = FALSE) 

  ggsave(plot=ScatterPlot, filename=paste0('./Results/Scatter/SmkInit_', protein,'.pdf'), width=10.01, height=6.29, units='in')
  ggsave(plot=LOOPlot[[1]], filename=paste0('./Results/LOO/SmkInit_', protein,'.pdf'), width=7.01, height=20.29, units='in')
  ggsave(plot=SingleSNPPlot[[1]], filename=paste0('./Results/SingleSNP/SmkInit_', protein,'.pdf'), width=7.01, height=20.29, units='in')
  ggsave(plot=FunnelPlot, filename=paste0('./Results/Funnel/SmkInit_', protein,'.pdf'), width=7.01, height=4.29, units='in')
    
  RadialPlots <- cowplot::plot_grid(radialIVWPlot + theme_bw() + 
                                    theme(legend.position = 'bottom', text = ggplot2::element_text(size = 10, family='Franklin Gothic Book')) + 
                                    scale_colour_manual(breaks=c("Variant","Outlier","IVW"),values=c("#EE4800","#009900","#B23600")), 
                                    radialEggerPlot + theme_bw() + 
                                    theme(legend.position = 'bottom', text = ggplot2::element_text(size = 10, family='Franklin Gothic Book'))+ 
                                    scale_colour_manual(breaks=c("variant","Outlier","MR-Egger"), values=c("#A1D5D2","#009900","#42968F"), labels=c("Variant","Outlier","MR-Egger")), rel_widths = c(6,6))

  ggsave(plot=RadialPlots, filename=paste0('./Results/Radial/SmkInit_', protein,'.pdf'), width=7.01, height=4.29, units='in')
  message('Completed\n')
  
  message(paste0('Completed: ', protein))
  endTime <- Sys.time()
  message(paste0('Time taken: ', endTime-startTime, ' minutes'))

}  

################################## FINISHED ####################################### 
