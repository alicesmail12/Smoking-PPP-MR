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

################################# OUTCOME #######################################
message('\nReading PPP file...')
options(datatable.fread.datatable=FALSE)
PPP <- data.table::fread('./Mendelian_Randomisation/GWAS_PPP/PPPGWAS.txt', header=TRUE)

# Set phenotype and size
PPP$phenotype <- 'PPP'
PPP$N <- 403506
message('Completed\n')

# Reformat outcome data
message('Formatting PPP file...')
PPPMR <- TwoSampleMR::format_data(
  snps = NULL,
  type='outcome',
  PPP,
  header=TRUE,
  phenotype_col = 'phenotype',
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  chr='CHR',
  samplesize='N',
  pos='BP',
  eaf_col = 'Freq1')
message('Completed\n')

################################ EXPOSURE ####################################### 
# Increase timeout to allow downloads
options(timeout=10000000)
setwd('./Mendelian_Randomisation/SmkInit_Proteome_PPP/MR_Proteome_PPP/')

# Get deCODE GWAS URL
message('\nGetting URLs...')
deCODEurls <- readLines("./deCODE_GWAS/deCODE_URLs.txt")
message('Completed\n')

for (i in 1:10){

  # Set start time
  startTime <- Sys.time()

  # Get protein name from URL
  if (grepl('OLINK2', deCODEurls[i], fixed = TRUE)){protein <- gsub("_adjAgeSexPC_InvNorm.txt.gz", "", gsub("^.{90}_.{3}_.{6}_.{8}_", "", deCODEurls[i]))}
  if (grepl('OLINK_', deCODEurls[i], fixed = TRUE)){protein <- gsub("_adjAgeSexBatPC_InvNorm.txt.gz", "", gsub("^.{90}_.{3}_.{5}_.{8}_", "", deCODEurls[i]))}
  protein <- sub("_.*", "", protein)

  # Display start time for each protein
  message(paste0('(', i, ') Performing MR for: ', protein))
  message(paste0('Start time: ', startTime))

  # Download file
  message('\nDownloading deCODE file...')
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
  if (file.exists(paste0('./deCODE_GWAS/', protein, '.txt.gz'))) {file.remove(paste0('./deCODE_GWAS/', protein, '.txt.gz'))}
  message('Completed\n')

  ######################### IDENTIFYING EXPOSURE IVS ##############################
  # Remove MHC region
  # GRCh38.p14 Location: chr6:28,510,120-33,480,577
  # https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC
  
  message('Removing MHC SNPs...')
  MHC <- deCODE %>%
    arrange(Chrom, Pos) %>%
    filter(Chrom == 6) %>%
    filter(Pos >= 28510120 & Pos <= 33480577)
  MHCSNP <- MHC %>% dplyr::select(rsids)
  deCODE <- deCODE %>% dplyr::filter(!rsids %in% unlist(MHCSNP))
  message('Completed\n')
  
  # Get significant SNPs
  message('Filtering SNPs using pval threshold...')
  pval <- 1e-8
  deCODE <- deCODE %>% filter(Pval <= pval)
  message('Completed\n')

  ######################### FORMATTING EXPOSURE DATA ##############################
  # Rename columns
  message('Formatting data...')
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
  deCODEMR <- TwoSampleMR::format_data(deCODE)
  message('Completed\n')

  ################################# CLUMPING ######################################
  # LD clumping
  message('LD Clumping...')
  options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
  deCODEMRClumped <- clump_data(deCODEMR, clump_kb=500, clump_r2=0.01)
  message('Completed\n')

  ############################## HARMONISATION ####################################
  # Harmonize exposure and outcome data
  message('Harmonising data...')
  harm <- harmonise_data(exposure_dat = deCODEMRClumped, 
                         outcome_dat = PPPMR) 

  # Save output to generate plots externally
  write.table(harm, file = paste0("./Results/Harm_", protein, "_PPP.txt"), sep = ",", quote = FALSE)
  message('Completed\n')

  ########################## IDENTIFYING PROXY SNPS ###############################
  # Get exposure SNPs present and absent in outcome GWAS
  message('Identifying proxy SNPs...')
  SNPsPresent <- PPPMR %>% filter(SNP %in% deCODEMRClumped$SNP) 
  SNPsAbsent <- deCODEMRClumped %>% filter(SNP %!in% PPPMR$SNP) %>% filter(mr_keep.exposure == 'TRUE')
  message('Absent SNPs:')
  print(SNPsAbsent)
  message('Completed\n')

  # If SNPs are absent, identify proxy SNPs
  if(nrow(SNPsAbsent)>0){

    # Submit request to LDlink using modified function (MR_Functions.R)
    proxy <- LDproxy_batch(SNPsAbsent$SNP, 
                           pop = "EUR",             
                           r2d = "r2", 
                           token = 'e9f2ef6a36ed', 
                           append = TRUE,           
                           genome_build = "grch37",
                           code=protein) 
    
    # Munge proxies with PPP outcome data
    PPPMRProxies <- munge_proxies(paste0(protein, "_combined_query_snp_list_grch37.txt"), PPPMR, SNPsPresent)
    message('Completed\n')

    # Remove proxy SNP file
    message('\nRemoving proxy SNPs file...')
    if (file.exists(paste0(protein, "_combined_query_snp_list_grch37.txt"))){
      file.remove(paste0(protein, "_combined_query_snp_list_grch37.txt"))}
    message('Completed\n')
    
    # Harmonize exposure and outcome data including proxy SNPs
    message('Harmonising data...')
    harm <- harmonise_data(exposure_dat = deCODEMRClumped, 
                           outcome_dat = PPPMRProxies) 

    # Write new harmonised file including proxies
    write.table(harm, file = paste0("./Results/Harm_", protein, "_PPP_Proxies.txt"), sep = ",", quote = FALSE)
    message('Completed\n')
    
  }

  ########################### IDENTIFYING OUTLIERS ################################
  # Format data for RadialMR
  message('\nRunning RadialMR...')
  radialdat <- harm %>% filter(mr_keep == TRUE) 
  radialdat <- radialdat %>% dat_to_RadialMR()
  rmr <- radialdat[[1]]
  
  # Run Radial MR for IVW
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
  
  # Set outliers removed value to NO
  outliers_removed <- 'NO'

  ############################ REMOVING OUTLIERS ##################################
  # If outliers are present, remove from exposure GWAS and re-harmonise data
  if (class(outliers) == 'data.frame'){

    # Remove outliers from clumped data
    message('\nOutliers detected: removing outliers...')
    outliers <- outliers$SNP
    deCODEMRClumpednoOutliers <- deCODEMRClumped %>% filter(!SNP %in% outliers)
    outliers_removed <- 'YES'

    # Harmonize exposure and outcome data with outliers removed
    message('\n')
    harm <- harmonise_data(exposure_dat = deCODEMRClumpednoOutliers, outcome_dat = PPPMR) 
    write.table(harm, file = paste0("./Results/Harm_", protein, "_PPP_OutliersRemoved.txt"), sep = ",", quote = FALSE)
    message('Completed\n')
  }

  ############################### RUNNING MR ###################################### 
  # Perform MR
  message('Performing MR...')
  message('Completed\n')
  res <- mr(harm, method_list = c("mr_egger_regression", "mr_simple_median", "mr_weighted_median", "mr_ivw_fe", "mr_ivw_mre"))
  message('Results:')
  print(res)
  
  # Save results
  write.table(res, file = paste0("./Results/Res_", protein, "_PPP.txt"), sep = ",", quote = FALSE)

  ######################### RUNNING ADDITIONAL TESTS ############################## 
  # Estimate heterogeneity
  message('\nPerforming Additional Tests...')
  message('\nHeterogeneity:')
  het <- mr_heterogeneity(harm) 
  print(het)
    
  # Estimate heterogeneity for MR-Egger
  hetRuckers <- mr_rucker(harm)
  message('\nMR-Egger Heterogeneity:')
  print(hetRuckers[[1]]$Q) 

  # Write Cochran's vs Rucker's Q to file
  write.table(hetRuckers[[1]]$Q, file = paste0("./Results/Ruckers_", protein, "_PPP.txt"), sep = ",", quote = FALSE)
  
  # Single SNP and I2
  SingleSNP <- mr_singlesnp(harm)
  SingleSNPbeta <- SingleSNP$b
  SingleSNPse <- SingleSNP$se
  I2 <- Isq(SingleSNPbeta, SingleSNPse)
    
  # Testing for pleiotropy
  message('\nPleiotropy:')
  pleio <- mr_pleiotropy_test(harm)
  print(pleio)
    
  # Remove IVW method that is not applicable (based on whether significant heterogeneity is detected)
  if(het$Q_pval[2]<0.05) res[4,] <- NA else res[5,] <- NA
 
  # Get F-stat for each SNP
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
    
  # Get F-stat from RadialMR (should match manually calculated value)
  message('\n(RadialMR F-stat:')
  F_STAT <- radialIVW$meanF
  print(F_STAT)
    
  # RadialMR
  message('\nRadialMR Outliers:')
  radialIVW$outliers
  print(radialIVW$outliers)

  ############################## CREATING PLOTS ################################### 
  # Radial plots (should not contain outliers)
  message('Creating plots...')
  radialIVWPlot <- plot_radial(radialIVW, radial_scale=F, show_outliers=T, scale_match=F)
  radialEggerPlot <- plot_radial(radialEgger, radial_scale=F, show_outliers=T, scale_match=F)
  
  RadialPlots <- cowplot::plot_grid(radialIVWPlot + theme_bw() + 
                                    theme(legend.position = 'bottom', text = ggplot2::element_text(size = 10, family='Franklin Gothic Book')) + 
                                    scale_colour_manual(breaks=c("Variant","Outlier","IVW"),values=c("#EE4800","#009900","#B23600")), 
                                    radialEggerPlot + theme_bw() + 
                                    theme(legend.position = 'bottom', text = ggplot2::element_text(size = 10, family='Franklin Gothic Book'))+ 
                                    scale_colour_manual(breaks=c("variant","Outlier","MR-Egger"), values=c("#A1D5D2","#009900","#42968F"), labels=c("Variant","Outlier","MR-Egger")), rel_widths = c(6,6))
  
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
                        'Het_Q'=het$Q[2], 
                        'Het_df'=het$Q_df[2], 
                        'Het_pval'=het$Q_pval[2], 
                        'I2'=I2,
                        'Pleio_pval'=pleio$pval[1], 
                        'Outliers'=paste(unlist(outliers), collapse='|'), 
                        'Mean_F_statistic'=F_STAT[1],
                        'SNP_pval'=pval,
                        'Outliers_Removed?'=outliers_removed)
  
  write.table(summary, file = "./Results/Results.txt", sep = ",", append = TRUE, quote = FALSE, 
              col.names = !file.exists("./Results/Results.txt"), row.names = FALSE) 
  
  ggsave(plot=ScatterPlot, filename=paste0('./Results/Scatter/', protein,'_PPP.pdf'), width=10, height=7, units='in')
  ggsave(plot=LOOPlot[[1]], filename=paste0('./Results/LOO/', protein,'_PPP.pdf'), width=7, height=7, units='in')
  ggsave(plot=SingleSNPPlot[[1]], filename=paste0('./Results/SingleSNP/', protein,'_PPP.pdf'), width=7, height=7, units='in')
  ggsave(plot=FunnelPlot, filename=paste0('./Results/Funnel/', protein,'_PPP.pdf'), width=7, height=5, units='in')
  ggsave(plot=RadialPlots, filename=paste0('./Results/Radial/', protein,'_PPP.pdf'), width=7, height=5, units='in')
  message('Completed')
  
  message(paste0('\nCompleted: ', protein, '\n'))

  endTime <- Sys.time()
  message(paste0('Time taken: ', endTime-startTime, ' minutes'))

}

################################## FINISHED ####################################### 
