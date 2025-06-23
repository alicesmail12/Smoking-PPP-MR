# TNFRSF4 → PPP Mendelian Randomisation: MHC Region & Pleiotropy
Further analysis of the apparent affect of TNFRSF4 on PPP. 

### MHC Region
* Following my first bulk analysis, which identified TNFRSF4 as a potential mediator of SmkInit → PPP, I removed 4 MHC SNPs from the original analysis as these may have been driving the result.
  * The MHC region has complex linkage disequilibrium patterns and can complicate MR analyses (in particular, after the clumping step of MR some SNPs that are in LD with each other remain).
  * The 4 SNPs were removed by re-clumping the 5 SNPs located in the MHC region (`MHC-SNPs-TNFRSF4.png`) - these 5 SNPs are in LD but passed the first round of clumping as they are quite distant from each other.
  * The TNFRSF4 → PPP relationship persisted using the remaining 14 SNPs (`TNFRSF4-MHCrm-PPP-Results.csv` & `Pleiotropic-SNPs-TNFRSF4.png`).
  * The 14 SNPs used here are detailed in `TNFRSF4-MHCrm-PPP-Harmonised.csv` (GWAS information removed as this is not my own data).

### Pleiotropy
* I then looked at sources of pleiotropy that might arise from the SNPs used as TNFRSF4 IVs, as other protein levels may be influencing this finding.
  * `TNFRSF4-PPP-IVs-Overlap.R` looks at which TNFRSF4 IVs are pleiotropic: which IVs are also IVs for other proteins measured in Step 2 of this analysis.
  * The result is displayed in `Pleiotropy-Summary-Plot.png`, where asterisks denote significant associations. SNPs rs7525284 and rs78037977 are located in proximity of the TNFRSF4 gene region. The remaining SNPs are trans-acting. 
  * `Pleiotropic-SNPs-TNFRSF4.png` is an MR plot including the 14 SNPs obtained after MHC-SNP removal. The 3 most pleiotropic SNPs are highlighted in light blue.
