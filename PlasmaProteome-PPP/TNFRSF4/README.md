# TNFRSF4 → PPP Mendelian Randomisation: MHC Region & Pleiotropy
Further analysis of the apparent affect of TNFRSF4 on PPP. 

* Following my first bulk analysis, which identified TNFRSF4 as a potential mediator of SmkInit → PPP, I removed 4 SNPs from the original analysis that are located in the MHC region.
  * This region has complex linkage disequilibrium patterns and can complicate MR analyses (in particular, after the clumping step of MR some SNPs that are in LD with each other remain).
  * The 4 SNPs were removed by re-clumping the 5 SNPs located in the MHC region (`MHC-SNPs-TNFRSF4.png`) - these 5 SNPs were in LD but missed the first round of clumping as they are quite distant from each other.
  * The TNFRSF4 → PPP relationship persisted using the remaining 14 SNPs (`TNFRSF4-MHCrm-PPP-Results.csv`).

* I then looked at sources of pleiotropy that might arise from the SNPs used as TNFRSF4 IVs: other protein levels may be influencing this finding.
  * `ProteomePPPIVsOverlap.R` looks at which TNFRSF4 IVs are pleiotropic: which IVs are also IVs for other proteins measured in Step 2 of this analysis.
  * The result is displayed in `Pleiotropy-Summary-Plot.png`.
  * `Pleiotropic-SNPs-TNFRSF4.png` is an MR plot including the 14 SNPs obtained after MHC-SNP removal, which highlights the 3 most pleiotropic SNPs.
