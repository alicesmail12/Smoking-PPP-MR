# Plasma Proteome → PPP Mendelian Randomisation 
PlasmaProteomePPP.R performs separate MR analyses to estimate the causal effect of different plasma protein levels on PPP.

### Input files
- Exposure GWAS: UK BioBank Plasma Proteome URLs (PlasmaProteomeURLs.txt)
  - Eldjarn et al. (2023) *Large-scale plasma proteomics comparisons through genetics and disease associations*.
    - Phenotype: plasma protein abundance measured using Olink.
    - Population: 46,218 British or Irish individuals from the UK BioBank.
  - URL list obtained from www.decode.com/summarydata/

* Outcome GWAS: PPP (PPPGWAS.txt)
  * Meta-analysis of 3 PPP GWAS datasets:
    * UK: 288 cases & 7,321 controls
    * HUNT: 225 cases & 64,050 controls
    * FinnGenn: 969 cases & 330,975 controls
  * Accession ID: GCST90624455

### Steps
1. PPP GWAS filtered & formatted for MR.
2. For each plasma protein GWAS:
	- File downloaded and read.
  	- SNPs with p<1e-8 identified & MHC region removed.
	- Clumped at 500kb and r<sup>2</sup> of 0.01.
	- SNPs harmonised prior to performing MR.
	- Additional tests performed & plots created.
3. Outputs saved.

* SNPs in the Plasma Proteome IVs that were absent from the outcome GWAS were identified using the *LDproxy_batch* function (https://www.rdocumentation.org/packages/LDlinkR/versions/1.4.0/topics/LDproxy_batch) and munged using the *munge_proxies* function from https://andrewslabucsf.github.io/MR-tutorial/scripts/mr_harmonization.html. These are saved in `MR_Functions.R`.
