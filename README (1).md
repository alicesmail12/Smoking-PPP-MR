# SmkInit → Plasma Proteome Mendelian Randomisation 
SmkInitPlasmaProteome.R performs separate MR analyses to estimate the causal effect of smoking initiation on different plasma protein levels.

### Input files
* Exposure GWAS: Smoking Initiation (GSCANSmkInit2022EUR.txt)
  * Saunders et al. (2022) *Genetic diversity fuels gene discovery for tobacco and alcohol use*.
    * Phenotype: ever or never being a regular smoker.
    * Population: 2,669,029 European individuals from the US, Europe & Australia.

- Outcome GWAS: deCODE Plasma Proteome URLs (PlasmaProteomeURLs.txt)
  - Eldjarn et al. (2023) *Large-scale plasma proteomics comparisons through genetics and disease associations*.
    - Phenotype: plasma protein abundance measured using Olink.
    - Population: 46,218 British or Irish individuals from the UK BioBank.
  - URL list obtained from www.decode.com/summarydata/

### Steps
1. SmkInit GWAS filtered & formatted for MR.
	- SNPs with p<1e-8 identified & MHC region removed.
	- Clumped at 500kb and r<sup>2</sup> of 0.01.
2. For each plasma protein GWAS:
	- File downloaded and read.
	- GWAS filtered using SmkInit IVs and formatted for MR.
	- SNPs harmonised prior to performing MR.
	- Additional tests performed & plots created.
3. Outputs saved.

* If SNPs in the SmkInit IVs are absent from the outcome GWAS:
	* Proxies can be identified using the *LDproxy_batch* function (https://www.rdocumentation.org/packages/LDlinkR/versions/1.4.0/topics/LDproxy_batch) and munged using the *munge_proxies* function from https://andrewslabucsf.github.io/MR-tutorial/scripts/mr_harmonization.html. These are saved in MR_Functions.R.
