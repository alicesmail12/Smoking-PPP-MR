############################## GET LD PROXIES ################################### 
#' This function is modified from the LDlinkR package (https://rdrr.io/cran/LDlinkR/src/R/LDproxy_batch.R)
#' Reference: Myers TA, Chanock SJ, Machiela MJ. 
#' LDlinkR: An R Package for Rapidly Calculating Linkage Disequilibrium Statistics in Diverse Populations. 
#' Front Genet. 2020 Feb 28; 11:157. 
#' doi: 10.3389/fgene.2020.00157. 
#' PMID: 32180801; PMCID: PMC7059597.

LDproxy_batch <- function(snp,
                          pop="CEU",
                          r2d="r2",
                          token=NULL,
                          append = FALSE,
                          genome_build = "grch37",
                          api_root="https://ldlink.nih.gov/LDlinkRest",
                          code=code) {
  
  snp <- as.data.frame(snp)
  
  if(append == FALSE) {
    for (i in 1:nrow(snp)) {
      myfile <- paste(snp[i,], "_", genome_build, ".txt", sep="")
      cat("\nSubmitting request for query variant ", snp[i,],".", sep = "")
      cat("\nChecking status of server...")
      df_proxy <- LDproxy(snp = snp[i,],
                          pop = pop,
                          r2d = r2d,
                          token = token,
                          genome_build = genome_build,
                          api_root = api_root)
      if(!(grepl("error", df_proxy[1,1])))
      {
        write.table(df_proxy, file = myfile,
                    append = FALSE,
                    quote = FALSE,
                    row.names = TRUE,
                    sep = "\t")
        file_path <- getwd()
        cat("File for query variant ", snp[i,], " saved to:\n ", file_path, "/", myfile,"\n", sep="")
      }
    }
  } else if (append == TRUE) {
    myfile <- paste0(code, "_combined_query_snp_list_", genome_build, ".txt")
    for (i in 1:nrow(snp)) {
      cat("\nSubmitting request for query variant ", snp[i,],".", sep = "")
      cat("\nChecking status of server...")
      df_proxy <- LDproxy(snp = snp[i,],
                          pop = pop,
                          r2d = r2d,
                          token = token,
                          genome_build = genome_build,
                          api_root = api_root)
      if(!(grepl("error", df_proxy[1,1])))
      {
        # add new column, query_snp
        df_proxy["query_snp"] <- rep(snp[i,], nrow(df_proxy))
        # rearrange by column index
        df_proxy <- df_proxy[, colnames(df_proxy)[c(ncol(df_proxy), 1:(ncol(df_proxy)-1))]]
        suppressWarnings(
          write.table(df_proxy,
                      file = myfile,
                      append = TRUE,
                      quote = FALSE,
                      row.names = TRUE,
                      col.names = !file.exists(myfile),
                      sep = "\t")
        )
      }
    }
    file_path <- getwd()
    cat("Combined file for all query variants saved to:\n", file_path, "/", myfile,"\n", sep="")
  }
}

############################# MUNGE LD PROXIES ################################## 
#' This function is defined at https://andrewslabucsf.github.io/MR-tutorial/scripts/mr_harmonization.html
#' Reference: Andrews SJ.
#' Mendelian Randomisation Tutorial.
#' 2023 Feb 2.

# Adds LD proxy variants to outcome dataset
munge_proxies <- function(LDLink_file, outcome, outcome_clump){
  LDLink_file_path <- LDLink_file
  proxy_snps <- read_tsv(LDLink_file_path, skip = 1, col_names = F) %>%
    rename(id = X1, func = X2, proxy_snp = X3, coord = X4, alleles = X5, maf = X6, 
           distance = X7, dprime = X8, rsq = X9, correlated_alleles = X10, FORGEdb = X11, RegulomeDB = X12) %>%
    separate(coord, c('chr', 'pos'), sep = ":") %>%
    mutate(snp = ifelse(id == 1, proxy_snp, NA), 
           chr = str_replace(chr, 'chr', ""), 
           chr = as.numeric(chr), 
           pos = as.numeric(pos)) %>%
    fill(snp, .direction = 'down') %>%
    relocate(snp, .before = proxy_snp) %>%
    dplyr::select(-id, -func, -FORGEdb, -RegulomeDB) %>%
    filter(rsq >= 0.8)
}
  
