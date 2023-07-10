#' GWAS_quantitative function
#' This function performs GWAS using plink2 and outputs the GWAS summary statistics file with additive SNP effects named Q_out.trd.sum
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param Qphe_discovery Name (with file extension) of the phenotype file containing family ID, individual ID and phenotype of the discovery dataset as columns, without heading
#' @param Qcov_discovery Name (with file extension) of the covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the discovery dataset as columns, without heading
#' @param thread Number of threads used
#' @param summary_output Name of the SNP effects of the GWAS summary statistics file specified by the user
#' @keywords gwas
#' @export 
#' @importFrom stats binomial fitted.values glm lm
#' @importFrom utils read.table write.table
#' @return This function will perform GWAS and output
#' \item{Q_out.trd.sum}{GWAS summary statistics file with additive SNP effects}
#' @examples \dontrun{
#' x <- GWAS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery, 
#'                        thread = 20, summary_output = "Q_out.trd.sum")
#' head(x) #to read the head of all columns in Q_out.trd.sum file
#' x$V1 #to extract the chromosome number (CHROM)
#' x$V2 #to extract the base pair position (POS)
#' x$V3 #to extract the SNP ID (ID)
#' x$V4 #to extract the reference allele (REF)
#' x$V5 #to extract the alternate allele (ALT)
#' x$V6 #to extract the minor allele (A1)
#' x$V7 #to extract the type of test performed (TEST)
#' x$V8 #to extract the nmber of allele observations (OBS_CT)
#' x$V9 #to extract the SNP effect (BETA)
#' x$V10 #to extract the standard error of each SNP effect (SE)
#' x$V11 #to extract the test statistic (T_STAT)
#' x$V12 #to extract the p value (P)
#' x$V13 #to extract the error code (ERRCODE)
#' }
GWAS_quantitative <- function(plink_path, b_file, Qphe_discovery, Qcov_discovery, thread = 20, summary_output = "Q_out.trd.sum"){  
  cov_file <- read.table(Qcov_discovery)
  n_confounders = ncol(cov_file) - 4
  if(n_confounders > 0){
    parameters <- c(1, (1:n_confounders)+3)
  }
  else{
    parameters <- 1
  }
  param_vec <- paste0(parameters, collapse = ", ")
  runPLINK <- function(PLINKoptions = "") system(paste(plink_path, PLINKoptions))
  runPLINK(paste0(" --bfile ", b_file, 
                  " --linear interaction --pheno ", 
                  Qphe_discovery, 
                  " --covar ", Qcov_discovery, 
                  " --parameters ", param_vec, 
                  " --allow-no-sex --threads ", 
                  thread,
                  " --out Q_gwas"))
  plink_output <- read.table("Q_gwas.PHENO1.glm.linear", header = FALSE)
  filtered_output <- plink_output[(plink_output$V7=="ADD"),]
  sink(summary_output)
  write.table(filtered_output, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
  return(filtered_output)
}
