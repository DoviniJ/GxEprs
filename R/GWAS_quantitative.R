#' GWAS_quantitative function
#' This function performs GWAS using plink2 and outputs the GWAS summary statistics with additive SNP effects. Users may save the output in a user-specified file (see example).
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param Qphe_discovery Name (with file extension) of the phenotype file containing family ID, individual ID and phenotype of the discovery dataset as columns, without heading
#' @param Qcov_discovery Name (with file extension) of the covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the discovery dataset as columns, without heading
#' @param thread Number of threads used
#' @keywords gwas
#' @export 
#' @importFrom stats binomial fitted.values glm lm 
#' @importFrom utils read.table 
#' @return This function will perform GWAS and output
#' \item{Q_out.trd.sum}{GWAS summary statistics with additive SNP effects}
#' @examples \dontrun{
#' x <- GWAS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery, 
#' thread = 20)
#' sink("Q_out.trd.sum") #to create a file in the working directory
#' write.table(x[c("ID", "A1", "BETA")], sep = " ", 
#' row.names = FALSE, quote = FALSE) #to write the output
#' sink() #to save the output
#' head(x) #to obtain the head of GWAS summary statistics of additive SNP effects
#' x$CHROM #to extract the chromosome number 
#' x$POS #to extract the base pair position 
#' x$ID #to extract the SNP ID 
#' x$REF #to extract the reference allele 
#' x$ALT #to extract the alternate allele 
#' x$A1 #to extract the minor allele 
#' x$OBS_CT #to extract the number of allele observations 
#' x$BETA #to extract the SNP effects
#' x$SE #to extract the standard errors of the SNP effects
#' x$T_STAT #to extract the test statistics 
#' x$P #to extract the p values 
#' }
GWAS_quantitative <- function(plink_path, b_file, Qphe_discovery, Qcov_discovery, thread = 20){  
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
  log_file <- runPLINK(paste0(" --bfile ", b_file, 
                  " --linear interaction --pheno ", 
                  Qphe_discovery, 
                  " --covar ", Qcov_discovery, 
                  " --parameters ", param_vec, 
                  " --allow-no-sex --threads ", 
                  thread,
                  " --out ", tempdir(),"/Q_gwas"))
  plink_output <- read.table(paste0(tempdir(), "/Q_gwas.PHENO1.glm.linear"), header = FALSE)
  filtered_output <- plink_output[(plink_output$V7=="ADD"),]
  colnames(filtered_output) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")
  Q_out.trd.sum <- filtered_output[c("CHROM", "POS", "ID", "REF", "ALT", "A1", "OBS_CT", "BETA", "SE", "T_STAT", "P")]
  rownames(Q_out.trd.sum) <- NULL
  return(Q_out.trd.sum)
}
