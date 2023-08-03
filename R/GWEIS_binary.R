#' GWEIS_binary function
#' This function performs GWEIS using plink2 and outputs the GWEIS summary statistics with additive SNP effects and interaction SNP effects. Users may save the outputs in separate user-specified files (see examples).
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param Bphe_discovery Phenotype file containing family ID, individual ID and phenotype of the discovery dataset as columns, without heading
#' @param Bcov_discovery Covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the discovery dataset as columns, without heading
#' @param thread Number of threads used
#' @keywords gwies interaction gxe
#' @export 
#' @importFrom stats binomial fitted.values glm lm
#' @importFrom utils read.table 
#' @return This function will perform GWEIS and output
#' \item{B_out.sum}{GWEIS summary statistics with additive and interaction SNP effects}
#' @examples \dontrun{ 
#' x <- GWEIS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery, 
#' thread = 20)
#' sink("B_out.add.sum") #to create a file in the working directory
#' write.table(x[c("ID", "A1", "ADD_OR")], sep = " ", 
#' row.names = FALSE, quote = FALSE) #to write the output
#' sink() #to save the output
#' sink("B_out.gxe.sum") #to create a file in the working directory
#' write.table(x[c("ID", "A1", "INTERACTION_OR")], sep = " ", 
#' row.names = FALSE, quote = FALSE) #to write the output
#' sink() #to save the output
#' head(x) #to extract the head of all columns in GWEIS summary statistics of 
#' additive and interaction SNP effects 
#' x$CHROM #to extract the chromosome number 
#' x$POS #to extract the base pair position
#' x$ID #to extract the SNP ID
#' x$REF #to extract the reference allele
#' x$ALT #to extract the alternate allele 
#' x$A1 #to extract the minor allele
#' x$OBS_CT #to extract the number of allele observations 
#' x$ADD_OR #to extract the odds ratios of additive SNP effects
#' x$ADD_LOG_OR_SE #to extract the standard errors of log odds of 
#' additive SNP effects
#' x$ADD_Z_STAT #to extract the test statistics of additive SNP effects
#' x$ADD_P #to extract the p values of additive SNP effects
#' x$INTERACTION_OR #to extract the odds ratios of the SNP effects of 
#' interaction SNP effects
#' x$INTERACTION_LOG_OR_SE #to extract the standard errors of log odds 
#' of interaction SNP effects
#' x$INTERACTION_Z_STAT #to extract the test statistics of interaction 
#' SNP effects
#' x$INTERACTION_P #to extract the p values of interaction SNP effects
#' }
GWEIS_binary <- function(plink_path, b_file, Bphe_discovery, Bcov_discovery, thread = 20){
  cov_file <- read.table(Bcov_discovery)
  n_confounders = ncol(cov_file) - 4
  if(n_confounders > 0){
    parameters <- c(1, 2, 3, (1:(n_confounders+1))+3)
  }
  else{
    parameters <- c(1, 2, 3, 4)
  }
  param_vec <- paste0(parameters, collapse = ", ")
  runPLINK <- function(PLINKoptions = "") system(paste(plink_path, PLINKoptions))
  log_file <- runPLINK(paste0(" --bfile ", b_file,
                " --glm interaction --pheno ", 
                Bphe_discovery, 
                " --covar ", Bcov_discovery, 
                " --parameters ", param_vec, 
                " --allow-no-sex --threads ", 
                thread,
                " --out ", tempdir(),"/B_gweis"))
  plink_output <- read.table(paste0(tempdir(), "/B_gweis.PHENO1.glm.logistic.hybrid"), header = FALSE)
  filtered_output <- as.data.frame(plink_output[(plink_output$V8=="ADD"),])
  filtered_output$V10 = log(filtered_output$V10)
  colnames(filtered_output) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "FIRTH", "TEST", "OBS_CT", "ADD_OR", "ADD_LOG_OR_SE", "ADD_Z_STAT", "ADD_P", "ADD_ERRCODE")
  filtered_output2 <- plink_output[(plink_output$V8=="ADDxCOVAR1"),]
  filtered_output2$V10 <- log(filtered_output2$V10)
  colnames(filtered_output2) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "FIRTH", "TEST", "OBS_CT", "INTERACTION_OR", "INTERACTION_LOG_OR_SE", "INTERACTION_Z_STAT", "INTERACTION_P", "INTERACTION_ERRCODE")
  out <- cbind(filtered_output[c("CHROM", "POS", "ID", "REF", "ALT", "A1", "OBS_CT", "ADD_OR", "ADD_LOG_OR_SE", "ADD_Z_STAT", "ADD_P")], filtered_output2[c("INTERACTION_OR", "INTERACTION_LOG_OR_SE", "INTERACTION_Z_STAT", "INTERACTION_P")])
  rownames(out) <- NULL 
  return(out)
 }

