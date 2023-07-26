#' GWEIS_quantitative function
#' This function performs GWEIS using plink2 and outputs the GWEIS summary statistics with additive SNP effects and interaction SNP effects separately. It is recommended to save the outputs in separate user-specified files (see examples).
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param Qphe_discovery Phenotype file containing family ID, individual ID and phenotype of the discovery dataset as columns, without heading
#' @param Qcov_discovery Covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the discovery dataset as columns, without heading
#' @param thread Number of threads used
#' @keywords gwies interaction gxe
#' @export 
#' @importFrom stats binomial fitted.values glm lm
#' @importFrom utils read.table 
#' @return This function will perform GWEIS and output
#' \item{Q_out.add.sum}{GWEIS summary statistics with additive SNP effects}
#' \item{Q_out.gxe.sum}{GWEIS summary statistics with interaction SNP effects}
#' @examples \dontrun{ 
#' x <- GWEIS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery, 
#' thread = 20)
#' sink("Q_out.add.sum") #to create a file in the working directory
#' write.table(x[[1]], sep = " ", row.names = FALSE, quote = FALSE) #to write the output
#' sink() #to save the output
#' sink("Q_out.gxe.sum") #to create a file in the working directory
#' write.table(x[[2]], sep = " ", row.names = FALSE, quote = FALSE) #to write the output
#' sink() #to save the output
#' head(x[[1]]) #to read the head of all columns in GWEIS summary statistics of  
#' additive SNP effects 
#' x[[1]]$V1 #to extract the chromosome number (CHROM)
#' x[[1]]$V2 #to extract the base pair position (POS)
#' x[[1]]$V3 #to extract the SNP ID (ID)
#' x[[1]]$V4 #to extract the reference allele (REF)
#' x[[1]]$V5 #to extract the alternate allele (ALT)
#' x[[1]]$V6 #to extract the minor allele (A1)
#' x[[1]]$V7 #to extract the type of test performed (TEST)
#' x[[1]]$V8 #to extract the nmber of allele observations (OBS_CT)
#' x[[1]]$V9 #to extract the SNP effect (BETA)
#' x[[1]]$V10 #to extract the standard error of each SNP effect (SE)
#' x[[1]]$V11 #to extract the test statistic (T_STAT)
#' x[[1]]$V12 #to extract the p value (P)
#' x[[1]]$V13 #to extract the error code (ERRCODE)
#' head(x[[2]]) #to read the head of all columns in GWEIS summary statistics of 
#' interaction SNP effects 
#' x[[2]]$V1 #to extract the chromosome number (CHROM)
#' x[[2]]$V2 #to extract the base pair position (POS)
#' x[[2]]$V3 #to extract the SNP ID (ID)
#' x[[2]]$V4 #to extract the reference allele (REF)
#' x[[2]]$V5 #to extract the alternate allele (ALT)
#' x[[2]]$V6 #to extract the minor allele (A1)
#' x[[2]]$V7 #to extract the type of test performed (TEST)
#' x[[2]]$V8 #to extract the number of allele observations (OBS_CT)
#' x[[2]]$V9 #to extract the SNP effect (BETA)
#' x[[2]]$V10 #to extract the standard error of each SNP effect (SE)
#' x[[2]]$V11 #to extract the test statistic (T_STAT)
#' x[[2]]$V12 #to extract the p value (P)
#' x[[2]]$V13 #to extract the error code (ERRCODE)
#' }
GWEIS_quantitative <- function(plink_path, b_file, Qphe_discovery, Qcov_discovery, thread = 20){
  cov_file <- read.table(Qcov_discovery)
  n_confounders = ncol(cov_file) - 4
  if(n_confounders > 0){
    parameters <- c(1, 2, (1:(n_confounders+1))+3)
  }
  else{
    parameters <- c(1, 2, 4)
  }
  param_vec <- paste0(parameters, collapse = ", ")
  runPLINK <- function(PLINKoptions = "") system(paste(plink_path, PLINKoptions))
  log_file <- runPLINK(paste0(" --bfile ", b_file,
                " --glm interaction --pheno ", 
                Qphe_discovery, 
                " --covar ", Qcov_discovery, 
                " --parameters ", param_vec, 
                " --allow-no-sex --threads ", 
                thread,
                " --out ", tempdir(),"/Q_gweis"))
  plink_output <- read.table(paste0(tempdir(), "/Q_gweis.PHENO1.glm.linear"), header = FALSE)
  filtered_output <- as.data.frame(plink_output[(plink_output$V7=="ADD"),])
  filtered_output2 <- plink_output[(plink_output$V7=="ADDxCOVAR1"),]
  Q_out.add.sum <- filtered_output
  Q_out.gxe.sum <- filtered_output2
 out <- list(Q_out.add.sum, Q_out.gxe.sum)
 return(out)
 }

