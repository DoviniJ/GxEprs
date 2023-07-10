#' GWEIS_binary function
#' This function performs GWEIS using plink2 and outputs the GWEIS summary statistics files with additive SNP effects named B_out.add.sum and interaction SNP effects named B_out.gxe.sum
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param Bphe_discovery Phenotype file containing family ID, individual ID and phenotype of the discovery dataset as columns, without heading
#' @param Bcov_discovery Covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the discovery dataset as columns, without heading
#' @param thread Number of threads used
#' @param summary_output Name (prefix) of the additive or interaction SNP effects of the GWEIS summary statistics file specified by the user
#' @keywords gwies, interaction, gxe
#' @export 
#' @importFrom stats binomial fitted.values glm lm
#' @importFrom utils read.table write.table
#' @return This function will perform GWEIS and output
#' \item{B_out.add.sum}{GWEIS summary statistics file with additive SNP effects}
#' \item{B_out.gxe.sum}{GWEIS summary statistics file with interaction SNP effects}
#' @examples \dontrun{ 
#' x <- GWEIS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery, 
#'                   thread = 20, summary_output = "B_out")
#' head(x[[1]]) #to extract the head of all columns in B_out.add.sum file
#' x[[1]]$V1 #to extract the chromosome number (CHROM)
#' x[[1]]$V2 #to extract the base pair position (POS)
#' x[[1]]$V3 #to extract the SNP ID (ID)
#' x[[1]]$V4 #to extract the reference allele (REF)
#' x[[1]]$V5 #to extract the alternate allele (ALT)
#' x[[1]]$V6 #to extract the minor allele (A1)
#' x[[1]]$V7 #to extract whether firth regression is used (FIRTH?)
#' x[[1]]$V8 #to extract the type of test performed (TEST)
#' x[[1]]$V9 #to extract the nmber of allele observations (OBS_CT)
#' x[[1]]$V10 #to extract the odds ration of the SNP effect (OR)
#' x[[1]]$V11 #to extract the standard error of log odds (LOG(OR)_SE)
#' x[[1]]$V12 #to extract the test statistic (Z_STAT)
#' x[[1]]$V13 #to extract the p value (P)
#' x[[1]]$V14 #to extract the error code (ERRCODE)
#' head(x[[2]]) #to extract the head of all columns in B_out.gxe.sum file
#' x[[2]]$V1 #to extract the chromosome number (CHROM)
#' x[[2]]$V2 #to extract the base pair position (POS)
#' x[[2]]$V3 #to extract the SNP ID (ID)
#' x[[2]]$V4 #to extract the reference allele (REF)
#' x[[2]]$V5 #to extract the alternate allele (ALT)
#' x[[2]]$V6 #to extract the minor allele (A1)
#' x[[2]]$V7 #to extract whether firth regression is used (FIRTH?)
#' x[[2]]$V8 #to extract the type of test performed (TEST)
#' x[[2]]$V9 #to extract the nmber of allele observations (OBS_CT)
#' x[[2]]$V10 #to extract the odds ration of the SNP effect (OR)
#' x[[2]]$V11 #to extract the standard error of log odds (LOG(OR)_SE)
#' x[[2]]$V12 #to extract the test statistic (Z_STAT)
#' x[[2]]$V13 #to extract the p value (P)
#' x[[2]]$V14 #to extract the error code (ERRCODE)
#' }
GWEIS_binary <- function(plink_path, b_file, Bphe_discovery, Bcov_discovery, thread = 20, summary_output = "B_out"){
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
  runPLINK(paste0(" --bfile ", b_file,
                " --glm interaction --pheno ", 
                Bphe_discovery, 
                " --covar ", Bcov_discovery, 
                " --parameters ", param_vec, 
                " --allow-no-sex --threads ", 
                thread,
                " --out B_gweis"))
  plink_output <- read.table("B_gweis.PHENO1.glm.logistic.hybrid", header = FALSE)
  filtered_output <- as.data.frame(plink_output[(plink_output$V8=="ADD"),])
  filtered_output$V10 = log(filtered_output$V10)
  filtered_output2 <- plink_output[(plink_output$V8=="ADDxCOVAR1"),]
  filtered_output2$V10 <- log(filtered_output2$V10)
  summary_output1 <- paste0(noquote(summary_output), ".add.sum")
  summary_output2 <- paste0(noquote(summary_output), ".gxe.sum")
  sink(summary_output1)
  write.table(filtered_output, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
  sink(summary_output2)
  write.table(filtered_output2, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
 out <- list(filtered_output, filtered_output2)
 return(out)
 }

