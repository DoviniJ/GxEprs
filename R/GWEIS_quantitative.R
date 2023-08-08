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
#' \item{Q_out.sum}{GWEIS summary statistics with additive and interaction SNP effects}
#' @examples \dontrun{ 
#' x <- GWEIS_quantitative (plink_path, DummyData, Qphe_discovery, Qcov_discovery, 
#' thread = 20)
#' sink("Q_out.add.sum") #to create a file in the working directory
#' write.table(x[c("ID", "A1", "ADD_BETA")], sep = " ", 
#' row.names = FALSE, quote = FALSE) #to write the output
#' sink() #to save the output
#' sink("Q_out.gxe.sum") #to create a file in the working directory
#' write.table(x[c("ID", "A1", "INTERACTION_BETA")], sep = " ", 
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
#' x$ADD_BETA #to extract the additive SNP effects
#' x$ADD_SE #to extract the standard errors of the 
#' additive SNP effects
#' x$ADD_T_STAT #to extract the test statistics of additive 
#' SNP effects
#' x$ADD_P #to extract the p values of additive SNP effects
#' x$INTERACTION_BETA #to extract the interaction SNP effects
#' x$INTERACTION_SE #to extract the standard errors of the 
#' interaction SNP effects
#' x$INTERACTION_T_STAT #to extract the test statistics of 
#' interaction SNP effects
#' x$INTERACTION_P #to extract the p values of interaction 
#' SNP effects
#' }
GWEIS_quantitative <- function(plink_path, b_file, Qphe_discovery, Qcov_discovery, thread = 20){
  os_name <- Sys.info()["sysname"]
   if (startsWith(os_name, "Win")) {
     slash <- paste0("\\")
   } else {
     slash <- paste0("/")
   }
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
                " --out ", tempdir(), slash, "Q_gweis"))
  first_line <- readLines(paste0(tempdir(), slash, "Q_gweis.PHENO1.glm.linear"), n = 1)
  col_names <- strsplit(first_line, "\t")[[1]]
  col_names[1] <- sub("#", "", col_names[1])
  plink_output <- read.table(paste0(tempdir(), slash, "Q_gweis.PHENO1.glm.linear"), skip = 1, col.names = col_names, sep = "\t")
  filtered_output <- as.data.frame(plink_output[(plink_output$TEST=="ADD"),])
  colnames(filtered_output)[c(grep("BETA", colnames(filtered_output)), grep("SE", colnames(filtered_output)), 
                           grep("T_STAT", colnames(filtered_output)), grep("^\\bP\\b$", colnames(filtered_output)), 
                           grep("ERRCODE", colnames(filtered_output)))] <- c("ADD_BETA", "ADD_SE", "ADD_T_STAT", "ADD_P", "ADD_ERRCODE")
  filtered_output2 <- plink_output[(plink_output$TEST=="ADDxCOVAR1"),]
  colnames(filtered_output2)[c(grep("BETA", colnames(filtered_output2)), grep("SE", colnames(filtered_output2)), 
                           grep("T_STAT", colnames(filtered_output2)), grep("^\\bP\\b$", colnames(filtered_output2)), 
                           grep("ERRCODE", colnames(filtered_output2)))] <- c("INTERACTION_BETA", "INTERACTION_SE", "INTERACTION_T_STAT", 
										"INTERACTION_P", "INTERACTION_ERRCODE")
  out <- cbind(filtered_output[c("CHROM", "POS", "ID", "REF", "ALT", "A1", "OBS_CT", "ADD_BETA", "ADD_SE", "ADD_T_STAT", "ADD_P")], filtered_output2[c("INTERACTION_BETA", "INTERACTION_SE", "INTERACTION_T_STAT", "INTERACTION_P")])
  rownames(out) <- NULL 
  return(out)
 }

