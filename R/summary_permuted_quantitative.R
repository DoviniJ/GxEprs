#' summary_permuted_quantitative function
#' This function uses plink2 and outputs the p value of permuted model in the target dataset, using pre-generated GWAS and GWEIS summary statistics files named Q_trd.sum, Q_add.sum and Q_gxe.sum
#' @param Qphe_target Phenotype file containing family ID, individual ID and phenotype of the target dataset as columns, without heading
#' @param Qcov_target Covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the target dataset as columns, without heading
#' @param iterations Number of iterations used in permutation
#' @param add_score The .sscore file generated using additive SNP effects of GWEIS summary statistics
#' @param gxe_score The .sscore file generated using interaction SNP effects of GWEIS summary statistics
#' @keywords regression summary
#' @export 
#' @importFrom stats binomial fitted.values glm lm
#' @importFrom utils read.table write.table
#' @return This function will output
#' \item{Q_permuted_p.txt}{the p value of the permuted model}
#' @examples \dontrun{ 
#' x <- summary_permuted_quantitative(Qphe_target, Qcov_target)
#' x 
#' }
summary_permuted_quantitative <- function(Qphe_target, Qcov_target, iterations = 1000, add_score = "Q_add.sscore", gxe_score = "Q_gxe.sscore"){
  cov_file <- read.table(Qcov_target)
  n_confounders = ncol(cov_file) - 4
  fam=read.table(Qphe_target, header=F) 
  colnames(fam) <- c("FID", "IID", "PHENOTYPE")
  dat=read.table(Qcov_target, header=F)
  colnames(dat)[1] <- "FID"
  colnames(dat)[2] <- "IID"
  prs1_all=read.table(add_score)
  colnames(prs1_all)[1] <- "FID"
  colnames(prs1_all)[2] <- "IID"
  prs1=merge(fam, prs1_all, by = "FID")
  prs2_all=read.table(gxe_score)
  colnames(prs2_all)[1] <- "FID"
  colnames(prs2_all)[2] <- "IID"
  prs2=merge(fam, prs2_all, by = "FID")
  m1 <- match(dat$IID, prs1$IID.x)
  out = fam$PHENOTYPE[m1]
  cov=scale(dat$V3[m1])
  ps1=scale(prs1$V5)
  ps2=scale(prs2$V5)
  xv1=scale(prs1$V5*cov)
  xv2=scale(prs2$V5*cov)
  cov2=scale(cov^2)
  if(n_confounders == 0){
    pn=iterations; pp_gxe_x_E=0
    df_regular_new <- as.data.frame(cbind(out, cov, ps1, ps2, xv2))
    regular_m = lm(out ~., data = df_regular_new)     
    regular_p = summary(regular_m)$coefficients[5,4]
    for(i in 1:pn){
    percentage <- (i / iterations) * 100
    percentage <- (i / iterations) * 100
    cat(sprintf("\rProgress: %3.0f%%", percentage))
      sv2=sample(seq(1, length(out)))
      xv2=scale(prs2$V5[sv2]*cov)
    df_new <- as.data.frame(cbind(out, cov, ps1, ps2, xv2))
    colnames(df_new)[1] <- "out"
    colnames(df_new)[2] <- "E"
    colnames(df_new)[3] <- "PRS_add"
    colnames(df_new)[4] <- "PRS_gxe"
    colnames(df_new)[5] <- "PRS_gxe x E"
    m = lm(out ~., data = df_new)
    if (regular_p < summary(m)$coefficients[5,4]) pp_gxe_x_E=pp_gxe_x_E+1
    }
    sink("Q_permuted_p.txt")
    print(pp_gxe_x_E/pn)
    sink()}else{
    pn=iterations; pp_gxe_x_E=0
    df_regular_new <- as.data.frame(cbind(out, cov, ps1, ps2, xv2))
    regular_m = lm(out ~., data = df_regular_new)     
    regular_p = summary(regular_m)$coefficients[5,4]
    conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
    for (k in 1:n_confounders) {
      conf_var[, k] <- as.numeric(dat[, k+4])
    }
    conf_var <- conf_var[m1,]
    for(i in 1:pn){
    percentage <- (i / iterations) * 100
    cat(sprintf("\rProgress: %3.0f%%", percentage))
      sv2=sample(seq(1, length(out)))
      xv2=scale(prs2$V5[sv2]*cov)
      df_new <- as.data.frame(cbind(out, cov, ps1, ps2, xv2, conf_var))
      colnames(df_new)[1] <- "out"
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_gxe"
      colnames(df_new)[5] <- "PRS_gxe x E"
      m = lm(out ~., data = df_new)
      if (regular_p < summary(m)$coefficients[5,4]) pp_gxe_x_E=pp_gxe_x_E+1
    }
  }
  cat("\n")
  return(pp_gxe_x_E/pn)
}
