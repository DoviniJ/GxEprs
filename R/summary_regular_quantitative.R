#' summary_regular_quantitative function
#' This function uses plink2 and outputs the summary of regular model in the target dataset, using pre-generated GWAS and GWEIS summary statistics files named Q_trd.sum, Q_add.sum and Q_gxe.sum
#' @param Qphe_target Phenotype file containing family ID, individual ID and phenotype of the target dataset as columns, without heading
#' @param Qcov_target Covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the target dataset as columns, without heading
#' @param trd_score The .sscore file generated using additive SNP effects of GWAS summary statistics
#' @param add_score The .sscore file generated using additive SNP effects of GWEIS summary statistics
#' @param gxe_score The .sscore file generated using interaction SNP effects of GWEIS summary statistics
#' @param Model Specify the model number (1: y = PRS_trd + E + PRS_trd x E + confounders, 2: y = PRS_add + E + PRS_add x E + confounders, 3: y = PRS_add + E + PRS_gxe x E + confounders, 4: y = PRS_add + E + PRS_gxe + PRS_gxe x E + confounders, where y is the outcome variable, E is the covariate of interest, PRS_trd and PRS_add are the polygenic risk scores computed using additive SNP effects of GWAS and GWEIS summary statistics respectively, and PRS_gxe is the polygenic risk scores computed using GxE interaction SNP effects of GWEIS summary statistics.)
#' @param summary_output Name of the model summary file specified by the user
#' @param risk_output Name of the file containing risk scores of the target individuals specified by the user
#' @keywords regression summary risk scores
#' @export 
#' @importFrom stats binomial fitted.values glm lm
#' @importFrom utils read.table write.table
#' @return This function will output
#' \item{Qsummary.txt}{the summary of the fitted model}
#' \item{Individual_risk_values.txt}{the estimated risk values of individuals in the target sample}
#' @examples \dontrun{
#' x <- summary_regular_quantitative(Qphe_target, Qcov_target, 
#'                                   add_score = "Q_add.sscore", 
#'                                   gxe_score = "Q_gxe.sscore", 
#'                                   Model = 4, 
#'                                   summary_output = "Qsummary.txt", 
#'                                   risk_output = "Individual_risk_values.txt")
#' x[[1]][[1]] #to obtain the content of the Qsummary.txt file
#' x[[1]][[2]] #to extract "Call" of the model summary
#' x[[1]][[3]] #to extract terms of the model summary
#' x[[1]][[4]] #to extract the residuals
#' x[[1]][[5]] #to extract regrerssion coefficients of the model summary
#' x[[1]][[6]] #to extract aliesed information of the model summary
#' x[[1]][[7]] #to extract "sigma" (residual standard error) information 
#'             #of the model summary
#' x[[1]][[8]] #to extract degrees of freedom of the model summary
#' x[[1]][[9]] #to extract the R squared value of the model summary
#' x[[1]][[10]] #to extract the adjusted R squared value of the model summary
#' x[[1]][[11]] #to extract the test statistic values of the model summary
#' x[[1]][[12]] #to extract unscaled variance covariance matrix of all variables
#' head(x[[2]]) #to read the head of all columns in Individual_risk_values.txt file
#' x[[2]][,1] #to extract the column containing family ID's 
#' x[[2]][,2] #to extract the column containing individual ID's 
#' x[[2]][,3] #to extract the column containing predicted risk scores 
#' }
summary_regular_quantitative <- function(Qphe_target, Qcov_target, trd_score = "Q_trd.sscore", add_score = "Q_add.sscore", gxe_score = "Q_gxe.sscore", Model, summary_output = "Qsummary.txt", risk_output = "Individual_risk_values.txt"){
  cov_file <- read.table(Qcov_target)
  n_confounders = ncol(cov_file) - 4
  fam=read.table(Qphe_target, header=F) 
  colnames(fam) <- c("FID", "IID", "PHENOTYPE")
  dat=read.table(Qcov_target ,header=F)
  colnames(dat)[1] <- "FID"
  colnames(dat)[2] <- "IID"
  if(file.exists(trd_score)){
    prs0_all=read.table(trd_score)
    colnames(prs0_all)[1] <- "FID"
    colnames(prs0_all)[2] <- "IID"
    prs0=merge(fam, prs0_all, by = "FID")
    m1 <- match(dat$IID, prs0$IID.x)
    ps0=scale(prs0$V5)
    out = scale(fam$PHENOTYPE[m1])
    cov=scale(dat$V3[m1])
    xv0=scale(prs0$V5*cov)
  }
  if(file.exists(add_score)){
    prs1_all=read.table(add_score)
    colnames(prs1_all)[1] <- "FID"
    colnames(prs1_all)[2] <- "IID"
    prs1=merge(fam, prs1_all, by = "FID")
    m1 <- match(dat$IID, prs1$IID.x)
    ps1=scale(prs1$V5)
    out = scale(fam$PHENOTYPE[m1])
    cov=scale(dat$V3[m1])
    xv1=scale(prs1$V5*cov)
  }
  if(file.exists(gxe_score)){
    prs2_all=read.table(gxe_score)
    colnames(prs2_all)[1] <- "FID"
    colnames(prs2_all)[2] <- "IID"
    prs2=merge(fam, prs2_all, by = "FID")
    m1 <- match(dat$IID, prs2$IID.x)
    ps2=scale(prs2$V5)
    out = scale(fam$PHENOTYPE[m1])
    cov=scale(dat$V3[m1])
    xv2=scale(prs2$V5*cov)
  }
  if(Model == 1){
    if(n_confounders == 0){
      df_new <- as.data.frame(cbind(out, cov, ps0, xv0))
      colnames(df_new)[1] <- "out"
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_trd"
      colnames(df_new)[4] <- "PRS_trd x E"
      m = lm(out ~., data = df_new)
      sink(summary_output)
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink(risk_output)
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }else{
      conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
      for (k in 1:n_confounders) {
        conf_var[, k] <- as.numeric(dat[, k+4])
      }
      conf_var <- conf_var[m1,]
      df_new <- as.data.frame(cbind(out, cov, ps0, xv0, conf_var))
      colnames(df_new)[1] <- "out"
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_trd"
      colnames(df_new)[4] <- "PRS_trd x E"
      m = lm(out ~., data = df_new)
      sink(summary_output)
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink(risk_output)
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }
    s <- summary(m)
    out1 <- list(s, s$call, s$terms, s$residuals, s$coefficients, s$aliesed, s$sigma, s$df, s$r.squared, s$adj.r.squared, s$fstatistics, s$cov.unscaled)
    out2 <- cbind(fam$FID[m1], fam$IID[m1], m_fit)
    out_all <- list(out1, out2)
  }
  if(Model == 2){
    if(n_confounders == 0){
      df_new <- as.data.frame(cbind(out, cov, ps1, xv1))
      colnames(df_new)[1] <- "out"
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_add x E"
      m = lm(out ~., data = df_new)
      sink(summary_output)
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink(risk_output)
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }else{
      conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
      for (k in 1:n_confounders) {
        conf_var[, k] <- as.numeric(dat[, k+4])
      }
      conf_var <- conf_var[m1,]
      df_new <- as.data.frame(cbind(out, cov, ps1, xv1, conf_var))
      colnames(df_new)[1] <- "out"
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_add x E"
      m = lm(out ~., data = df_new)
      sink(summary_output)
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink(risk_output)
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }
    s <- summary(m)
    out1 <- list(s, s$call, s$terms, s$residuals, s$coefficients, s$aliesed, s$sigma, s$df, s$r.squared, s$adj.r.squared, s$fstatistics, s$cov.unscaled)
    out2 <- cbind(fam$FID[m1], fam$IID[m1], m_fit)
    out_all <- list(out1, out2)
  }
  if(Model == 3){
    if(n_confounders == 0){
      df_new <- as.data.frame(cbind(out, cov, ps1, xv2))
      colnames(df_new)[1] <- "out"
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_gxe x E"
      m = lm(out ~., data = df_new)
      sink(summary_output)
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink(risk_output)
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }else{
      conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
      for (k in 1:n_confounders) {
        conf_var[, k] <- as.numeric(dat[, k+4])
      }
      conf_var <- conf_var[m1,]
      df_new <- as.data.frame(cbind(out, cov, ps1, xv2, conf_var))
      colnames(df_new)[1] <- "out"
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_gxe x E"
      m = lm(out ~., data = df_new)
      sink(summary_output)
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink(risk_output)
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }
    s <- summary(m)
    out1 <- list(s, s$call, s$terms, s$residuals, s$coefficients, s$aliesed, s$sigma, s$df, s$r.squared, s$adj.r.squared, s$fstatistics, s$cov.unscaled)
    out2 <- cbind(fam$FID[m1], fam$IID[m1], m_fit)
    out_all <- list(out1, out2)
  }
  if(Model == 4){
    if(n_confounders == 0){
      df_new <- as.data.frame(cbind(out, cov, ps1, ps2, xv2))
      colnames(df_new)[1] <- "out"
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_gxe"
      colnames(df_new)[5] <- "PRS_gxe x E"
      m = lm(out ~., data = df_new)
      sink(summary_output)
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink(risk_output)
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }else{
      conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
      for (k in 1:n_confounders) {
        conf_var[, k] <- as.numeric(dat[, k+4])
      }
      conf_var <- conf_var[m1,]
      df_new <- as.data.frame(cbind(out, cov, ps1, ps2, xv2, conf_var))
      colnames(df_new)[1] <- "out"
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_gxe"
      colnames(df_new)[5] <- "PRS_gxe x E"
      m = lm(out ~., data = df_new)
      sink(summary_output)
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink(risk_output)
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }
    s <- summary(m)
    out1 <- list(s, s$call, s$terms, s$residuals, s$coefficients, s$aliesed, s$sigma, s$df, s$r.squared, s$adj.r.squared, s$fstatistics, s$cov.unscaled)
    out2 <- cbind(fam$FID[m1], fam$IID[m1], m_fit)
    out_all <- list(out1, out2)
  }
   return(out_all)
}
