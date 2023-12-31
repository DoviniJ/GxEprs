% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_permuted_quantitative.R
\name{summary_permuted_quantitative}
\alias{summary_permuted_quantitative}
\title{summary_permuted_quantitative function
This function outputs the p value of permuted model in the target dataset, using pre-generated Polygenic Risk Scores (PRSs) of all the individuals. Note that the input used in this function can be generated by using PRS_quantitative functions. It is recommended to run this function, if you choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 4) when generating risk scores. If the 'PRS_gxe x E' term is significant in Model 4, and insignificant in Model 4* (permuted p value), consider that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model).}
\usage{
summary_permuted_quantitative(
  Qphe_target,
  Qcov_target,
  iterations = 1000,
  add_score,
  gxe_score
)
}
\arguments{
\item{Qphe_target}{Phenotype file containing family ID, individual ID and phenotype of the target dataset as columns, without heading}

\item{Qcov_target}{Covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the target dataset as columns, without heading}

\item{iterations}{Number of iterations used in permutation}

\item{add_score}{PRSs generated using additive SNP effects of GWEIS summary statistics}

\item{gxe_score}{PRSs generated using interaction SNP effects of GWEIS summary statistics}
}
\value{
This function will output
\item{Q_permuted_p}{the p value of the permuted model}
}
\description{
summary_permuted_quantitative function
This function outputs the p value of permuted model in the target dataset, using pre-generated Polygenic Risk Scores (PRSs) of all the individuals. Note that the input used in this function can be generated by using PRS_quantitative functions. It is recommended to run this function, if you choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 4) when generating risk scores. If the 'PRS_gxe x E' term is significant in Model 4, and insignificant in Model 4* (permuted p value), consider that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model).
}
\examples{
\dontrun{ 
a <- GWEIS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery)
add <- a[c("ID", "A1", "ADD_BETA")]
gxe <- a[c("ID", "A1", "INTERACTION_BETA")]
p <- PRS_quantitative(plink_path, DummyData, summary_input = add)
q <- PRS_quantitative(plink_path, DummyData, summary_input = gxe)
x <- summary_permuted_quantitative(Qphe_target, Qcov_target, iterations = 1000, 
add_score = p, gxe_score = q)
x 
}
}
\keyword{model}
\keyword{permuted}
\keyword{pvalue}
