#binary files: mydata.fam, mydata.bim and mydata.bed
#phenotype file of discovery sample (quantitative outcome): Qpd.txt
#covariate file of discovery sample (quantitative outcome): Qcd.txt
#phenotype file of target sample (quantitative outcome): Qpt.txt
#covariate file of the target sample (quantitative outcome): Qct.txt

plink_path <- "<plink_path>/plink2" #Give the path (<plink_path>) where plink executable file is located

library(devtools)
install_github("DoviniJ/GxEprs") 
library(GxEprs)

###########
#Example 1: Steps to estimate risk values for individuals in the target dataset only by using functions available in GxEprs package
###########

#outcome variable: Quantitative trait
#target model: 4

a <- GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt") #perform GWEIS
add <- a[c("ID", "A1", "ADD_BETA")] #store the additive SNP effects from GWEIS
gxe <- a[c("ID", "A1", "INTERACTION_BETA")] #store the interaction SNP effects from GWEIS
p <- PRS_quantitative(plink_path, "mydata", summary_input = add]) #obtain PRSs by using additive SNP effects from GWEIS
q <- PRS_quantitative(plink_path, "mydata", summary_input = gxe) #obtain PRSs by using interaction effects from GWEIS

x <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = p, gxe_score = q, Model = 4) 
x$summary #contains the target regular model summary output, when the outcome is quantitative.
x$risk.values #contains all the calculated individual risk scores using the target dataset (e.g. Model 4), when the outcome is quantitative. 

summary_permuted_quantitative("qpt.txt", "qct.txt", iterations = 1000, add_score = p, gxe_score = q)

#Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models 
#(using summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = q, gxe_score = r, Model = 4) and 
#summary_permuted_quantitative("qpt.txt", "qct.txt", iterations = 1000, add_score = q, gxe_score = r)), if you 
#choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 4) when generating risk scores. 
#If the 'PRS_gxe x E' term is significant in Model 4, and insignificant in Model 4* (permuted p value), consider 
#that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model).



###########
#Example 2: Steps to estimate risk values for individuals in the target dataset by using your own GWEIS summary statistics
###########

#outcome variable: Quantitative trait
#target model: 4

a <- read.table("add_sum.txt", header = T) #read the additive SNP effects from GWEIS summary statistics (columns should contain Snp ID, A1 allele and BETA only)
b <- read.table("gxe_sum.txt", header = T) #read the interaction effects from GWEIS summary statistics (columns should contain Snp ID, A1 allele and BETA only)
p <- PRS_quantitative(plink_path, "mydata", summary_input = a)
q <- PRS_quantitative(plink_path, "mydata", summary_input = b)

x <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = p, gxe_score = q, Model = 4) 
x$summary 
x$risk.values  

summary_permuted_quantitative("qpt.txt", "qct.txt", iterations = 1000, add_score = p, gxe_score = q)



###########
#Example 3: Steps to estimate risk values for individuals in the target dataset by using your own PRSs
###########

#outcome variable: Quantitative trait
#target model: 4

p <- read.table("add_prs.txt", header = T) #read the PRSs generated using additive SNP effects (columns should contain family IDs, individual IDs and PRSs and the column names should be "FID", "IID" and "PRS")
q <- read.table("gxe_prs.txt", header = T) #read the PRSs generated using interaction SNP effects (columns should contain family IDs, individual IDs and PRSs and the column names should be "FID", "IID" and "PRS")

x <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = p, gxe_score = q, Model = 4) 
x$summary 
x$risk.values  

summary_permuted_quantitative("qpt.txt", "qct.txt", iterations = 1000, add_score = p, gxe_score = q)
