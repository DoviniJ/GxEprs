---
title: "GxEprs"
author: "Dovini Jayasinghe, Md Moksedul Momin and Hong Lee"
date: "13-07-2023"
output: pdf_document
---

# GxEprs
The 'GxEprs' is an R package to detect and estimate GxE. It uses a novel PRS model that can enhance the prediction accuracy by utilising GxE effects. Firstly it performs Genome Wide Association Studies (GWAS)  and Genome Wide Environment Interaction Studies (GWEIS) using the discovery dataset (see functions ```GWAS_binary()```,```GWAS_quantitative()```, ```GWEIS_binary()```, ```GWEIS_quantitative()```). See the section $\color{red}{IMPORTANT}$ for the discovery models used. Secondly, it uses the GWAS and GWEIS summary statistics generated from the fucntions above to obtain polygenic risk scores (PRSs) (see functions ```PRS_binary()``` and ```PRS_quantitative()```) for the target sample. Finally it predicts the risk values of each individual in the target sample (see functions ```summary_regular_binary()``` and ```summary_regular_quantitative()```). Note that the users can fit 4 different models when the outcome is a quantitative trait, and 5 different models when the outcome is a binary disease trait. See the section $\color{red} {IMPORTANT}$ for the target models used. Finally, it is recommended to check the p-value from permutations using Model 4 (see function ```summary_permuted_quantitative()```), and Model 5 (see function```summary_permuted_binary()```), to make sure that the significance of GxE is not spurious due to model misspecification (see references).

# Data preparation

## File formats
### Input files
1) mydata.fam - This is a file associated with the PLINK binary format file which contains the following columns in order. The example dataset has 1,000 individuals. Note that the file has no column headings. This follows the PLINK .fam file format.
* family ID (FID) 
* individual ID (IID) 
* father's ID 
* mother's ID 
* sex 
* phenotype value

```
ID_1 ID_1 0 0 1 -9
ID_2 ID_2 0 0 2 -9
ID_3 ID_3 0 0 2 -9
ID_4 ID_4 0 0 2 -9
ID_5 ID_5 0 0 1 -9
```
  
2) mydata.bim - This is is a file associated with the PLINK binary format file which contains the following columns in order. The example dataset has 1,000 SNPs. Note that the file has no column headings. This follows the PLINK .bim file format.
* chromosome code 
* SNP ID 
* position of centimorgans 
* base-pair coordinate 
* minor allele  
* reference allele 

```
1	SNP_1	0	768448	A	G
1	SNP_2	0	853954	C	A
1	SNP_3	0	880390	A	C
1	SNP_4	0	940203	A	G
1	SNP_5	0	987670	T	G
```

3) mydata.bed - This is the PLINK binary format file which includes genotype information. This follows the PLINK .bed file format.
4) Bpd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 800 individuals. Note that the file has no column headings.
* FID 
* IID  
* binary phenotype (1=controls, 2=cases) of the discovery sample

```
ID_1 ID_1 1
ID_2 ID_2 2
ID_3 ID_3 1
ID_4 ID_4 1
ID_5 ID_5 1
```

5) Bcd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 800 individuals. Note that the file has no column headings.    
* FID 
* IID 
* standardized covariate 
* square of the standardized covariate (Note: This column is required to address model misspecification, if present. See reference for details.) 
* 14 confounders of the discovery sample

```
ID_1 ID_1 0.787403812314451 0.620004763647331 -3.04026 45 -12.048 2.17634 -0.940322 -0.446351 -5.45685 -2.53161 -2.13435 -1.95623 -3.82792 -0.380636 0 10
ID_2 ID_2 -0.119722138532781 0.0143333904548625 -5.1054 64 -14.5169 6.01889 -3.85803 3.62625 5.10717 -3.54574 0.393994 3.64275 4.42975 -2.26704 1 19
ID_3 ID_3 0.173372721351375 0.0300581005087816 -1.91044 59 -12.7462 5.95244 0.0738914 1.80523 4.76284 0.130369 -1.05615 0.316777 0.988783 -1.76502 1 7
ID_4 ID_4 -0.699321184695051 0.48905011936329 -1.83526 68 -10.3349 4.71264 -1.84521 -0.524855 -3.80275 0.837965 0.265233 2.10903 -0.210259 0.71504 0 20
ID_5 ID_5 3.69300366651739 13.6382760809109 -3.15649 69 -8.56737 4.78248 -1.49547 -7.49413 -5.39887 1.85316 4.07476 1.05351 0.825942 -2.09669 1 20
```

6) Bpt.txt - This is a .txt file which contains the following columns in order. The target dataset has 200 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID  
* binary phenotype (0=controls, 1=cases) of the target sample

```
ID_801 ID_801 0
ID_802 ID_802 1
ID_803 ID_803 0
ID_804 ID_804 0
ID_805 ID_805 1
```

7) Bct.txt - This is a .txt file which contains the following columns in order. The target dataset has 200 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID 
* standardized covariate 
* square of the standardized covariate  (Note: This column is required to address model misspecification, if present. See reference for details.)
* 14 confounders of the target sample

```
ID_801 ID_801 -0.420822976931972 0.177091977913887 -4.12263 57 -13.5185 5.40198 -4.81994 0.664494 -4.92217 -0.451329 3.14677 0.42704 0.821306 -2.77705 1 7
ID_802 ID_802 -0.0805583280660987 0.00648964422080519 2.92534 56 -13.6236 3.21643 -0.856048 0.750187 -2.01798 -0.350832 5.10141 2.1807 -6.04343 1.78928 1 19
ID_803 ID_803 -1.32752644870189 1.76232647200306 -3.09118 61 -9.94475 3.60562 -0.917639 0.905664 -5.09843 -1.16329 -1.88102 -1.24154 0.699574 2.2442 0 20
ID_804 ID_804 0.698007239555549 0.487214106471958 4.58829 49 -12.5471 4.09467 -2.58951 6.06898 12.9822 -0.704179 2.90357 -0.334968 5.04274 0.66175 0 10
ID_805 ID_805 -0.657981606980219 0.432939795124272 -3.53948 56 -12.795 2.91524 -2.72794 3.61555 3.92957 -2.93899 -0.454737 2.31013 2.51783 -4.15592 0 7
```

8) Qpd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 800 individuals. Note that the file has no column headings.
* FID 
* IID  
* quantitative phenotype of the discovery sample

```
ID_1 ID_1 31.6534
ID_2 ID_2 25.5035
ID_3 ID_3 26.7391
ID_4 ID_4 25.5271
ID_5 ID_5 26.7165
```

9) Qcd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 800 individuals. Note that the file has no column headings.    
* FID 
* IID 
* standardized covariate 
* square of the standardized covariate  (Note: Although squared covariate term is not involved when analyzing quantitative traits, this file can/may be used as the covariate file when modulating with a binary outcome. In that event, this column is required to address model misspecification, if present. See reference for details.)
* 14 confounders of the discovery sample

```
ID_1 ID_1 -0.644020461494502 0.414762354823591 -3.83142 64 -14.0364 5.51742 0.0714337 5.66263 0.865562 -2.26957 -0.0965859 -2.35497 1.05889 0.195302 0 7
ID_2 ID_2 -0.0278698075809014 0.00077672617459647 0.614044 66 -10.8505 2.11998 -0.882883 -0.441662 -2.64177 2.78944 0.524586 2.67134 -2.63724 -0.998764 1 20
ID_3 ID_3 2.1286574811167 4.53118267191409 -0.237792 55 -9.75369 3.18343 -2.09793 6.87345 11.3777 2.96961 -1.11879 0.873649 3.35523 -4.57831 1 10
ID_4 ID_4 2.1286574811167 4.53118267191409 6.69866 47 -9.07045 0.956878 -2.48407 1.06359 -3.13247 2.1232 -0.00976751 0.820582 0.0305345 1.6303 1 20
ID_5 ID_5 -0.952095788451302 0.906486390386706 -1.61423 59 -12.9379 1.29461 -1.79973 1.44404 -6.82898 -2.96795 -2.91577 -1.82881 7.15892 2.10916 1 20
```

10) Qpt.txt - This is a .txt file which contains the following columns in order. The target dataset has 200 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID  
* quantitative phenotype of the target sample

```
ID_801 ID_801 26.5723
ID_802 ID_802 20.2632
ID_803 ID_803 27.7365
ID_804 ID_804 18.75
ID_805 ID_805 23.3025
```

11) Qct.txt - This is a .txt file which contains the following columns in order. The target dataset has 200 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID 
* standardized covariate 
* square of the standardized covariate (Note: Although squared covariate term is not involved when analyzing quantitative traits, this file can/may be used as the covariate file when modulating with a binary outcome. In that event, this column is required to address model misspecification, if present. See reference for details.) 
* 14 confounders of the target sample

```
ID_801 ID_801 -0.644020461494502 0.414762354823591 -3.82659 69 -13.8514 3.9608 -1.78805 0.0692473 -6.32556 2.85359 1.08516 -1.30304 3.41659 1.41577 0 7
ID_802 ID_802 -0.952095788451302 0.906486390386706 2.06515 60 -12.2438 4.04169 -0.905739 5.9656 8.35545 -1.43576 -0.618153 0.746918 5.11019 -0.207188 1 19
ID_803 ID_803 -0.0278698075809014 0.00077672617459647 -0.795863 62 -10.9195 6.91985 -2.92088 1.26019 -5.56624 -0.552624 -0.0756095 -0.910047 -1.33896 1.72636 0 7
ID_804 ID_804 -0.644020461494502 0.414762354823591 -2.62088 67 -9.9271 4.1096 -2.35454 0.719021 -1.82806 -1.82107 1.21574 -3.56693 -7.91232 2.71011 0 10
ID_805 ID_805 0.280205519375899 0.0785151330887172 -3.33164 67 -11.8637 5.88272 1.07288 2.74488 -7.32776 -2.39477 -3.07983 -1.43625 2.08822 1.42939 1 15
```

# Package installation
The current GitHub version of **GxEprs** can be installed via:
```
library(devtools)
install_github("DoviniJ/GxEprs") 
```
# Load the library
```
library(GxEprs)
```

# Quick start
We will illustrate the usage of **GxEprs** using a few example datasets downloaded and white labelled from UK biobank. Follow the step wise process:

##### Step 1: Download and install plink2 in your machine. The package supports both Linux and Windows version.
Link: https://www.cog-genomics.org/plink/2.0/

##### Step 2: Obtain the path to the executable plink application <plink_path>

##### Step 3.1: Run the following code (functions should be run in the given order) to obtain the risk scores of individuals in the target dataset:

###### Step 3.1.1 Give the path where plink executable file is located
```
plink_path <- "<plink_path>/plink2" 
```
###### Step 3.1.2 It is always recommended to check how the files look like before using them in functions, for better understanding. You may directly use the data files embedded in the package as a trial. Note that, for convenience, we have used identical names for the embedded data object, and for the corresponding function argument. You can check the top proportion of each data file using the following code:
```
head(Bphe_discovery) #phenotype file of the discovery sample when the outcome is binary
head(Bcov_discovery) #covariate file of the discovery sample when the outcome is binary
head(Bphe_target) #phenotype file of the target sample when the outcome is binary
head(Bcov_target) #covariate file of the target sample when the outcome is binary
head(Qphe_discovery) #phenotype file of the discovery sample when the outcome is quantitative
head(Qcov_discovery) #covariate file of the discovery sample when the outcome is quantitative
head(Qphe_target) #phenotype file of the target sample when the outcome is quantitative
head(Qcov_target) #covariate file of the target sample when the outcome is quantitative
```

###### Step 3.1.3 To call the data files saved in "inst" directory, you can follow the following code to obtain the path of each data file. 
```
inst_path <- system.file(package = "GxEprsDummy") 
DummyData <- paste0(inst_path, "/DummyData") #this contains all .fam, .bed and .bim files. They can be accessed by a direct call of prefix "DummyData"
Bphe_discovery <- paste0(inst_path, "/Bphe_discovery.txt")
Bcov_discovery <- paste0(inst_path, "/Bcov_discovery.txt")
Bphe_target <- paste0(inst_path, "/Bphe_target.txt")
Bcov_target <- paste0(inst_path, "/Bcov_target.txt")
Qphe_discovery <- paste0(inst_path, "/Qphe_discovery.txt")
Qcov_discovery <- paste0(inst_path, "/Qcov_discovery.txt")
Qphe_target <- paste0(inst_path, "/Qphe_target.txt")
Qcov_target <- paste0(inst_path, "/Qcov_target.txt")
```
Note that the step 3.1.3 described above is to call the embedded data files in this package itself. However, when users have to call their own data, they can follow the same approach. It is more convenient if the users can store all their data files in the same working directory. For example, assume that the file names are as follows 
(Refer to 'File formats' section of this document to view the formatting details of each of the following input file): 
* binary files: **mydata.fam**, **mydata.bim** and **mydata.bed**
* phenotype file of discovery sample (binary outcome): **Bpd.txt**
* covariate file of discovery sample (binary outcome): **Bcd.txt**
* phenotype file of target sample (binary outcome): **Bpt.txt**
* covariate file of the target sample (binary outcome): **Bct.txt**
* phenotype file of discovery sample (quantitative outcome): **Qpd.txt** 
* covariate file of discovery sample (quantitative outcome): **Qcd.txt** 
* phenotype file of target sample (quantitative outcome): **Qpt.txt**
* covariate file of the target sample (quantitative outcome): **Qct.txt**

###### Additional note:
_Note that, all these files can be placed in a separate location. It is always upto the users choice. In that case remember to give the full path to the file location since R identifies files by name, only when they are in the same directory._
```
b_file <- "<path>/mydata"
Bphe_discovery <- "<path>/Bpd.txt"
Bcov_discovery <- "<path>/Bcd.txt"
Bphe_target <- "<path>/Bpt.txt"
Bcov_target <- "<path>/Bct.txt"
Qphe_discovery <- "<path>/Qpd.txt"
Qcov_discovery <- "<path>/Qcd.txt"
Qphe_target <- "<path>/Qpt.txt"
Qcov_target <- "<path>/Qct.txt"
```

###### Step 3.1.4 Set the working directory and run the following R functions in the given order
```
setwd("<path to working directory>") #set the working directory where you need to save the output files
```
$\color{red}{NOTE:}$ Read **manual.pdf** document for descriptions of arguments passed for each function.

###### When the outcome variable is binary
**Command**
```
x <- GWAS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20)
```
As explained above, “mydata” is the prefix of the PLINK format files, “Bpd.txt” is binary phenotype file of the discovery sample, "Bcd.txt" is the covariate file of the discovery sample, thread indicates the number of CPUs used to run the command which can be optionally specified by the user (default is 20). This command performs GWAS using a logistic regression, and outputs GWAS summary statistics of all additive SNP effects.

**Output**
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
1 768448 SNP_1 G A A N ADD 800 -0.0912494357966473 0.421016 -0.216738 0.828413 .
1 853954 SNP_2 A C C N ADD 800 0.580991407122803 0.266219 2.18239 0.0290808 .
1 880390 SNP_3 C A A N ADD 800 0.446184697385182 0.744786 0.599076 0.549122 .
1 940203 SNP_4 G A A N ADD 800 0.476159642057072 0.456313 1.0435 0.296718 .
```
```x```  contains GWAS summary statistics of all additive SNP effects, when the outcome is binary. V1 to V14 denote the following columns in order. Note that all summary statistics follow the same structure.
* chromosome 
* base pair position 
* SNP ID 
* reference allele 
* alternate allele 
* counted allele A1 (in regression) 
* firth regression status 
* test identifier 
* number of samples in regression 
* odds ratio for A1 allele 
* standard error of log odds 
* test statistic 
* p-value  
* error code 

See the topic GWAS_binary (page 6) in manual.pdf for examples.

**Command**
```
x <- GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20)
```
This performs GWEIS using a logistic regression, and outputs GWEIS summary statistics of both additive and interaction SNP effects.

**Output**
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
1 768448 SNP_1 G A A N ADD 800 -0.0429377728682119 0.546246 -0.0786059 0.937346 .
1 853954 SNP_2 A C C N ADD 800 0.691871367055893 0.338615 2.04325 0.0410281 .
1 880390 SNP_3 C A A N ADD 800 -0.0243828605684827 1.27677 -0.0190971 0.984764 .
1 940203 SNP_4 G A A N ADD 800 0.532238718269307 0.57697 0.922471 0.356283 .
```
```x[[1]]``` contains GWEIS summary statistics of all additive SNP effects, when the outcome is binary. 

```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
1 768448 SNP_1 G A A N ADDxCOVAR1 800 -0.15253365458422 0.37862 -0.402868 0.687045 .
1 853954 SNP_2 A C C N ADDxCOVAR1 800 -0.114150142610688 0.211225 -0.540423 0.588905 .
1 880390 SNP_3 C A A N ADDxCOVAR1 800 0.612983215535033 0.936868 0.654292 0.512924 .
1 940203 SNP_4 G A A N ADDxCOVAR1 800 -0.311823137172646 0.390234 -0.799069 0.424251 .
```
```x[[2]]```  contains GWEIS summary statistics of all interaction SNP effects, when the outcome is binary. 

See the topic GWEIS_binary (page 9) in manual.pdf for examples.

**Command**
```
a <- GWAS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")
b <- GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")

x <- PRS_binary(plink_path, "mydata", summary_input = a)
y <- PRS_binary(plink_path, "mydata", summary_input = b[[1]])
z <- PRS_binary(plink_path, "mydata", summary_input = b[[2]])
```
As explained above, “mydata” is the prefix of the PLINK format files, a, b[[1]] and b[[2]] are summary statistics generated from previous functions and used as an input for this function to construct PRS. These commands compute polygenic risk scores for each individual in the target dataset and outputs the PRSs of all individuals.

**Output**
```
V1 V2 V3 V4 V5
ID_1 ID_1 1970 459 -0.0112131
ID_2 ID_2 1990 535 0.00675556
ID_3 ID_3 1970 432 -0.00961485
ID_4 ID_4 1976 469 -0.00708284
```
```x```  contains the following columns in order.
* FID 
* IID 
* number of alleles across scored variants (ALLELE_CT)
* sum of named allele dosages (NAMED_ALLELE_DOSAGE_SUM)
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the additive effects of GWAS summary statistics), of the full dataset


```
V1 V2 V3 V4 V5
ID_1 ID_1	1970 459 -0.00728841
ID_2 ID_2	1990 535 0.0298431
ID_3 ID_3	1970 432 -0.000156035
ID_4 ID_4	1976 469 -0.00745457
```
```y```  This contains the the following columns in order.
* FID 
* IID 
* ALLELE_CT
* NAMED_ALLELE_DOSAGE_SUM
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the additive effects of GWEIS summary statistics), of the full dataset 


```
V1 V2 V3 V4 V5
ID_1 ID_1	1970 459 -0.00572028
ID_2 ID_2	1990 535 -0.0262972
ID_3 ID_3	1970 432 -0.0126608
ID_4 ID_4 1976 469 -0.00128371
```
```z```  contains the the following columns in order.
* FID 
* IID 
* ALLELE_CT
* NAMED_ALLELE_DOSAGE_SUM
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the interaction effects of GWEIS summary statistics), of the full dataset

See the topic PRS_binary (page 12) in manual.pdf for examples.

**Command**
```
a <- GWAS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")
b <- GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")
p <- PRS_binary(plink_path, "mydata", summary_input = a)
q <- PRS_binary(plink_path, "mydata", summary_input = b[[1]])
r <- PRS_binary(plink_path, "mydata", summary_input = b[[2]])

v <- summary_regular_binary("Bpt.txt", "Bct.txt", trd_score = p, Model = 1)
w <- summary_regular_binary("Bpt.txt", "Bct.txt", add_score = q, Model = 2)
x <- summary_regular_binary("Bpt.txt", "Bct.txt", add_score = q, gxe_score = r, Model = 3)
y <- summary_regular_binary("Bpt.txt", "Bct.txt", add_score = q, gxe_score = r, Model = 4)
z <- summary_regular_binary("Bpt.txt", "Bct.txt", add_score = q, gxe_score = r, Model = 5)
```
“Bpt.txt” is binary phenotype file of the target sample, "Bct.txt" is the covariate file of the target sample, p, q and r are the PRSs generated. Depending on the model used, the input should be varied. (See section $\color{red}{IMPORTANT}$ for the target models.) This function outputs (for demonstration, we used Model = 5 situation) both the summary of the fitted **regular** model for **binary** outcome and all the calculated individual risk scores of the fitted **regular** model for **binary** outcome. 

##### Refer to the section $\color{red}{IMPORTANT}$ at the end of this document for details about models fitted at this step.

**Output**
```
Call:
glm(formula = out ~ ., family = binomial(link = logit), data = df_new)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.2110  -0.3911  -0.2364  -0.1357   2.5362  

Coefficients:
               Estimate Std. Error z value Pr(>|z|)  
(Intercept)    4.058401   3.526630   1.151   0.2498  
E              0.112498   0.424978   0.265   0.7912  
`E squared`   -0.255133   0.629549  -0.405   0.6853  
PRS_add       -1.323038   0.867841  -1.525   0.1274  
PRS_gxe       -1.504617   1.076980  -1.397   0.1624  
`PRS_gxe x E`  0.224745   0.502532   0.447   0.6547  
V7             0.080369   0.101398   0.793   0.4280  
V8            -0.039235   0.036984  -1.061   0.2887  
V9             0.239702   0.238018   1.007   0.3139  
V10            0.107022   0.232814   0.460   0.6457  
V11            0.096563   0.198627   0.486   0.6269  
V12           -0.077362   0.160028  -0.483   0.6288  
V13            0.044817   0.085248   0.526   0.5991  
V14            0.025107   0.206558   0.122   0.9033  
V15            0.227534   0.167280   1.360   0.1738  
V16            0.025793   0.165044   0.156   0.8758  
V17            0.003261   0.081670   0.040   0.9681  
V18            0.082463   0.159361   0.517   0.6048  
V19            1.027611   0.697780   1.473   0.1408  
V20           -0.189688   0.079375  -2.390   0.0169 *
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 106.55  on 199  degrees of freedom
Residual deviance:  85.53  on 180  degrees of freedom
AIC: 125.53

Number of Fisher Scoring iterations: 7
```
```z[[1]][[1]]``` contains the target regular model summary output, when the outcome is binary. 

```
ID_802 ID_802 0.142584993327248
ID_803 ID_803 0.00944707783426451
ID_804 ID_804 0.314374300730685
ID_805 ID_805 0.0401077327720606
ID_806 ID_806 0.0272302426618824
```
```z[[2]]``` contains all the calculated individual risk scores using the target dataset (e.g. Model 5), when the outcome is binary. The columns denote the following in order.
* FID
* IID 
* estimated risk value

Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models (using ```summary_regular_binary("Bpt.txt", "Bct.txt", add_score = q, gxe_score = r, Model = 5)``` and ```summary_permuted_binary("Bpt.txt", "Bct.txt", iterations = 1000, add_score = q, gxe_score = r)```), if you choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 5) when generating risk scores. If the 'PRS_gxe x E' term is significant in Model 5, and insignificant in Model 5* (permuted p value), consider that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model). 

See the topic summary_regular_binary (page 20) in manual.pdf for examples.

**Command**
```
a <- GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")
p <- PRS_binary(plink_path, "mydata", summary_input = a[[1]])
q <- PRS_binary(plink_path, "mydata", summary_input = a[[2]])
summary_permuted_binary("Bpt.txt", "Bct.txt", iterations = 1000, add_score = p, gxe_score = q)
```
This outputs the p value of the fitted **permuted** model for **binary** outcome.




###### When the outcome variable is quantitative
**Command**
```
x <- GWAS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20)
```
As explained above, “mydata” is the prefix of the PLINK format files, “Qpd.txt” is quantitative phenotype file of the discovery sample, "Qcd.txt" is the covariate file of the discovery sample, thread indicates the number of CPUs used to run the command which can be optionally specified by the user (default is 20). This command performs GWAS using a linear regression, and outputs GWAS summary statistics of all additive SNP effects.

**Output**
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13
1 768448 SNP_1 G A A ADD 800 0.183597 0.35942 0.510815 0.609624 .
1 853954 SNP_2 A C C ADD 800 -0.336911 0.234538 -1.43649 0.151262 .
1 880390 SNP_3 C A A ADD 800 0.139377 0.683708 0.203855 0.83852 .
1 940203 SNP_4 G A A ADD 800 -0.896864 0.460195 -1.94888 0.0516666 .
```
```x``` contains GWAS summary statistics of all additive SNP effects, when the outcome is quantitative. V1 to V13 denote the following columns in order. Note that all summary statistics follow the same structure.
* chromosome 
* base pair position 
* SNP ID 
* reference allele 
* alternate allele 
* counted allele A1 (in regression) 
* test identifier 
* number of samples in regression 
* odds ratio for A1 allele 
* standard error of log odds 
* test statistic 
* p-value  
* error code 

See the topic GWAS_quantitative (page 7) in manual.pdf for examples.

**Command**
```
x <- GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20)
```
This performs GWEIS using a linear regression, and outputs GWEIS summary statistics of all additive and interaction SNP effects. 

**Output**
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13
1 768448 SNP_1 G A A ADD 800 0.170496 0.360579 0.472841 0.636459 .
1 853954 SNP_2 A C C ADD 800 -0.334001 0.235261 -1.4197 0.156093 .
1 880390 SNP_3 C A A ADD 800 0.255825 0.743858 0.343916 0.731002 .
1 940203 SNP_4 G A A ADD 800 -0.954298 0.471745 -2.02291 0.043422 .
```
```x[[1]]``` contains GWEIS summary statistics of all additive SNP effects, when the outcome is quantitative. 

```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13
1 768448 SNP_1 G A A ADDxCOVAR1 800 0.349176 0.370182 0.943256 0.345841 .
1 853954 SNP_2 A C C ADDxCOVAR1 800 0.065845 0.229713 0.28664 0.774464 .
1 880390 SNP_3 C A A ADDxCOVAR1 800 0.364167 0.883745 0.412073 0.680399 .
1 940203 SNP_4 G A A ADDxCOVAR1 800 -0.262215 0.481748 -0.544299 0.586391 .
```
```x[[2]]``` contains GWEIS summary statistics of all interaction SNP effects, when the outcome is quantitative. 

See the topic GWEIS_quantitative (page 10) in manual.pdf for examples.

**Command**
```
a <- GWAS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt")
b <- GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt")

x <- PRS_quantitative(plink_path, "mydata", summary_input = a)
y <- PRS_quantitative(plink_path, "mydata", summary_input = b[[1]])
z <- PRS_quantitative(plink_path, "mydata", summary_input = b[[2]])
```
As explained above, “mydata” is the prefix of the PLINK format files, a, b[[1]] and b[[2]] are summary statistics generated from previous functions and used as an input for this function to construct PRS. These commands compute polygenic risk scores for each individual in the target dataset and outputs the PRSs of all individuals.

**Output**
```
V1 V2 V3 V4 V5
ID_1 ID_1 1970 459 -0.00436674
ID_2 ID_2 1990 535 0.00120973
ID_3 ID_3 1970 432 0.00158428
ID_4 ID_4 1976 469 -0.000431784
```
```x``` contains the following columns in order.
* FID 
* IID 
* number of alleles across scored variants (ALLELE_CT)
* sum of named allele dosages (NAMED_ALLELE_DOSAGE_SUM)
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the additive effects of GWAS summary statistics), of the full dataset

```
V1 V2 V3 V4 V5
ID_1 ID_1	1970 459 -0.00421272
ID_2 ID_2	1990 535 0.00137912
ID_3 ID_3	1970 432 0.00160997
ID_4 ID_4	1976 469 -0.000112026
```
```y``` contains the the following columns in order.
* FID 
* IID 
* ALLELE_CT
* NAMED_ALLELE_DOSAGE_SUM
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the additive effects of GWEIS summary statistics), of the full dataset 

```
V1 V2 V3 V4 V5
ID_1 ID_1	1970 459 0.00180745
ID_2 ID_2	1990 535 -0.00172071
ID_3 ID_3	1970 432 0.000248408
ID_4 ID_4	1976 469 -0.00144963
```
```z``` contains the the following columns in order.
* FID 
* IID 
* ALLELE_CT
* NAMED_ALLELE_DOSAGE_SUM
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the interaction effects of GWEIS summary statistics), of the full dataset


See the topic PRS_quantitative (page 13) in manual.pdf for examples.

**Command**
```
a <- GWAS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt")
b <- GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt")
p <- PRS_quantitative(plink_path, "mydata", summary_input = a)
q <- PRS_quantitative(plink_path, "mydata", summary_input = b[[1]])
r <- PRS_quantitative(plink_path, "mydata", summary_input = b[[2]])

v <- summary_regular_quantitative("Qpt.txt", "Qct.txt", trd_score = p, Model = 1)
w <- summary_regular_quantitative(""Qpt.txt", "Qct.txt", add_score = q, Model = 2)
x <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = q, gxe_score = r, Model = 3)
y <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = q, gxe_score = r, Model = 4)
```
“Qpt.txt” is quantitative phenotype file of the target sample, "Qct.txt" is the covariate file of the target sample, p, q and r are the PRSs generated. Depending on the model used, the input should be varied. (See section $\color{red}{IMPORTANT}$ for the target models.) This function outputs (for demonstration, we used Model = 4 situation) both the summary of the fitted **regular** model for **quantitative** outcome and all the calculated individual risk scores of the fitted **regular** model for **quantitative** outcome. 

##### Refer to the section $\color{red}{IMPORTANT}$ at the end of this document for details about models fitted at this step.

**Output**
```
Call:
lm(formula = out ~ ., data = df_new)

Residuals:
    Min      1Q  Median      3Q     Max 
-2.1824 -0.6360 -0.0455  0.4398  3.9346 

Coefficients:
                Estimate Std. Error t value Pr(>|t|)   
(Intercept)    0.3965247  0.7666961   0.517  0.60566   
E             -0.0157672  0.0744757  -0.212  0.83257   
PRS_add       -0.0488255  0.0711937  -0.686  0.49371   
PRS_gxe        0.0746073  0.0725750   1.028  0.30532   
`PRS_gxe x E`  0.0009011  0.0786609   0.011  0.99087   
V6            -0.0679999  0.0277280  -2.452  0.01514 * 
V7            -0.0053487  0.0086598  -0.618  0.53758   
V8             0.0423173  0.0491093   0.862  0.38999   
V9             0.0882473  0.0478252   1.845  0.06664 . 
V10           -0.0830620  0.0465820  -1.783  0.07624 . 
V11           -0.0293785  0.0331433  -0.886  0.37657   
V12            0.0205769  0.0147545   1.395  0.16484   
V13           -0.0521727  0.0481595  -1.083  0.28010   
V14           -0.0619072  0.0451513  -1.371  0.17204   
V15           -0.0441643  0.0389190  -1.135  0.25797   
V16           -0.0035800  0.0182968  -0.196  0.84509   
V17            0.0133222  0.0394064   0.338  0.73570   
V18            0.4620297  0.1431690   3.227  0.00148 **
V19           -0.0239217  0.0145944  -1.639  0.10293   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.9697 on 181 degrees of freedom
Multiple R-squared:  0.1447,	Adjusted R-squared:  0.05965 
F-statistic: 1.701 on 18 and 181 DF,  p-value: 0.04244
```
```y[[1]][[1]]``` contains the target regular model summary output, when the outcome is quantitative. 

```
ID_802 ID_802 -0.180057820294863
ID_803 ID_803 0.402506322498664
ID_804 ID_804 0.424579119997805
ID_805 ID_805 0.583828103567573
ID_806 ID_806 0.313724826923073
```
```y[[2]]``` contains all the calculated individual risk scores using the target dataset (e.g. Model 4), when the outcome is quantitative. The columns denote the following in order.
* FID
* IID 
* estimated risk value

Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models (using ```summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = q, gxe_score = r, Model = 4)``` and ```summary_permuted_quantitative("qpt.txt", "qct.txt", iterations = 1000, add_score = q, gxe_score = r)```), if you choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 4) when generating risk scores. If the 'PRS_gxe x E' term is significant in Model 4, and insignificant in Model 4* (permuted p value), consider that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model). 

See the topic summary_regular_quantitative (page 22) in manual.pdf for examples.

**Command**
```
a <- GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt")
p <- PRS_quantitative(plink_path, "mydata", summary_input = a[[1]])
q <- PRS_quantitative(plink_path, "mydata", summary_input = a[[2]])
summary_permuted_quantitative("Qpt.txt", "Qct.txt", iterations = 1000, add_score = p, gxe_score = q)
```
This outputs the p value of the fitted **permuted** model for **quantitative** outcome.


## $$\color{red}{IMPORTANT}$$
The discovery model used in ```GWAS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20, summary_output = "B_trd.sum")``` or ```GWAS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20)``` is as follows:
* y = b_trd.W + error
 where y is the outcome variable, b_trd is the estimated SNP effect and W is the SNP genotype.
 
The discovery model used in ```GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20)``` is as follows:
 * y = b_add.W + b_cov.E + b_cov2.E^2 + b_gxe.(WxE) + error
where y is the outcome variable, b_add is the estimated additive SNP effect, E is the covariate, W is the SNP genotype, b_cov is the estimated effect of the covariate, b_cov2 is the estimated effect of the squared covariate and b_gxe is the estimated effect of the Hadamard product of WxE.

The discovery model used in ```GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20)``` is as follows:
 * y = b_add.W + b_cov.E + b_gxe.(WxE) + error
where y is the outcome variable, b_add is the estimated additive SNP effect, E is the covariate, W is the SNP genotype, b_cov is the estimated effect of the covariate and b_gxe is the estimated effect of the Hadamard product of WxE.


The fitted (target) models in ```summary_regular_binary("Bpt.txt", "Bct.txt", trd_score = NULL, add_score = NULL, gxe_score = NULL, Model)``` or ```summary_regular_quantitative("Qpt.txt", "Qct.txt", trd_score = NULL, add_score = NULL, gxe_score = NULL, Model)``` are as follows:

* Model 1: y = PRS_trd + E + PRS_trd x E + confounders + error
* Model 2: y = PRS_add + E + PRS_add x E + confounders + error
* Model 3: y = PRS_add + E + PRS_gxe x E + confounders + error
* Model 4: y = PRS_add + E + PRS_gxe + PRS_gxe x E + confounders + error
* Model 4*: permuted Model 4
* Model 5: y = PRS_add + E + E^2 + PRS_gxe + PRS_gxe x E + confounders + error
* Model 5*: permuted Model 5

where y is the outcome variable, E is the covariate of interest, PRS_trd and PRS_add are the polygenic risk scores computed using additive SNP effects of GWAS and GWEIS summary statistics respectively, and PRS_gxe is the polygenic risk scores computed using GxE interaction SNP effects of GWEIS summary statistics.

When deciding on the number of permutations to be used, always get an idea from the p value obtained from either Model 4 or 5 (accordingly). If that p value is insignificant, you can use any number of permututations (e.g. 1000), but if that p value is highly significant (p value is very small and very close to zero), it is recommended to increase the number of permutations for better results. For example assume that the p value obtained from Model 4 is 0.000154, then it is recommended to use at least 10,000 iterations in the permutation model.


# Contact 
Please contact Dovini Jayasinghe (dovini.jayasinghe@mymail.unisa.edu.au) or Md Moksedul Momin (md_moksedul.momin@mymail.unisa.edu.au) or Hong Lee (hong.lee@unisa.edu.au) for queries.
