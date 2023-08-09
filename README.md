---
Title: "GxEprs"
Authors: "Dovini Jayasinghe, Md Moksedul Momin and Hong Lee"
Last updated: "09-08-2023"
---

# GxEprs
The 'GxEprs' is an R package to detect and estimate GxE. It uses a novel PRS model that can enhance the prediction accuracy by utilising GxE effects. Firstly it performs Genome Wide Association Studies (GWAS)  and Genome Wide Environment Interaction Studies (GWEIS) using the discovery dataset (see functions ```GWAS_binary()```,```GWAS_quantitative()```, ```GWEIS_binary()```, ```GWEIS_quantitative()```). See the section $\color{red}{IMPORTANT}$ for the discovery models used. Secondly, it uses the GWAS and GWEIS summary statistics generated from the fucntions above to obtain polygenic risk scores (PRSs) (see functions ```PRS_binary()``` and ```PRS_quantitative()```) for the target sample. Finally it predicts the risk values of each individual in the target sample (see functions ```summary_regular_binary()``` and ```summary_regular_quantitative()```). Note that the users can fit 4 different models when the outcome is a quantitative trait, and 5 different models when the outcome is a binary disease trait. See the section $\color{red} {IMPORTANT}$ for the target models used. Finally, it is recommended to check the p-value from permutations using Model 4 (see function ```summary_permuted_quantitative()```), and Model 5 (see function```summary_permuted_binary()```), to make sure that the significance of GxE is not spurious due to model misspecification (Jayasinghe et. al.).

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
###### Step 3.1.2 It is always recommended to check how the files look like before using them in functions, for better understanding. You may directly use the data files embedded in the package as a trial. You can check the top proportion of each data file using the following code:
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
inst_path <- system.file(package = "GxEprs") 
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
* covariate file of target sample (binary outcome): **Bct.txt**
* phenotype file of discovery sample (quantitative outcome): **Qpd.txt** 
* covariate file of discovery sample (quantitative outcome): **Qcd.txt** 
* phenotype file of target sample (quantitative outcome): **Qpt.txt**
* covariate file of target sample (quantitative outcome): **Qct.txt**

###### Additional note:
_Note that, all these files (user's own data files) can be placed in a separate location. It is always upto the user's choice. In that case remember to give the full path to the file location since R identifies files by name, only when they are in the same directory. For example, if your data files are in a different directory;_

* instead of ```"mydata"``` use ```"<path>/mydata"``` as DummyData
* instead of ```"Bpd.txt"``` use ```"<path>/Bpd.txt"``` as Bphe_discovery
* instead of ```"Bcd.txt"``` use ```"<path>/Bcd.txt"``` as Bcov_discovery
* instead of ```"Bpt.txt"``` use ```"<path>/Bpt.txt"``` as Bphe_target
* instead of ```"Bct.txt"``` use ```"<path>/Bct.txt"``` as Bcov_target
* instead of ```"Qpd.txt"``` use ```"<path>/Qpd.txt"``` as Qphe_discovery
* instead of ```"Qcd.txt"``` use ```"<path>/Qcd.txt"``` as Qcov_discovery
* instead of ```"Qpt.txt"``` use ```"<path>/Qpt.txt"``` as Qphe_target
* instead of ```"Qct.txt"``` use ```"<path>/Qct.txt"``` as Qcov_target


###### Step 3.1.4 Set the working directory and run the following R functions
```
setwd("<path to working directory>") #set the working directory
```
$\color{red}{NOTE:}$ Read **manual.pdf** document for descriptions of arguments passed for each function.

###### When the outcome variable is binary
**Command**
```
x <- GWAS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20)
```
As explained above, “mydata” is the prefix of the PLINK format files, “Bpd.txt” is binary phenotype file of the discovery sample, "Bcd.txt" is the covariate file of the discovery sample, thread indicates the number of CPUs used to run the command which can be optionally specified by the user (default is 20). This command performs GWAS using a logistic regression, and outputs GWAS summary statistics of all additive SNP effects.

To perform the same on embedded data
```
x <- GWAS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery, thread = 20)
```

**Output**
```
  CHROM POS ID REF ALT A1 OBS_CT OR LOG_OR_SE Z_STAT P
1   1 768448 SNP_1 G A A 800 -0.0912494357966473 0.421016 -0.216738 0.828413
2   1 853954 SNP_2 A C C 800 0.580991407122803 0.266219 2.18239 0.0290808
3   1 880390 SNP_3 C A A 800 0.446184697385182 0.744786 0.599076 0.549122
4   1 940203 SNP_4 G A A 800 0.476159642057072 0.456313 1.0435 0.296718
```
```x```  contains GWAS summary statistics of all additive SNP effects, when the outcome is binary. 

* x$CHROM : chromosome number 
* x$POS : base pair position
* x$ID : SNP ID
* x$REF : reference allele
* x$ALT : alternate allele 
* x$A1 : minor allele
* x$OBS_CT : number of allele observations 
* x$OR : odds ratios of SNP effects
* x$LOG_OR_SE : standard errors of log odds of SNP effects
* x$Z_STAT : test statistics of SNP effects
* x$P : p values of SNP effects

See the topic GWAS_binary (page 6) in manual.pdf for examples.

**Command**
```
x <- GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20)
```
This performs GWEIS using a logistic regression, and outputs GWEIS summary statistics of both additive and interaction SNP effects.


To perform the same on embedded data
```
x <- GWEIS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery, thread = 20)
```

**Output**
```
  CHROM POS ID REF ALT A1 OBS_CT ADD_OR ADD_LOG_OR_SE ADD_Z_STAT ADD_P INTERACTION_OR INTERACTION_LOG_OR_SE INTERACTION_Z_STAT INTERACTION_P
1   1 768448 SNP_1 G A A 800 -0.0429377728682119 0.546246 -0.0786059 0.937346 -0.15253365458422 0.37862 -0.402868 0.687045
2   1 853954 SNP_2 A C C 800 0.691871367055893 0.338615 2.04325 0.0410281 -0.114150142610688 0.211225 -0.540423 0.588905
3   1 880390 SNP_3 C A A 800 -0.0243828605684827 1.27677 -0.0190971 0.984764 0.612983215535033 0.936868 0.654292 0.512924
4   1 940203 SNP_4 G A A 800 0.532238718269307 0.57697 0.922471 0.356283 -0.311823137172646 0.390234 -0.799069 0.424251
```
```x``` contains GWEIS summary statistics of all additive and interaction SNP effects, when the outcome is binary. 

* x$CHROM : chromosome number 
* x$POS : base pair position
* x$ID : SNP ID
* x$REF : reference allele
* x$ALT : alternate allele 
* x$A1 : minor allele
* x$OBS_CT : number of allele observations 
* x$ADD_OR : odds ratios of additive SNP effects
* x$ADD_LOG_OR_SE : standard errors of log odds of additive SNP effects
* x$ADD_Z_STAT : test statistics of additive SNP effects
* x$ADD_P : p values of additive SNP effects
* x$INTERACTION_OR : odds ratios of the SNP effect of interaction SNP effects
* x$INTERACTION_LOG_OR_SE : standard errors of log odds of interaction SNP effects
* x$INTERACTION_Z_STAT : test statistics of interaction SNP effects
* x$INTERACTION_P : p values of interaction SNP effects

See the topic GWEIS_binary (page 9) in manual.pdf for examples.

**Command**
```
a <- GWAS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")
trd <- a[c("ID", "A1", "OR")]
b <- GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")
add <- b[c("ID", "A1", "ADD_OR")]
gxe <- b[c("ID", "A1", "INTERACTION_OR")]

x <- PRS_binary(plink_path, "mydata", summary_input = trd)
y <- PRS_binary(plink_path, "mydata", summary_input = add)
z <- PRS_binary(plink_path, "mydata", summary_input = gxe)
```
As explained above, “mydata” is the prefix of the PLINK format files, , add and gxe are summary statistics generated from previous functions and used as an input for this function to construct PRS. These commands compute polygenic risk scores for each individual in the target dataset and outputs the PRSs of all individuals.


To perform the same on embedded data
```
a <- GWAS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery)
trd <- a[c("ID", "A1", "OR")]
b <- GWEIS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery)
add <- b[c("ID", "A1", "ADD_OR")]
gxe <- b[c("ID", "A1", "INTERACTION_OR")]

x <- PRS_binary(plink_path, DummyData, summary_input = trd)
y <- PRS_binary(plink_path, DummyData, summary_input = add)
z <- PRS_binary(plink_path, DummyData, summary_input = gxe)
```

**Output**
```
   FID  IID         PRS
1 ID_1 ID_1 -0.01121310
2 ID_2 ID_2  0.00675556
3 ID_3 ID_3 -0.00961485
4 ID_4 ID_4 -0.00708284
```
```x```  contains the following columns in order.
* x$FID : family IDs of the full dataset
* x$IID : individual IDs of the full dataset
* x$PRS : polygenic risk scores (PRSs), computed from the additive effects of GWAS summary statistics, of the full dataset


```
   FID  IID          PRS
1 ID_1 ID_1 -0.007288410
2 ID_2 ID_2  0.029843100
3 ID_3 ID_3 -0.000156035
4 ID_4 ID_4 -0.007454570
```
```y``` contains the the following columns in order.
* y$FID : family IDs of the full dataset
* y$IID : individual IDs of the full dataset
* y$PRS : polygenic risk scores (PRSs), computed from the additive effects of GWEIS summary statistics), of the full dataset 


```
   FID  IID         PRS
1 ID_1 ID_1 -0.00572028
2 ID_2 ID_2 -0.02629720
3 ID_3 ID_3 -0.01266080
4 ID_4 ID_4 -0.00128371
```
```z```  contains the the following columns in order.
* z$FID : family IDs of the full dataset
* z$IID : individual IDs of the full dataset
* z$PRS : polygenic risk scores (PRSs), computed from the interaction effects of GWEIS summary statistics), of the full dataset

See the topic PRS_binary (page 12) in manual.pdf for examples.

**Command**
```
a <- GWAS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")
trd <- a[c("ID", "A1", "OR")]
b <- GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")
add <- b[c("ID", "A1", "ADD_OR")]
gxe <- b[c("ID", "A1", "INTERACTION_OR")]
p <- PRS_binary(plink_path, "mydata", summary_input = trd)
q <- PRS_binary(plink_path, "mydata", summary_input = add)
r <- PRS_binary(plink_path, "mydata", summary_input = gxe)

u <- summary_regular_binary("Bpt.txt", "Bct.txt", add_score = p, Model = 0)
v <- summary_regular_binary("Bpt.txt", "Bct.txt", add_score = p, Model = 1)
w <- summary_regular_binary("Bpt.txt", "Bct.txt", add_score = q, Model = 2)
x <- summary_regular_binary("Bpt.txt", "Bct.txt", add_score = q, gxe_score = r, Model = 3)
y <- summary_regular_binary("Bpt.txt", "Bct.txt", add_score = q, gxe_score = r, Model = 4)
z <- summary_regular_binary("Bpt.txt", "Bct.txt", add_score = q, gxe_score = r, Model = 5)
```
“Bpt.txt” is binary phenotype file of the target sample, "Bct.txt" is the covariate file of the target sample, p, q and r are the PRSs generated. If users choose to use their own PRSs stored previously, then they should create their own files (having column names as FID, IID and PRS) to use as inputs to ```summary_regular_binary()``` function (see examples.R script file attached to this repository). Depending on the model used, the input should be varied. (See section $\color{red}{IMPORTANT}$ for the target models.) This function outputs both the summary of the fitted **regular** model for **binary** outcome and all the calculated individual risk scores of the fitted **regular** model for **binary** outcome. 

##### Refer to the section $\color{red}{IMPORTANT}$ at the end of this document for details about models fitted at this step.


To perform the same on embedded data:
```
a <- GWAS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery)
trd <- a[c("ID", "A1", "OR")]
b <- GWEIS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery)
add <- b[c("ID", "A1", "ADD_OR")]
gxe <- b[c("ID", "A1", "INTERACTION_OR")]
p <- PRS_binary(plink_path, DummyData, summary_input = trd)
q <- PRS_binary(plink_path, DummyData, summary_input = add)
r <- PRS_binary(plink_path, DummyData, summary_input = gxe)

u <- summary_regular_binary(Bphe_target, Bcov_target, add_score = p, Model = 0)
v <- summary_regular_binary(Bphe_target, Bcov_target, add_score = p, Model = 1)
w <- summary_regular_binary(Bphe_target, Bcov_target, add_score = q, Model = 2)
x <- summary_regular_binary(Bphe_target, Bcov_target, add_score = q, gxe_score = r, Model = 3)
y <- summary_regular_binary(Bphe_target, Bcov_target, add_score = q, gxe_score = r, Model = 4)
z <- summary_regular_binary(Bphe_target, Bcov_target, add_score = q, gxe_score = r, Model = 5)
```
Note: For demonstration, we used Model = 5:

**Output**
```
            Coefficient Std.Error Test.Statistic    pvalue
E            0.10324092 0.4322434     0.23884906 0.8112226
E squared    0.04255961 0.4434101     0.09598248 0.9235345
PRS_add     -1.32972343 0.9151400    -1.45302727 0.1462162
PRS_gxe     -1.79547662 1.1971366    -1.49980932 0.1336638
PRS_gxe x E  0.31121844 0.5955648     0.52256013 0.6012804

```
```z$summary``` contains the target regular model summary output, when the outcome is binary. 

* z$summary[,"Coefficient"] : regression coefficients of all model components
* z$summary[,"Std.Error"] : standard errors of all regression coefficients
* z$summary[,"Test.Statistic"] : test statistic value of all regression coefficients
* z$summary[,"pvalue"] : p-value of all regression coefficients

```
  FID      IID      Risk.Values
1 "ID_801" "ID_801" "0.328602143945709"
2 "ID_802" "ID_802" "0.178590231862973"
3 "ID_803" "ID_803" "0.00732537258032875"
4 "ID_804" "ID_804" "0.285956066373666"
```
```z$risk.values``` contains all the calculated individual risk scores using the target dataset (e.g. Model 5), when the outcome is binary. The columns denote the following in order.
* z$risk.values[,"FID"] : family IDs of the target dataset
* z$risk.values[,"IID"] : individual IDs of the target dataset
* z$risk.values[,"Risk.Values] : estimated risk values of the target dataset

Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models (using ```summary_regular_binary("Bpt.txt", "Bct.txt", add_score = q, gxe_score = r, Model = 5)``` and ```summary_permuted_binary("Bpt.txt", "Bct.txt", iterations = 1000, add_score = q, gxe_score = r)```), if you choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 5) when generating risk scores. If the 'PRS_gxe x E' term is significant in Model 5, and insignificant in Model 5* (permuted p value), consider that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model). See Figure 5 in Jayasinghe et.al. for simulation verification.

See the topic summary_regular_binary (page 19) in manual.pdf for examples.

**Command**
```
a <- GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt")
add <- a[c("ID", "A1", "ADD_OR")]
gxe <- a[c("ID", "A1", "INTERACTION_OR")]
p <- PRS_binary(plink_path, "mydata", summary_input = add)
q <- PRS_binary(plink_path, "mydata", summary_input = gxe)
x <- summary_permuted_binary("Bpt.txt", "Bct.txt", iterations = 1000, add_score = p, gxe_score = q)
```

To perform the same on embedded data:
```
a <- GWEIS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery)
add <- a[c("ID", "A1", "ADD_OR")]
gxe <- a[c("ID", "A1", "INTERACTION_OR")]
p <- PRS_binary(plink_path, DummyData, summary_input = add)
q <- PRS_binary(plink_path, DummyData, summary_input = gxe)
x <- summary_permuted_binary(Bphe_target, Bcov_target, iterations = 1000, add_score = p, gxe_score = q)
```

```x``` outputs the p value of the fitted **permuted** model for **binary** outcome.




###### When the outcome variable is quantitative
**Command**
```
x <- GWAS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20)
```
As explained above, “mydata” is the prefix of the PLINK format files, “Qpd.txt” is quantitative phenotype file of the discovery sample, "Qcd.txt" is the covariate file of the discovery sample, thread indicates the number of CPUs used to run the command which can be optionally specified by the user (default is 20). This command performs GWAS using a linear regression, and outputs GWAS summary statistics of all additive SNP effects.

To perform the same on embedded data:
```
x <- GWAS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery, thread = 20)
```

**Output**
```
   CHROM     POS    ID REF ALT A1 OBS_CT       BETA       SE     T_STAT         P
1     1  768448 SNP_1   G   A  A    800  0.1835970 0.359420  0.5108150 0.6096240
2     1  853954 SNP_2   A   C  C    800 -0.3369110 0.234538 -1.4364900 0.1512620
3     1  880390 SNP_3   C   A  A    800  0.1393770 0.683708  0.2038550 0.8385200
4     1  940203 SNP_4   G   A  A    800 -0.8968640 0.460195 -1.9488800 0.0516666
```
```x``` contains GWAS summary statistics of all additive SNP effects, when the outcome is quantitative. V1 to V13 denote the following columns in order. Note that all summary statistics follow the same structure.
* x$CHROM : chromosome number 
* x$POS : base pair position 
* x$ID : SNP ID 
* x$REF : reference allele 
* x$ALT : alternate allele 
* x$A1 : minor allele 
* x$OBS_CT : number of allele observations 
* x$BETA : SNP effects
* x$SE : standard errors of SNP effects
* x$T_STAT : test statistics of SNP effects
* x$P : p values of SNP effects


See the topic GWAS_quantitative (page 7) in manual.pdf for examples.

**Command**
```
x <- GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20)
```
This performs GWEIS using a linear regression, and outputs GWEIS summary statistics of all additive and interaction SNP effects. 

To perform the same on embedded data
```
x <- GWEIS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery, thread = 20)
```

**Output**
```
  CHROM POS ID REF ALT A1 OBS_CT ADD_BETA ADD_SE ADD_T_STAT ADD_P INTERACTION_BETA INTERACTION_SE INTERACTION_T_STAT INTERACTION_P
1   1 768448 SNP_1 G A A 800 0.170496 0.360579 0.472841 0.636459 0.349176 0.370182 0.943256 0.345841
2   1 853954 SNP_2 A C C 800 -0.334001 0.235261 -1.4197 0.156093 0.065845 0.229713 0.28664 0.774464
3   1 880390 SNP_3 C A A 800 0.255825 0.743858 0.343916 0.731002 0.364167 0.883745 0.412073 0.680399
4   1 940203 SNP_4 G A A 800 -0.954298 0.471745 -2.02291 0.043422 -0.262215 0.481748 -0.544299 0.586391
```
```x``` contains GWEIS summary statistics of all additive and interaction SNP effects, when the outcome is quantitative. 


See the topic GWEIS_quantitative (page 10) in manual.pdf for examples.

**Command**
```
a <- GWAS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt")
trd <- a[c("ID", "A1", "OR")]
b <- GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt")
add <- b[c("ID", "A1", "ADD_OR")]
gxe <- b[c("ID", "A1", "INTERACTION_OR")]

x <- PRS_quantitative(plink_path, "mydata", summary_input = trd)
y <- PRS_quantitative(plink_path, "mydata", summary_input = add)
z <- PRS_quantitative(plink_path, "mydata", summary_input = gxe)
```
As explained above, “mydata” is the prefix of the PLINK format files, , add and gxe are summary statistics generated from previous functions and used as an input for this function to construct PRS. These commands compute polygenic risk scores for each individual in the target dataset and outputs the PRSs of all individuals.


To perform the same on embedded data
```
a <- GWAS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery)
trd <- a[c("ID", "A1", "OR")]
b <- GWEIS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery)
add <- b[c("ID", "A1", "ADD_OR")]
gxe <- b[c("ID", "A1", "INTERACTION_OR")]

x <- PRS_quantitative(plink_path, DummyData, summary_input = trd)
y <- PRS_quantitative(plink_path, DummyData, summary_input = add)
z <- PRS_quantitative(plink_path, DummyData, summary_input = gxe)
```

**Output**
```
   FID  IID          PRS
1 ID_1 ID_1 -0.004366740
2 ID_2 ID_2  0.001209730
3 ID_3 ID_3  0.001584280
4 ID_4 ID_4 -0.000431784
```
```x```  contains the following columns in order.
* x$FID : family IDs of the full dataset
* x$IID : individual IDs of the full dataset
* x$PRS : polygenic risk scores (PRSs), computed from the additive effects of GWAS summary statistics, of the full dataset


```
   FID  IID          PRS
1 ID_1 ID_1 -0.007288410
2 ID_2 ID_2  0.029843100
3 ID_3 ID_3 -0.000156035
4 ID_4 ID_4 -0.007454570
```
```y``` contains the the following columns in order.
* y$FID : family IDs of the full dataset
* y$IID : individual IDs of the full dataset
* y$PRS : polygenic risk scores (PRSs), computed from the additive effects of GWEIS summary statistics), of the full dataset 


```
   FID  IID          PRS
1 ID_1 ID_1  1.80745e-03
2 ID_2 ID_2 -1.72071e-03
3 ID_3 ID_3  2.48408e-04
4 ID_4 ID_4 -1.44963e-03
```
```z```  contains the the following columns in order.
* z$FID : family IDs of the full dataset
* z$IID : individual IDs of the full dataset
* z$PRS : polygenic risk scores (PRSs), computed from the interaction effects of GWEIS summary statistics), of the full dataset


See the topic PRS_quantitative (page 13) in manual.pdf for examples.

**Command**
```
a <- GWAS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt")
trd <- a[c("ID", "A1", "BETA")]
b <- GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt")
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]
p <- PRS_quantitative(plink_path, "mydata", summary_input = trd)
q <- PRS_quantitative(plink_path, "mydata", summary_input = add)
r <- PRS_quantitative(plink_path, "mydata", summary_input = gxe)

u <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = p, Model = 0)
v <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = p, Model = 1)
w <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = q, Model = 2)
x <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = q, gxe_score = r, Model = 3)
y <- summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = q, gxe_score = r, Model = 4)
```
“Qpt.txt” is quantitative phenotype file of the target sample, "Qct.txt" is the covariate file of the target sample, p, q and r are the PRSs generated. If users choose to use their own PRSs stored previously, then they should create their own files (having column names as FID, IID and PRS) to use as inputs to ```summary_regular_quantitative()``` function (see examples.R script file attached to this repository). Depending on the model used, the input should be varied. (See section $\color{red}{IMPORTANT}$ for the target models.) This function outputs both the summary of the fitted **regular** model for **quantitative** outcome and all the calculated individual risk scores of the fitted **regular** model for **quantitative** outcome. 

##### Refer to the section $\color{red}{IMPORTANT}$ at the end of this document for details about models fitted at this step.


To perform the same on embedded data:
```
a <- GWAS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery)
trd <- a[c("ID", "A1", "BETA")]
b <- GWEIS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery)
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]
p <- PRS_quantitative(plink_path, DummyData, summary_input = trd)
q <- PRS_quantitative(plink_path, DummyData, summary_input = add)
r <- PRS_quantitative(plink_path, DummyData, summary_input = gxe)

u <- summary_regular_quantitative(Qphe_target, Qcov_target, add_score = p, Model = 0)
v <- summary_regular_quantitative(Qphe_target, Qcov_target, add_score = p, Model = 1)
w <- summary_regular_quantitative(Qphe_target, Qcov_target, add_score = q, Model = 2)
x <- summary_regular_quantitative(Qphe_target, Qcov_target, add_score = q, gxe_score = r, Model = 3)
y <- summary_regular_quantitative(Qphe_target, Qcov_target, add_score = q, gxe_score = r, Model = 4)
```
Note: For demonstration, we used Model = 4:

**Output**
```
            Coefficient  Std.Error Test.Statistic    pvalue
E           -0.01197670 0.07255330     -0.1650745 0.8690696
PRS_add     -0.04674501 0.07235187     -0.6460788 0.5190466
PRS_gxe     -0.02716207 0.07516492     -0.3613663 0.7182471
PRS_gxe x E -0.01731705 0.07376728     -0.2347525 0.8146662
```
```y$summary``` contains the target regular model summary output, when the outcome is quantitative. 

* y$summary[,"Coefficient"] : regression coefficients of all model components
* y$summary[,"Std.Error"] : standard errors of all regression coefficients
* y$summary[,"Test.Statistic"] : test statistic value of all regression coefficients
* y$summary[,"pvalue"] : p-value of all regression coefficients
  
```
  FID      IID      Risk.Values
1 "ID_801" "ID_801" "-0.250074266216069"
2 "ID_802" "ID_802" "-0.0219033515829476"
3 "ID_803" "ID_803" "0.188098265142081"
4 "ID_804" "ID_804" "0.214732615239066"
```
```y$risk.values``` contains all the calculated individual risk scores using the target dataset (e.g. Model 4), when the outcome is quantitative. The columns denote the following in order.
* y$risk.values[,"FID"] : family IDs of the target dataset
* y$risk.values[,"IID"] : individual IDs of the target dataset
* y$risk.values[,"Risk.Values] : estimated risk values of the target dataset

Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models (using ```summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = q, gxe_score = r, Model = 4)``` and ```summary_permuted_quantitative("Qpt.txt", "Qct.txt", iterations = 1000, add_score = q, gxe_score = r)```), if you choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 4) when generating risk scores. If the 'PRS_gxe x E' term is significant in Model 4, and insignificant in Model 4* (permuted p value), consider that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model). See Figure 3 in Jayasinghe et. al. for simulation verification.

See the topic summary_regular_quantitative (page 21) in manual.pdf for examples.

**Command**
```
a <- GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt")
add <- a[c("ID", "A1", "ADD_BETA")]
gxe <- a[c("ID", "A1", "INTERACTION_BETA")]
p <- PRS_quantitative(plink_path, "mydata", summary_input = add)
q <- PRS_quantitative(plink_path, "mydata", summary_input = gxe)
x <- summary_permuted_quantitative("Qpt.txt", "Qct.txt", iterations = 1000, add_score = p, gxe_score = q)
```

To perform the same on embedded data:
```
a <- GWEIS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery)
add <- a[c("ID", "A1", "ADD_BETA")]
gxe <- a[c("ID", "A1", "INTERACTION_BETA")]
p <- PRS_quantitative(plink_path, DummyData, summary_input = add)
q <- PRS_quantitative(plink_path, DummyData, summary_input = gxe)
x <- summary_permuted_quantitative(Qphe_target, Qcov_target, iterations = 1000, add_score = p, gxe_score = q)
```
```x``` outputs the p value of the fitted **permuted** model for **quantitative** outcome.


## $$\color{red}{IMPORTANT}$$
The discovery model used in ```GWAS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20)``` or ```GWAS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20)``` is as follows:
* y = b_trd.W + error
 where y is the outcome variable, b_trd is the estimated SNP effect and W is the SNP genotype.
 
The discovery model used in ```GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20)``` is as follows:
 * y = b_add.W + b_cov.E + b_cov2.E^2 + b_gxe.(WxE) + error
where y is the outcome variable, b_add is the estimated additive SNP effect, E is the covariate, W is the SNP genotype, b_cov is the estimated effect of the covariate, b_cov2 is the estimated effect of the squared covariate and b_gxe is the estimated effect of the Hadamard product of WxE.

The discovery model used in ```GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20)``` is as follows:
 * y = b_add.W + b_cov.E + b_gxe.(WxE) + error
where y is the outcome variable, b_add is the estimated additive SNP effect, E is the covariate, W is the SNP genotype, b_cov is the estimated effect of the covariate and b_gxe is the estimated effect of the Hadamard product of WxE.


The fitted (target) models in ```summary_regular_binary("Bpt.txt", "Bct.txt", add_score = NULL, gxe_score = NULL, Model)``` or ```summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = NULL, gxe_score = NULL, Model)``` are as follows:

* Model 0: y = PRS_trd + confounders + error
* Model 1: y = PRS_trd + E + PRS_trd x E + confounders + error
* Model 2: y = PRS_add + E + PRS_add x E + confounders + error
* Model 3: y = PRS_add + E + PRS_gxe x E + confounders + error
* Model 4: y = PRS_add + E + PRS_gxe + PRS_gxe x E + confounders + error
* Model 4*: permuted Model 4
* Model 5: y = PRS_add + E + E^2 + PRS_gxe + PRS_gxe x E + confounders + error
* Model 5*: permuted Model 5

where y is the outcome variable, E is the covariate of interest, PRS_trd and PRS_add are the polygenic risk scores computed using additive SNP effects of GWAS and GWEIS summary statistics respectively, and PRS_gxe is the polygenic risk scores computed using GxE interaction SNP effects of GWEIS summary statistics.


To address the potential issue of spurious GxE signals that may arise from unknown relationships between the quantitative outcome variable (y) and covariate (E), we implemented a permutation procedure on the term  PRS_gxe in the interaction component in Model 4 and denoted it as Model 4*. The purpose of this permutation was to maintain
the correlation structure between the outcome variable and other model components while specifically permuting the interaction component. Similarly, in Model 5 for binary outcome, we performed the same permutation on the PRS_gxe term in the interaction component and denoted it as Model 5*. It’s important to note that the number of permutations
performed in Models 4* or 5* was determined based on the p value obtained in Models 4 or 5 to ensure an adequate number of permutations (with a minimum of 1000). By selectively permuting the interaction term while preserving the correlation structure, we aimed to assess the robustness of the GxE effects in our models and mitigate the risk of spurious signals. This approach allowed us to obtain more reliable results and evaluate the significance of the observed GxE interactions.



# Contact 
Please contact Dovini Jayasinghe (dovini.jayasinghe@mymail.unisa.edu.au) or Md Moksedul Momin (md_moksedul.momin@mymail.unisa.edu.au) or Hong Lee (hong.lee@unisa.edu.au) for queries.

# Reference
Dovini Jayasinghe, Md. Moksedul Momin, Kerri Beckmann, Elina Hypponen, Beben Benyamin, S. Hong Lee _bioRxiv_ 2023.07.20.549816; doi: https://doi.org/10.1101/2023.07.20.549816
