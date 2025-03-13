##########################################################################################
##########################################################################################

####################################################################################################################################################################
# The following R code will allow you to complete all the EWAS requested in thechild blood DNA methylation analysis plan.
# The code also produces files summarising the variables included in the EWAS.
# You shouldn't have to rewrite or add to the following code, unless otherwise stated.
# There are just two inputs required for this analysis:
# 1) pheno: a dataframe containing all the "phenotype" data needed for this project. 
#    Each row is a sample (individual) and each column is a different variable: 
#    - Main covariates: "SampleID", "sex","child_age","mat_edu3","mat_smok2" (2-level maternal pregnancy smoking status) or "mat_smok3" (3-level maternal pregnancy smoking status),"NK","Gran","Bcell","CD8T","CD4T","Mono"
#    - Optional covariates: "anc" (ancestry), "batch" (technical covariates), "sel" (selection factor)
#    If these columns are named differently in your dataset, please rename the columns accordingly
#    Details on how to code these variables are provided in the analysis plan.
# 2) beta_matrix: a matrix of methylation illumina beta values. Each column is a sample and each row is a probe on the array (450k or EPIC). 
#    Column names must correspond to the SampleID in pheno.
####################################################################################################################################################################

##########################################################################################
### Go to the directory where to save results (change to your own working directory)
setwd("/home/user/data/WS_LIX/LIX_analyses/PdP_MC/results2/") 

##########################################################################################
### Load required packages 
#(if these are not already installed, you will have to install them as the first step) --> install.packages("openxlsx") 
library(openxlsx) # to open excel
library(readxl) # to get data out of excel and into R
library(data.table) # to process results
library(sandwich) # to estimate the standard error
library(lmtest) # to use coeftest
library(parallel) # to use multicore approach - part of base R
library(R.utils) # utility functions
library(Hmisc) # to describe function
library(psych) # to run correlation analysis
library(minfi) # to get DNA methylation data
library(xlsx) # to open excel
library(qqman) # to get QQ plots


##########################################################################################
### Set initial parameters
cohort <- "lix.EUR" #define cohort name, if several ethnic groups are tested, indicate this in the name (ie. "lix.EUR" for European ancestry children in lix)

sample.name <- "SampleID" 
cell.names <- c("NK","Gran","Bcell","CD8T","CD4T","Mono")
traits.and.covariates <- c("postSHS","postSHS_cont","postSHS_mat","postSHS_part","sex","child_age","mat_edu3","mat_smok3","NK","Gran","Bcell","CD8T","CD4T","Mono") # "batch" and "sel" are optional , # change to "mat_smok2" (2-level maternal pregnancy smoking status) or "mat_smok3" (3-level maternal pregnancy smoking status) according to your available data


##########################################################################################
### Load and check phenotype data 
# pheno: a dataframe containing all the "phenotype" data needed for this project. 
# Each row is a sample(individual) and each column is a different variable.
# Ensure all traits and covariates have been derived as specified in the analysis plan and that exclusions have been made.
# Remember to stratify by ancestry (work with only one ancestry at a time). Minimal sample size is 15 exposed children in each ancestry strata.

# Read pheno file
pheno <-read.csv("/home/isglobal.lan/user/data/WS_lix/lix_analyses/PACE_postSHS_PdP_MC/data/data/pheno_EUR_postsmoking_16072020.csv", stringsAsFactors = FALSE) ##change to your file name

# Make the ID the rownames, if needed
row.names(pheno) = pheno$SampleID

# Make sure variables are defined properly
pheno$SampleID <- as.character(pheno$SampleID)
pheno$postSHS <- as.factor(pheno$postSHS)
table(pheno$postSHS)
pheno$postSHS_mat <- as.factor(pheno$postSHS_mat)
table(pheno$postSHS_mat)
pheno$postSHS_part <- as.factor(pheno$postSHS_part)
table(pheno$postSHS_part)
pheno$postSHS_cont <- as.numeric(pheno$postSHS_cont)  ##check that it has the proper numbers, because if there was a conversion of the variable from factor to numeric in some cases it could change 0, 1, 2 to 1, 2, 3, and we want 0, 1, 2
head(pheno$postSHS_cont)
pheno$sex <- as.factor(pheno$sex)
table(pheno$sex)
pheno$child_age <- as.numeric(pheno$child_age)
names(pheno)[names(pheno) == "mat_edu"] <- "mat_edu3"
pheno$mat_edu3 <- as.factor(pheno$mat_edu3)
table(pheno$mat_edu3)
names(pheno)[names(pheno) == "pregSmok2"] <- "mat_smok3"
pheno$mat_smok3 <- as.factor (pheno$mat_smok3) ##change to "mat_smok2" (2-level maternal pregnancy smoking status) or "mat_smok3" (3-level maternal pregnancy smoking status) according to your available data
table(pheno$mat_smok3)
pheno$NK <- as.numeric(pheno$NK)
pheno$Bcell <- as.numeric(pheno$Bcell)
pheno$CD4T <- as.numeric(pheno$CD4T)
pheno$Gran <- as.numeric(pheno$Gran)
pheno$CD8T <- as.numeric(pheno$CD8T)
pheno$Mono <- as.numeric(pheno$Mono)


# Optional covariates: "anc", "batch", "sel"

# Round decimals
pheno$child_age <- round(pheno$child_age, digits=1)
pheno$NK <- round(pheno$NK, digits=4)
pheno$Bcell <- round(pheno$Bcell, digits=4)
pheno$CD4T <- round(pheno$CD4T, digits=4)
pheno$Gran <- round(pheno$Gran, digits=4)
pheno$CD8T <- round(pheno$CD8T, digits=4)
pheno$Mono <- round(pheno$Mono, digits=4)

# Check if all needed variables are present in phenotype file. If not, please correct.
for(i in 1:length(c(sample.name,traits.and.covariates,cell.names))) {
  print(ifelse(c(sample.name,traits.and.covariates,cell.names)[i] %in% colnames(pheno)==FALSE,
               paste("CAUTION: the variable called",c("sample.name",traits.and.covariates,cell.names)[i],"is missing from pheno"),
               paste("variable called",c("sample.name",traits.and.covariates,cell.names)[i],"is present in pheno")))
}


##########################################################################################
### Load methylation data
# beta_matrix: a matrix of methylation illumina beta values. Each column is a sample and each row is a probe on the array (450k or EPIC). 
# Column names must correspond to the Sample_ID in pheno.

# Load methylation data
load("/home/isglobal.lan/user/data/WS_lix/lix_preproc/methylation/Final_data/methylome_subcohort_ComBatSlide_6cells_notfitr_v4.Rdata")    #change to your methylation data

# extract beta matrix from the GenomicRatioSet, in case you data has this format.
beta_matrix <- getBeta(methylome_subcohort_ComBatSlide_notfitr)
dim(beta_matrix)


##########################################################################################
### Select subset of children for the analysis and match IDs

# If not done before stratify by ancestry origin 

# Select complete cases (within ancestry strata) to do the analysis
# If you have >5% missing for any of the variables, check with us before proceeding
pheno <- pheno[complete.cases(pheno),]
dim(pheno)
row.names(pheno) = pheno$SampleID

# Filter beta_matrix to have the same individuals as in the pheno file
beta_matrix <- beta_matrix[, pheno$SampleID]
dim(beta_matrix)

# Filter pheno file to have only children with methylation data
pheno = pheno[pheno$SampleID %in% intersect(pheno$SampleID, colnames(beta_matrix)),]
dim(pheno)

# Match IDs in beta_matrix with pheno file
table(ifelse(rownames(pheno)==colnames(beta_matrix), "MATCHED", "Not Matched"))
head(rownames(pheno))
head(colnames(beta_matrix))

# Save new pheno file
write.csv(pheno,"my-new-results/PACE_lix.EUR_pheno.SHS_20210215.csv", row.names = FALSE) ##change file name 


##########################################################################################
### Descriptive of pheno data


n <- nrow(pheno)

meanpostSHS_cont <- mean(pheno$postSHS_cont) 
sdpostSHS_cont <-  sd(pheno$postSHS_cont)
minpostSHS_cont <- min(pheno$postSHS_cont)
maxpostSHS_cont <- max(pheno$postSHS_cont)

meanchild_age <- mean(pheno$child_age) 
sdchild_age <-  sd(pheno$child_age)
minchild_age <- min(pheno$child_age)
maxchild_age <- max(pheno$child_age)

meanNK <- mean(pheno$NK)
sdNK <- sd(pheno$NK)
minNK <- min(pheno$NK)
maxNK <- max(pheno$NK)

meanBcell <- mean(pheno$Bcell)
sdBcell <- sd(pheno$Bcell)
minBcell <- min(pheno$Bcell)
maxBcell <- max(pheno$Bcell)

meanCD4T <- mean(pheno$CD4T)
sdCD4T <- sd(pheno$CD4T)
minCD4T <- min(pheno$CD4T)
maxCD4T <- max(pheno$CD4T)

meanGran <- mean(pheno$Gran)
sdGran <- sd(pheno$Gran)
minGran <- min(pheno$Gran)
maxGran <- max(pheno$Gran)

meanCD8T <- mean(pheno$CD8T)
sdCD8T <- sd(pheno$CD8T)
minCD8T <- min(pheno$CD8T)
maxCD8T <- max(pheno$CD8T)

meanMono <- mean(pheno$Mono)
sdMono <- sd(pheno$Mono)
minMono <- min(pheno$Mono)
maxMono <- max(pheno$Mono)

n_maternal_pregnancy_smoking_no <- sum(pheno$mat_smok3==0,na.rm=TRUE) ##change to "mat_smok2" (2-level maternal pregnancy smoking status) or "mat_smok3" (3-level maternal pregnancy smoking status) according to your available data
n_maternal_pregnancy_smoking_1st_trim <- sum(pheno$mat_smok3==1,na.rm=TRUE) ## if using a 2-level maternal pregnancy status then change this to: n_maternal_pregnancy_smoking_any <- sum(pheno$mat_smok2==1,na.rm=TRUE)
n_maternal_pregnancy_smoking_after_1st_trim <- sum(pheno$mat_smok3==2,na.rm=TRUE) ## if using a 2-level maternal pregnancy status then remove this line

n_anypostnatal_smoking_no <- sum(pheno$postSHS==0,na.rm=TRUE)
n_anypostnatal_smoking_yes <-sum(pheno$postSHS==1,na.rm=TRUE)

# (to know the numbers per each cathegory of variable "postSHS_cont", it must be converted to a factor as the new variable "postSHS3)
pheno$postSHS3<-as.factor(pheno$postSHS_cont)
table(pheno$postSHS3) ##check that it has the proper numbers, because if there was a conversion of the variable from numeric to factor in some cases it could change 0, 1, 2 to 1, 2, 3, and we want 0, 1, 2

n_household_smoking_no <- sum(pheno$postSHS3==0,na.rm=TRUE)
n_household_smoking_one <-sum(pheno$postSHS3==1,na.rm=TRUE)
n_household_smoking_two_or_more <-sum(pheno$postSHS3==2,na.rm=TRUE)

n_maternal_postnatal_smoking_no <- sum(pheno$postSHS_mat==0,na.rm=TRUE)
n_maternal_postnatal_smoking_yes <- sum(pheno$postSHS_mat==1,na.rm=TRUE)

n_partner_postnatal_smoking_no <- sum(pheno$postSHS_part==0,na.rm=TRUE)
n_partner_postnatal_smoking_yes <- sum(pheno$postSHS_part==1,na.rm=TRUE)

n_males <-sum(pheno$sex=="male",na.rm=TRUE)
n_females <-sum(pheno$sex=="female",na.rm=TRUE)

n_mat_education_primary <-sum(pheno$mat_edu3==1,na.rm=TRUE)
n_mat_education_secondary <-sum(pheno$mat_edu3==2,na.rm=TRUE)
n_mat_education_university <-sum(pheno$mat_edu3==3,na.rm=TRUE)

Tabledescriptives <- rbind(n, meanpostSHS_cont, sdpostSHS_cont, minpostSHS_cont,  maxpostSHS_cont, meanchild_age, sdchild_age, minchild_age, maxchild_age, meanNK, sdNK, minNK, maxNK, meanBcell, sdBcell, minBcell, maxBcell, meanCD4T, sdCD4T, minCD4T, maxCD4T, meanGran, sdGran, minGran, maxGran, meanCD8T, sdCD8T, minCD8T, maxCD8T, meanMono, sdMono, minMono, maxMono, n_males, n_females, n_maternal_pregnancy_smoking_no, n_maternal_pregnancy_smoking_1st_trim, n_maternal_pregnancy_smoking_after_1st_trim, n_anypostnatal_smoking_no, n_anypostnatal_smoking_yes, n_household_smoking_no, n_household_smoking_one, n_household_smoking_two_or_more, n_maternal_postnatal_smoking_no, n_maternal_postnatal_smoking_yes, n_partner_postnatal_smoking_no, n_partner_postnatal_smoking_yes, n_mat_education_primary, n_mat_education_secondary, n_mat_education_university) # if using mat_smok2, remove n_maternal_pregnancy_smoking_1st_trim and n_maternal_pregnancy_smoking_after_1st_trim, and add n_maternal_pregnancy_smoking_any
 
Tabledescriptives <- round(Tabledescriptives,2)
Tabledescriptives 

# Save results

write.xlsx(Tabledescriptives, file=paste0("PACE_postSHS_", cohort,"_Descriptives.postSHS_",format(Sys.Date(), 
                                                                                         "%d%m%Y"), ".xlsx"), row.names = TRUE,  append = FALSE,  showNA = TRUE,  password = NULL)


##########################################################################################
### Correlation among exposure variables within the same period and across periods (pregnancy and postnatal) - to test for potential multicollinearity

# If you have mat_smok3 then proceed:

# Tetrachoric correlation between dichotomic variables
# Polychoric correlation between smoking dichotomic variables and 3-level variables "mat_smok3" and "postSHS_cont"


tetra_1<- table(pheno$postSHS, pheno$postSHS_mat)
corr_any_mat<-tetrachoric(tetra_1, y=NULL,correct=0.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,
                       delete=TRUE)

tetra_2<- table(pheno$postSHS, pheno$postSHS_part)
corr_any_part<-tetrachoric(tetra_2, y=NULL,correct=0.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,
                          delete=TRUE)

tetra_3<- table(pheno$postSHS_mat, pheno$postSHS_part)
corr_mat_part<-tetrachoric(tetra_3, y=NULL,correct=0.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,
                           delete=TRUE)

polyc_1<- table(pheno$postSHS, pheno$mat_smok3) 
corr_any_pregsmk3<-polychoric(polyc_1,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE,  
           weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

polyc_2<- table(pheno$postSHS_mat, pheno$mat_smok3)
corr_mat_pregsmk3<-polychoric(polyc_2,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
                          weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

polyc_3<- table(pheno$postSHS_part, pheno$mat_smok3)
corr_part_pregsmk3<-polychoric(polyc_3,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
                             weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

polyc_4<- table(pheno$postSHS_cont, pheno$mat_smok3) 
corr_cont_smk_pregsmk3<-polychoric(polyc_4,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
                              weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

polyc_5<- table(pheno$postSHS_cont, pheno$postSHS_mat)
corr_cont_smk_mat<-polychoric(polyc_5,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
                              weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

polyc_6<- table(pheno$postSHS_cont, pheno$postSHS_part)
corr_cont_smk_part<-polychoric(polyc_6,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
                               weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

polyc_7<- table(pheno$postSHS_cont, pheno$postSHS)
corr_cont_smk_any<-polychoric(polyc_7,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
                              weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

  
# Save results

exposures_correlations<- rbind(corr_any_mat$rho, corr_any_part$rho, corr_mat_part$rho, corr_any_pregsmk3$rho, corr_mat_pregsmk3$rho, corr_part_pregsmk3$rho, corr_cont_smk_pregsmk3$rho, corr_cont_smk_mat$rho, corr_cont_smk_part$rho, corr_cont_smk_any$rho)
rownames(exposures_correlations)<-c("corr_anypost_matpost", "corr_anypost_partpost", "corr_matpost_partpost", "corr_anypost_pregsmk3", "corr_matpost_pregsmk3", "corr_partpost_pregsmk3", "corr_cont_postsmk_pregsmk3", "corr_cont_postsmk_matpost", "corr_cont_postsmk_partpost", "corr_cont_postsmk_anypost")
write.csv(exposures_correlations,"PACE_postSHS_lix.EUR_Correlations.postSHS_20082020.csv", row.names = TRUE) #change file name

#ALTERNATIVELY: If you have "mat_smok2", then proceed: 

  # Tetrachoric correlation between dichotomic variables
  # Polychoric correlation between smoking dichotomic variables and variables "mat_smok2" and "postSHS_cont"
  
tetra_1<- table(pheno$postSHS, pheno$postSHS_mat)
corr_any_mat<-tetrachoric(tetra_1, y=NULL,correct=0.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,
                          delete=TRUE)

tetra_2<- table(pheno$postSHS, pheno$postSHS_part)
corr_any_part<-tetrachoric(tetra_2, y=NULL,correct=0.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,
                           delete=TRUE)

tetra_3<- table(pheno$postSHS_mat, pheno$postSHS_part)
corr_mat_part<-tetrachoric(tetra_3, y=NULL,correct=0.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,
                           delete=TRUE)

tetra_4<- table(pheno$postSHS, pheno$mat_smok2)
corr_any_pregsmk2<-tetrachoric(tetra_4, y=NULL,correct=0.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,
                               delete=TRUE)

tetra_5<- table(pheno$postSHS_mat, pheno$mat_smok2)
corr_mat_pregsmk2<-tetrachoric(tetra_5, y=NULL,correct=0.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,
                               delete=TRUE)

tetra_6<- table(pheno$postSHS_part, pheno$mat_smok2)
corr_part_pregsmk2<-tetrachoric(tetra_6, y=NULL,correct=0.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,
                                delete=TRUE)

polyc_1<- table(pheno$postSHS_cont, pheno$mat_smok2) ## change to "mat_smok2" according to your data
corr_cont_smk_pregsmk2<-polychoric(polyc_1,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
                                   weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

polyc_2<- table(pheno$postSHS_cont, pheno$postSHS_mat)
corr_cont_smk_mat<-polychoric(polyc_2,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
                              weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

polyc_3<- table(pheno$postSHS_cont, pheno$postSHS_part)
corr_cont_smk_part<-polychoric(polyc_3,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
                               weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

polyc_4<- table(pheno$postSHS_cont, pheno$postSHS)
corr_cont_smk_any<-polychoric(polyc_4,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
                              weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)

# Save results

exposures_correlations<- rbind(corr_any_mat$rho, corr_any_part$rho, corr_mat_part$rho, corr_any_pregsmk2$rho, corr_mat_pregsmk2$rho, corr_part_pregsmk2$rho, corr_cont_smk_pregsmk2$rho, corr_cont_smk_mat$rho, corr_cont_smk_part$rho, corr_cont_smk_any$rho)
rownames(exposures_correlations)<-c("corr_anypost_matpost", "corr_anypost_partpost", "corr_matpost_partpost", "corr_anypost_pregsmk2", "corr_matpost_pregsmk2", "corr_partpost_pregsmk2", "corr_cont_postsmk_pregsmk2", "corr_cont_postsmk_matpost", "corr_cont_postsmk_partpost", "corr_cont_postsmk_anypost")
write.csv(exposures_correlations,"PACE_postSHS_lix.EUR_Correlations.postSHS_20082020.csv", row.names = TRUE) #change file name 


##########################################################################################
### Association between blood cell type proportion and postnatal SHS

# Linear regression models

NK_postSHS <- lm (pheno$NK ~ pheno$postSHS)
NK_postSHS<- summary(NK_postSHS)$coefficients[2,]

Bcell_postSHS <- lm (pheno$Bcell ~ pheno$postSHS)
Bcell_postSHS<- summary(Bcell_postSHS)$coefficients[2,]

CD4T_postSHS <- lm (pheno$CD4T ~ pheno$postSHS)
CD4T_postSHS<- summary(CD4T_postSHS)$coefficients[2,]

Gran_postSHS <- lm (pheno$Gran ~ pheno$postSHS)
Gran_postSHS<- summary(Gran_postSHS)$coefficients[2,]

Mono_postSHS <- lm (pheno$Mono ~ pheno$postSHS)
Mono_postSHS<- summary(Mono_postSHS)$coefficients[2,]

CD8T_postSHS <- lm (pheno$CD8T ~ pheno$postSHS)
CD8T_postSHS<- summary(CD8T_postSHS)$coefficients[2,]

NK_postSHS_cont <- lm (pheno$NK ~ pheno$postSHS_cont)
NK_postSHS_cont<- summary(NK_postSHS_cont)$coefficients[2,]

Bcell_postSHS_cont <- lm (pheno$Bcell ~ pheno$postSHS_cont)
Bcell_postSHS_cont<- summary(Bcell_postSHS_cont)$coefficients[2,]

CD4T_postSHS_cont <- lm (pheno$CD4T ~ pheno$postSHS_cont)
CD4T_postSHS_cont<- summary(CD4T_postSHS_cont)$coefficients[2,]

Gran_postSHS_cont <- lm (pheno$Gran ~ pheno$postSHS_cont)
Gran_postSHS_cont<- summary(Gran_postSHS_cont)$coefficients[2,]

Mono_postSHS_cont <- lm (pheno$Mono ~ pheno$postSHS_cont)
Mono_postSHS_cont<- summary(Mono_postSHS_cont)$coefficients[2,]

CD8T_postSHS_cont <- lm (pheno$CD8T ~ pheno$postSHS_cont)
CD8T_postSHS_cont<- summary(CD8T_postSHS_cont)$coefficients[2,]

NK_postSHS_mat <- lm (pheno$NK ~ pheno$postSHS_mat)
NK_postSHS_mat<- summary(NK_postSHS_mat)$coefficients[2,]

Bcell_postSHS_mat <- lm (pheno$Bcell ~ pheno$postSHS_mat)
Bcell_postSHS_mat<- summary(Bcell_postSHS_mat)$coefficients[2,]

CD4T_postSHS_mat <- lm (pheno$CD4T ~ pheno$postSHS_mat)
CD4T_postSHS_mat<- summary(CD4T_postSHS_mat)$coefficients[2,]

Gran_postSHS_mat <- lm (pheno$Gran ~ pheno$postSHS_mat)
Gran_postSHS_mat<- summary(Gran_postSHS_mat)$coefficients[2,]

Mono_postSHS_mat <- lm (pheno$Mono ~ pheno$postSHS_mat)
Mono_postSHS_mat<- summary(Mono_postSHS_mat)$coefficients[2,]

CD8T_postSHS_mat <- lm (pheno$CD8T ~ pheno$postSHS_mat)
CD8T_postSHS_mat<- summary(CD8T_postSHS_mat)$coefficients[2,]

NK_postSHS_part <- lm (pheno$NK ~ pheno$postSHS_part)
NK_postSHS_part<- summary(NK_postSHS_part)$coefficients[2,]

Bcell_postSHS_part <- lm (pheno$Bcell ~ pheno$postSHS_part)
Bcell_postSHS_part<- summary(Bcell_postSHS_part)$coefficients[2,]

CD4T_postSHS_part <- lm (pheno$CD4T ~ pheno$postSHS_part)
CD4T_postSHS_part<- summary(CD4T_postSHS_part)$coefficients[2,]

Gran_postSHS_part <- lm (pheno$Gran ~ pheno$postSHS_part)
Gran_postSHS_part<- summary(Gran_postSHS_part)$coefficients[2,]

Mono_postSHS_part <- lm (pheno$Mono ~ pheno$postSHS_part)
Mono_postSHS_part<- summary(Mono_postSHS_part)$coefficients[2,]

CD8T_postSHS_part <- lm (pheno$CD8T ~ pheno$postSHS_part)
CD8T_postSHS_part<- summary(CD8T_postSHS_part)$coefficients[2,]

# Save results

assoc_celltypes_smk<-rbind(NK_postSHS, Bcell_postSHS, CD4T_postSHS, Gran_postSHS, Mono_postSHS, CD8T_postSHS, NK_postSHS_cont, Bcell_postSHS_cont, CD4T_postSHS_cont, Gran_postSHS_cont, Mono_postSHS_cont, CD8T_postSHS_cont, NK_postSHS_mat, Bcell_postSHS_mat, CD4T_postSHS_mat, Gran_postSHS_mat, Mono_postSHS_mat, CD8T_postSHS_mat, NK_postSHS_part, Bcell_postSHS_part, CD4T_postSHS_part, Gran_postSHS_part, Mono_postSHS_part, CD8T_postSHS_part)
row.names(assoc_celltypes_smk)<-c("NK_postSHS", "Bcell_postSHS", "CD4T_postSHS", "Gran_postSHS", "Mono_postSHS", "CD8T_postSHS", "NK_postSHS_cont", "Bcell_postSHS_cont", "CD4T_postSHS_cont", "Gran_postSHS_cont", "Mono_postSHS_cont", "CD8T_postSHS_cont", "NK_postSHS_mat", "Bcell_postSHS_mat", "CD4T_postSHS_mat", "Gran_postSHS_mat", "Mono_postSHS_mat", "CD8T_postSHS_mat", "NK_postSHS_part", "Bcell_postSHS_part", "CD4T_postSHS_part", "Gran_postSHS_part", "Mono_postSHS_part", "CD8T_postSHS_part")
write.csv(assoc_celltypes_smk,"PACE_postSHS_lix.EUR_CellProp.postSHS_20082020.csv", row.names = TRUE) #change file name 

##########################################################################################
### Remove outliers based on IQR*3 method (in case not done before in your methylation data)

# Load function

removeOutliers<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR
  maskU <- probes > row2575[,2] + 3 * rowIQR
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}

# Remove outliers from beta_matrix and save log file

system.time(OutlierResults<-removeOutliers(beta_matrix)) 
beta_matrix<-OutlierResults[[1]]
Log <- OutlierResults[[2]]
save(Log, file = "my_LOG_FILES/PACE_postSHS_lix.EUR_Outlier_log.postSHS_20082020.Rdata") ## change file name
Log1 <- Log[order(Log$N_for_probe),]
write.table(Log1, file=paste0("RESULTS/results4collaborators/feb2019/PACE_postSHS_",cohort,"TrimmingLog.postSHS_",format(Sys.Date(), 
        "%d%m%Y"), ".txt"),sep = "\t", row.names = TRUE,  col.names= TRUE) ## change file name 

##########################################################################################
### EWAS Main models (all cohorts)

### EWAS Model 1
### Methylation = postSHS + mat_smok3 (or mat_smok2) + cov + optional cov

### Initial checks
# Check that rows are samples and columns are probes on the array (450k or EPIC)
dim(beta_matrix)
# If not, then transpose betas so that rows are samples and columns are probes on the array (450k or EPIC)
beta_matrix <- t(beta_matrix)


# Check that all IDs are the same - need to be in the same order!!
pheno <- pheno[match(rownames(beta_matrix), pheno$SampleID),]
all(rownames(beta_matrix)==rownames(pheno))
#TRUE
head(row.names(beta_matrix))
head(row.names(pheno))

### If you have maternal smoking  during pregnancy coded in 3 levels (mat_smok3)

# Load functions for running the model adjusted for mat_smok3. 
# Adapt it to the number of covariables needed in your study (add optional covariables).

LModel1 <- function (data)
{
  lm(CpGi~postSHS + mat_smok3 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}



apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel1(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel1')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_beta = estimate[,"postSHS1"],
                     postSHS_SE = SE[,"postSHS1"],
                     postSHS_pval = pval[,"postSHS1"],
                     postSHS_N = N[,1],
                     mat_smok3_1_beta = estimate[,"mat_smok31"],
                     mat_smok3_1_SE = SE[,"mat_smok31"],
                     mat_smok3_1_pval = pval[,"mat_smok31"],
                     mat_smok3_1_N = N[,1],
                     mat_smok3_2_beta = estimate[,"mat_smok32"],
                     mat_smok3_2_SE = SE[,"mat_smok32"],
                     mat_smok3_2_pval = pval[,"mat_smok32"],
                     mat_smok3_2_N = N[,1])

model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
rownames(model.data)<-c()


#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok3_1_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda3 <- qchisq(median(model.data$mat_smok3_2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)

# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model1_13112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model1_13112020.txt") #change file name

lambda <- rbind(lambda1, lambda2, lambda3)
write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model1_Lambda_",format(Sys.Date(), 
                                                                              "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot
jpeg("PACE_postSHS_lix.EUR_Model1_QQPlot1_13112020.jpg")
qq(model.data$postSHS_pval,main="PACE_postSHS_lix.EUR_Model1_QQPlot1_13112020") # please cha> dev.off()me
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model1_QQPlot2_13112020.jpg")
qq(model.data$mat_smok3_1_pval,main="PACE_postSHS_lix.EUR_Model1_QQPlot2_13112020") # please> dev.off()e name
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model1_QQPlot3_13112020.jpg")
qq(model.data$mat_smok3_2_pval,main="PACE_postSHS_lix.EUR_Model1_QQPlot3_13112020") # please> dev.off()e name
dev.off()

### ALTERNATIVE: If you have maternal smoking   during pregnancy coded in 2 levels (mat_smok2)

# Load function

LModel1 <- function (data)
{
  lm(CpGi~postSHS + mat_smok2 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}

apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel1(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel1')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_beta = estimate[,"postSHS1"],
                     postSHS_SE = SE[,"postSHS1"],
                     postSHS_pval = pval[,"postSHS1"],
                     postSHS_N = N[,1],
                     mat_smok2_beta = estimate[,"mat_smok21"],
                     mat_smok2_SE = SE[,"mat_smok21"],
                     mat_smok2_pval = pval[,"mat_smok21"],
                     mat_smok2_N = N[,1])

model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(9, 1, 2, 3, 4, 5, 6, 7, 8)]
rownames(model.data)<-c()


#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)


# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model1_14112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model1_14112020.txt") #change file name

lambda <- rbind(lambda1, lambda2)
write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model1_Lambda_matsmok2",format(Sys.Date(), 
                                                                                      "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot
jpeg("PACE_postSHS_lix.EUR_Model1_QQPlot1_14112020.jpg")
qq(model.data$postSHS_pval,main="PACE_postSHS_lix.EUR_Model1_QQPlot1_14112020") # please cha> dev.off()me
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model1_QQPlot2_14112020.jpg")
qq(model.data$mat_smok2_pval,main="PACE_postSHS_lix.EUR_Model1_QQPlot2_14112020") # please> dev.off()e name
dev.off()

##########################################################################################

### EWAS Model 2

### = postSHS_cont + mat_smok3 (or mat_smok2) + cov + optional cov


### Initial checks

# Check that rows are samples and columns are probes on the array (450k or EPIC)
dim(beta_matrix)
# If not, then transpose betas so that rows are samples and columns are probes on the array (450k or EPIC)
beta_matrix <- t(beta_matrix)


# Check that all IDs are the same - need to be in the same order!!
pheno <- pheno[match(rownames(beta_matrix), pheno$SampleID),]
all(rownames(beta_matrix)==rownames(pheno))
#TRUE
head(row.names(beta_matrix))
head(row.names(pheno))

### If you have maternal smoking  during pregnancy coded in 3 levels (mat_smok3)

# Load functions for running the model adjusted for mat_smok3. 
# Adapt it to the number of covariables needed in your study (add optional covariables).

LModel2 <- function (data)
{
  lm(CpGi~postSHS_cont + mat_smok3 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}



apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel2(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel2')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_cont_beta = estimate[,"postSHS_cont"],
                     postSHS_cont_SE = SE[,"postSHS_cont"],
                     postSHS_cont_pval = pval[,"postSHS_cont"],
                     postSHS_cont_N = N[,1],
                     mat_smok3_1_beta = estimate[,"mat_smok31"],
                     mat_smok3_1_SE = SE[,"mat_smok31"],
                     mat_smok3_1_pval = pval[,"mat_smok31"],
                     mat_smok3_1_N = N[,1],
                     mat_smok3_2_beta = estimate[,"mat_smok32"],
                     mat_smok3_2_SE = SE[,"mat_smok32"],
                     mat_smok3_2_pval = pval[,"mat_smok32"], 
                     mat_smok3_2_N = N[,1])


model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
rownames(model.data)<-c()


#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_cont_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok3_1_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda3 <- qchisq(median(model.data$mat_smok3_2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)

# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model2_13112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model2_13112020.txt") #change file name

lambda <- rbind(lambda1, lambda2, lambda3)
write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model2_Lambda_",format(Sys.Date(), 
                                                                              "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot
jpeg("PACE_postSHS_lix.EUR_Model2_QQPlot1_13112020.jpg")
qq(model.data$postSHS_cont_pval,main="PACE_postSHS_lix.EUR_Model2_QQPlot1_13112020") # please cha> dev.off()me
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model2_QQPlot2_13112020.jpg")
qq(model.data$mat_smok3_1_pval,main="PACE_postSHS_lix.EUR_Model2_QQPlot2_13112020") # please> dev.off()e name
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model2_QQPlot3_13112020.jpg")
qq(model.data$mat_smok3_2_pval,main="PACE_postSHS_lix.EUR_Model2_QQPlot3_13112020") # please> dev.off()e name
dev.off()

### ALTERNATIVE: If you have maternal smoking   during pregnancy coded in 2 levels (mat_smok2)

# Load function

LModel2 <- function (data)
{
  lm(CpGi~postSHS_cont + mat_smok2 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}



apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel2(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel2')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_cont_beta = estimate[,"postSHS_cont"],
                     postSHS_cont_SE = SE[,"postSHS_cont"],
                     postSHS_cont_pval = pval[,"postSHS_cont"],
                     postSHS_cont_N = N[,1],
                     mat_smok2_beta = estimate[,"mat_smok21"],
                     mat_smok2_SE = SE[,"mat_smok21"],
                     mat_smok2_pval = pval[,"mat_smok21"],
                     mat_smok2_N = N[,1])


model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(9, 1, 2, 3, 4, 5, 6, 7, 8)]
rownames(model.data)<-c()

#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_cont_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)


# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model2_14112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model2_14112020.txt") #change file name

lambda <- rbind(lambda1, lambda2)
write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model2_Lambda_matsmok2",format(Sys.Date(), 
                                                                                      "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot
jpeg("PACE_postSHS_lix.EUR_Model2_QQPlot1_14112020.jpg")
qq(model.data$postSHS_cont_pval,main="PACE_postSHS_lix.EUR_Model2_QQPlot1_14112020") # please cha> dev.off()me
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model2_QQPlot2_14112020.jpg")
qq(model.data$mat_smok2_pval,main="PACE_postSHS_lix.EUR_Model2_QQPlot2_14112020") # please> dev.off()e name
dev.off()


##########################################################################################
### EWAS Model 3
### Methylation = postSHS_mat + mat_smok3 (or mat_smok2) + cov + optional cov

### Initial checks

# Check that rows are samples and columns are probes on the array (450k or EPIC)
dim(beta_matrix)
# If not, then transpose betas so that rows are samples and columns are probes on the array (450k or EPIC)
beta_matrix <- t(beta_matrix)


# Check that all IDs are the same - need to be in the same order!!
pheno <- pheno[match(rownames(beta_matrix), pheno$SampleID),]
all(rownames(beta_matrix)==rownames(pheno))
#TRUE
head(row.names(beta_matrix))
head(row.names(pheno))

### If you have maternal smoking  during pregnancy coded in 3 levels (mat_smok3)

# Load functions for running the model adjusted for mat_smok3. 
# Adapt it to the number of covariables needed in your study (add optional covariables).

LModel3 <- function (data)
{
  lm(CpGi~postSHS_mat + mat_smok3 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}



apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel3(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel3')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_mat_beta = estimate[,"postSHS_mat1"],
                     postSHS_mat_SE = SE[,"postSHS_mat1"],
                     postSHS_mat_pval = pval[,"postSHS_mat1"],
                     postSHS_mat_N = N[,1],
                     mat_smok3_1_beta = estimate[,"mat_smok31"],
                     mat_smok3_1_SE = SE[,"mat_smok31"],
                     mat_smok3_1_pval = pval[,"mat_smok31"],
                     mat_smok3_1_N = N[,1],
                     mat_smok3_2_beta = estimate[,"mat_smok32"],
                     mat_smok3_2_SE = SE[,"mat_smok32"],
                     mat_smok3_2_pval = pval[,"mat_smok32"],
                     mat_smok3_2_N = N[,1])


model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
rownames(model.data)<-c()

#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_mat_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok3_1_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda3 <- qchisq(median(model.data$mat_smok3_2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)

# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model3_13112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model3_13112020.txt") #change file name

lambda <- rbind(lambda1, lambda2, lambda3)
write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model3_Lambda_",format(Sys.Date(), 
                                                                              "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot
jpeg("PACE_postSHS_lix.EUR_Model3_QQPlot1_13112020.jpg")
qq(model.data$postSHS_mat_pval,main="PACE_postSHS_lix.EUR_Model3_QQPlot1_13112020") # please cha> dev.off()me
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model3_QQPlot2_13112020.jpg")
qq(model.data$mat_smok3_1_pval,main="PACE_postSHS_lix.EUR_Model3_QQPlot2_13112020") # please> dev.off()e name
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model3_QQPlot3_13112020.jpg")
qq(model.data$mat_smok3_2_pval,main="PACE_postSHS_lix.EUR_Model3_QQPlot3_13112020") # please> dev.off()e name
dev.off()

### ALTERNATIVE: If you have maternal smoking   during pregnancy coded in 2 levels (mat_smok2)

# Load function

LModel3 <- function (data)
{
  lm(CpGi~postSHS_mat + mat_smok2 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}



apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel3(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel3')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_mat_beta = estimate[,"postSHS_mat1"],
                     postSHS_mat_SE = SE[,"postSHS_mat1"],
                     postSHS_mat_pval = pval[,"postSHS_mat1"],
                     postSHS_mat_N = N[,1],
                     mat_smok2_beta = estimate[,"mat_smok21"],
                     mat_smok2_SE = SE[,"mat_smok21"],
                     mat_smok2_pval = pval[,"mat_smok21"],
                     mat_smok2_N = N[,1])

model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(9, 1, 2, 3, 4, 5, 6, 7, 8)]
rownames(model.data)<-c()

#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_mat_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)


# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model3_14112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model3_14112020.txt") #change file name

lambda <- rbind(lambda1, lambda2)
write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model3_Lambda_matsmok2",format(Sys.Date(), 
                                                                                      "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot
jpeg("PACE_postSHS_lix.EUR_Model3_QQPlot1_14112020.jpg")
qq(model.data$postSHS_mat_pval,main="PACE_postSHS_lix.EUR_Model3_QQPlot1_14112020") # please cha> dev.off()me
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model3_QQPlot2_14112020.jpg")
qq(model.data$mat_smok2_pval,main="PACE_postSHS_lix.EUR_Model3_QQPlot2_14112020") # please> dev.off()e name
dev.off()

##########################################################################################
### EWAS Model 4
### Methylation = postSHS_part + mat_smok3 (or mat_smok2) + cov + optional cov

### Initial checks

# Check that rows are samples and columns are probes on the array (450k or EPIC)
dim(beta_matrix)
# If not, then transpose betas so that rows are samples and columns are probes on the array (450k or EPIC)
beta_matrix <- t(beta_matrix)


# Check that all IDs are the same - need to be in the same order!!
pheno <- pheno[match(rownames(beta_matrix), pheno$SampleID),]
all(rownames(beta_matrix)==rownames(pheno))
#TRUE
head(row.names(beta_matrix))
head(row.names(pheno))

### If you have maternal smoking  during pregnancy coded in 3 levels (mat_smok3)

# Load functions for running the model adjusted for mat_smok3. 
# Adapt it to the number of covariables needed in your study (add optional covariables).

LModel4 <- function (data)
{
  lm(CpGi~postSHS_part + mat_smok3 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}



apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel4(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel4')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_part_beta = estimate[,"postSHS_part1"],
                     postSHS_part_SE = SE[,"postSHS_part1"],
                     postSHS_part_pval = pval[,"postSHS_part1"],
                     postSHS_part_N = N[,1],
                     mat_smok3_1_beta = estimate[,"mat_smok31"],
                     mat_smok3_1_SE = SE[,"mat_smok31"],
                     mat_smok3_1_pval = pval[,"mat_smok31"],
                     mat_smok3_1_N = N[,1],
                     mat_smok3_2_beta = estimate[,"mat_smok32"],
                     mat_smok3_2_SE = SE[,"mat_smok32"],
                     mat_smok3_2_pval = pval[,"mat_smok32"],
                     mat_smok3_2_N = N[,1])

model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
rownames(model.data)<-c()


#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_part_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok3_1_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda3 <- qchisq(median(model.data$mat_smok3_2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)

# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model4_13112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model4_13112020.txt") #change file name

lambda <- rbind(lambda1, lambda2, lambda3)
write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model4_Lambda_",format(Sys.Date(), 
                                                                              "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot
jpeg("PACE_postSHS_lix.EUR_Model4_QQPlot1_13112020.jpg")
qq(model.data$postSHS_part_pval,main="PACE_postSHS_lix.EUR_Model4_QQPlot1_13112020") # please cha> dev.off()me
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model4_QQPlot2_13112020.jpg")
qq(model.data$mat_smok3_1_pval,main="PACE_postSHS_lix.EUR_Model4_QQPlot2_13112020") # please> dev.off()e name
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model4_QQPlot3_13112020.jpg")
qq(model.data$mat_smok3_2_pval,main="PACE_postSHS_lix.EUR_Model4_QQPlot3_13112020") # please> dev.off()e name
dev.off()

### ALTERNATIVE: If you have maternal smoking   during pregnancy coded in 2 levels (mat_smok2)

# Load function

LModel4 <- function (data)
{
  lm(CpGi~postSHS_part + mat_smok2 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}



apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel4(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel4')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_part_beta = estimate[,"postSHS_part1"],
                     postSHS_part_SE = SE[,"postSHS_part1"],
                     postSHS_part_pval = pval[,"postSHS_part1"],
                     postSHS_part_N = N[,1],
                     mat_smok2_beta = estimate[,"mat_smok21"],
                     mat_smok2_SE = SE[,"mat_smok21"],
                     mat_smok2_pval = pval[,"mat_smok21"],
                     mat_smok2_N = N[,1])

model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(9, 1, 2, 3, 4, 5, 6, 7, 8)]
rownames(model.data)<-c()

#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_part_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)


# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model4_14112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model4_14112020.txt") #change file name

lambda <- rbind(lambda1, lambda2)
write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model4_Lambda_matsmok2",format(Sys.Date(), 
                                                                                      "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot
jpeg("PACE_postSHS_lix.EUR_Model4_QQPlot1_14112020.jpg")
qq(model.data$postSHS_part_pval,main="PACE_postSHS_lix.EUR_Model4_QQPlot1_14112020") # please cha> dev.off()me
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model4_QQPlot2_14112020.jpg")
qq(model.data$mat_smok2_pval,main="PACE_postSHS_lix.EUR_Model4_QQPlot2_14112020") # please> dev.off()e name
dev.off()

##########################################################################################

### EWAS Model 5
### Methyl = xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# This model is inciduded in the code for Cotinine, whcih might have a different sample size!!!

##########################################################################################
### Interaction models (all cohorts)

### EWAS Model 6
### Methylation = postSHS*mat_smok3 (or mat_smok2) + cov + optional cov

### Initial checks

# Check that rows are samples and columns are probes on the array (450k or EPIC)
dim(beta_matrix)
# If not, then transpose betas so that rows are samples and columns are probes on the array (450k or EPIC)
beta_matrix <- t(beta_matrix)


# Check that all IDs are the same - need to be in the same order!!
pheno <- pheno[match(rownames(beta_matrix), pheno$SampleID),]
all(rownames(beta_matrix)==rownames(pheno))
#TRUE
head(row.names(beta_matrix))
head(row.names(pheno))

### If you have maternal smoking  during pregnancy coded in 3 levels (mat_smok3)

# Load functions for running the model adjusted for mat_smok3. 
# Adapt it to the number of covariables needed in your study (add optional covariables).

LModel6 <- function (data)
{
  lm(CpGi~postSHS*mat_smok3 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}



apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel6(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel6')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_beta = estimate[,"postSHS1"],
                     postSHS_SE = SE[,"postSHS1"],
                     postSHS_pval = pval[,"postSHS1"],
                     postSHS_N = N[,1],
                     mat_smok3_1_beta = estimate[,"mat_smok31"],
                     mat_smok3_1_SE = SE[,"mat_smok31"],
                     mat_smok3_1_pval = pval[,"mat_smok31"],
                     mat_smok3_1_N = N[,1],
                     mat_smok3_2_beta = estimate[,"mat_smok32"],
                     mat_smok3_2_SE = SE[,"mat_smok32"],
                     mat_smok3_2_pval = pval[,"mat_smok32"],
                     mat_smok3_2_N = N[,1],
                     int_postSHS_mat_smok3_1_beta = estimate[,"postSHS1:mat_smok31"],
                     int_postSHS_mat_smok3_1_SE = SE[,"postSHS1:mat_smok31"],
                     int_postSHS_mat_smok3_1_pval = pval[,"postSHS1:mat_smok31"],
                     int_postSHS_mat_smok3_1_N = N[,1],
                     int_postSHS_mat_smok3_2_beta = estimate[,"postSHS1:mat_smok32"],
                     int_postSHS_mat_smok3_2_SE = SE[,"postSHS1:mat_smok32"],
                     int_postSHS_mat_smok3_2_pval = pval[,"postSHS1:mat_smok32"],
                     int_postSHS_mat_smok3_2_N = N[,1])

model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(21, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)]
rownames(model.data)<-c()


#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok3_1_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda3 <- qchisq(median(model.data$mat_smok3_2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda4 <- qchisq(median(model.data$int_postSHS_mat_smok3_1_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda5 <- qchisq(median(model.data$int_postSHS_mat_smok3_2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)


# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model6_13112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model6_13112020.txt") #change file name

lambda <- rbind(lambda1, lambda2, lambda3, lambda4, lambda5)

write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model6_Lambda_",format(Sys.Date(), 
                                                                              "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot

jpeg("PACE_postSHS_lix.EUR_Model6_QQPlot1_20082020.jpg")
qq(model.data$postSHS_pval,main="PACE_postSHS_lix.EUR_Model6_QQPlot1_20082020") # please change file name 
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model6_QQPlot2_20082020.jpg")
qq(model.data$mat_smok3_1_pval,main="PACE_postSHS_lix.EUR_Model6_QQPlot2_20082020") # please change file name 
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model6_QQPlot3_20082020.jpg")
qq(model.data$mat_smok3_2_pval,main="PACE_postSHS_lix.EUR_Model6_QQPlot3_20082020") # please change file name 
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model6_QQPlot4_20082020.jpg")
qq(model.data$int_postSHS_mat_smok3_1_pval,main="PACE_postSHS_lix.EUR_Model6_QQPlot4_20082020") # please change file name 
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model6_QQPlot5_20082020.jpg")
qq(model.data$int_postSHS_mat_smok3_2_pval,main="PACE_postSHS_lix.EUR_Model6_QQPlot5_20082020") # please change file name 
dev.off()


### ALTERNATIVE: If you have maternal smoking   during pregnancy coded in 2 levels (mat_smok2)

# Load function

LModel6 <- function (data)
{
  lm(CpGi~postSHS*mat_smok2 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}



apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel6(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel6')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_beta = estimate[,"postSHS1"],
                     postSHS_SE = SE[,"postSHS1"],
                     postSHS_pval = pval[,"postSHS1"],
                     postSHS_N = N[,1],
                     mat_smok2_beta = estimate[,"mat_smok21"],
                     mat_smok2_SE = SE[,"mat_smok21"],
                     mat_smok2_pval = pval[,"mat_smok21"],
                     mat_smok2_N = N[,1],
                     int_postSHS_mat_smok2_beta = estimate[,"postSHS1:mat_smok21"],
                     int_postSHS_mat_smok2_SE = SE[,"postSHS1:mat_smok21"],
                     int_postSHS_mat_smok2_pval = pval[,"postSHS1:mat_smok21"],
                     int_postSHS_mat_smok2_N = N[,1])

model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
rownames(model.data)<-c()


#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda3 <- qchisq(median(model.data$int_postSHS_mat_smok2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)

# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model6_14112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model6_14112020.txt") #change file name

lambda <- rbind(lambda1, lambda2, lambda3)
write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model6_Lambda_matsmok2",format(Sys.Date(), 
                                                                                      "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot
jpeg("PACE_postSHS_lix.EUR_Model6_QQPlot1_20082020.jpg")
qq(model.data$postSHS_pval,main="PACE_postSHS_lix.EUR_Model6_QQPlot1_20082020") # please change file name 
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model6_QQPlot2_20082020.jpg")
qq(model.data$mat_smok2_pval,main="PACE_postSHS_lix.EUR_Model6_QQPlot2_20082020") # please change file name 
dev.off()


jpeg("PACE_postSHS_lix.EUR_Model6_QQPlot4_20082020.jpg")
qq(model.data$int_postSHS_mat_smok2_pval,main="PACE_postSHS_lix.EUR_Model6_QQPlot4_20082020") # please change file name 
dev.off()


##########################################################################################
### EWAS Model 7
### Methylation = postSHS_mat*mat_smok3 (or mat_smok2) + cov + optional cov


### Initial checks

# Check that rows are samples and columns are probes on the array (450k or EPIC)
dim(beta_matrix)
# If not, then transpose betas so that rows are samples and columns are probes on the array (450k or EPIC)
beta_matrix <- t(beta_matrix)


# Check that all IDs are the same - need to be in the same order!!
pheno <- pheno[match(rownames(beta_matrix), pheno$SampleID),]
all(rownames(beta_matrix)==rownames(pheno))
#TRUE
head(row.names(beta_matrix))
head(row.names(pheno))

### If you have maternal smoking  during pregnancy coded in 3 levels (mat_smok3)

# Load functions for running the model adjusted for mat_smok3. 
# Adapt it to the number of covariables needed in your study (add optional covariables).

LModel7 <- function (data)
{
  lm(CpGi~postSHS_mat*mat_smok3 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}



apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel7(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel7')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_mat_beta = estimate[,"postSHS_mat1"],
                     postSHS_mat_SE = SE[,"postSHS_mat1"],
                     postSHS_mat_pval = pval[,"postSHS_mat1"],
                     postSHS_mat_N = N[,1],
                     mat_smok3_1_beta = estimate[,"mat_smok31"],
                     mat_smok3_1_SE = SE[,"mat_smok31"],
                     mat_smok3_1_pval = pval[,"mat_smok31"],
                     mat_smok3_1_N = N[,1],
                     mat_smok3_2_beta = estimate[,"mat_smok32"],
                     mat_smok3_2_SE = SE[,"mat_smok32"],
                     mat_smok3_2_pval = pval[,"mat_smok32"],
                     mat_smok3_2_N = N[,1],
                     int_postSHS_mat_mat_smok3_1_beta = estimate[,"postSHS_mat1:mat_smok31"],
                     int_postSHS_mat_mat_smok3_1_SE = SE[,"postSHS_mat1:mat_smok31"],
                     int_postSHS_mat_mat_smok3_1_pval = pval[,"postSHS_mat1:mat_smok31"],
                     int_postSHS_mat_mat_smok3_1_N = N[,1],
                     int_postSHS_mat_mat_smok3_2_beta = estimate[,"postSHS_mat1:mat_smok32"],
                     int_postSHS_mat_mat_smok3_2_SE = SE[,"postSHS_mat1:mat_smok32"],
                     int_postSHS_mat_mat_smok3_2_pval = pval[,"postSHS_mat1:mat_smok32"],
                     int_postSHS_mat_mat_smok3_2_N = N[,1])


model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(21, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)]
rownames(model.data)<-c()


#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_mat_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok3_1_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda3 <- qchisq(median(model.data$mat_smok3_2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda4 <- qchisq(median(model.data$int_postSHS_mat_mat_smok3_1_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda5 <- qchisq(median(model.data$int_postSHS_mat_mat_smok3_2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)


# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model7_13112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model7_13112020.txt") #change file name

lambda <- rbind(lambda1, lambda2, lambda3, lambda4, lambda5)

write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model7_Lambda_",format(Sys.Date(), 
                                                                              "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot

jpeg("PACE_postSHS_lix.EUR_Model7_QQPlot1_20082020.jpg")
qq(model.data$postSHS_mat_pval,main="PACE_postSHS_lix.EUR_Model7_QQPlot1_20082020") # please change file name 
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model7_QQPlot2_20082020.jpg")
qq(model.data$mat_smok3_1_pval,main="PACE_postSHS_lix.EUR_Model7_QQPlot2_20082020") # please change file name 
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model7_QQPlot3_20082020.jpg")
qq(model.data$mat_smok3_2_pval,main="PACE_postSHS_lix.EUR_Model7_QQPlot3_20082020") # please change file name 
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model7_QQPlot4_20082020.jpg")
qq(model.data$int_postSHS_mat_mat_smok3_1_pval,main="PACE_postSHS_lix.EUR_Model7_QQPlot4_20082020") # please change file name 
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model7_QQPlot5_20082020.jpg")
qq(model.data$int_postSHS_mat_mat_smok3_2_pval,main="PACE_postSHS_lix.EUR_Model7_QQPlot5_20082020") # please change file name 
dev.off()


### ALTERNATIVE: If you have maternal smoking   during pregnancy coded in 2 levels (mat_smok2)

# Load function

LModel7 <- function (data)
{
  lm(CpGi~postSHS_mat*mat_smok2 + sex + mat_edu3 + child_age + NK + Bcell + CD4T + Gran + CD8T + Mono,data=data)
}

apply_model <- function( beta, phen, modelname )
{
  
  phen$CpGi <- beta
  # N[i,] <- length(which(!is.na(phen$CpGi)))
  coeff <- coef(summary(LModel7(phen)))
  modelname <- paste0(modelname,"(phen)")
  coef(summary( eval(parse(text=modelname)) ))
  
  return(list(estimate <- coeff[,1],
              SE <- coeff[,2],
              tval <- coeff[,3],
              pval <-  coeff[,4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows )
{
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if(dim(mat)[2] >1){
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}


Ncpgs <- ncol(beta_matrix)

#Run EWAS

results_list <- apply(beta_matrix, 2, apply_model, phen = pheno, modelname = 'LModel7')

#Process Results

estimate <- get_model_data( 1, results_list, Ncpgs )
SE <- get_model_data( 2, results_list, Ncpgs )
tval <- get_model_data( 3, results_list, Ncpgs )
pval <- get_model_data( 4, results_list, Ncpgs )
N <- get_model_data( 5, results_list, Ncpgs )

model.data <- cbind( postSHS_mat_beta = estimate[,"postSHS_mat1"],
                     postSHS_mat_SE = SE[,"postSHS_mat1"],
                     postSHS_mat_pval = pval[,"postSHS_mat1"],
                     postSHS_mat_N = N[,1],
                     mat_smok2_beta = estimate[,"mat_smok21"],
                     mat_smok2_SE = SE[,"mat_smok21"],
                     mat_smok2_pval = pval[,"mat_smok21"],
                     mat_smok2_N = N[,1],
                     int_postSHS_mat_mat_smok2_beta = estimate[,"postSHS_mat1:mat_smok21"],
                     int_postSHS_mat_mat_smok2_SE = SE[,"postSHS_mat1:mat_smok21"],
                     int_postSHS_mat_mat_smok2_pval = pval[,"postSHS_mat1:mat_smok21"],
                     int_postSHS_mat_mat_smok2_N = N[,1])

model.data<-data.frame(model.data)
model.data$probeID <- rownames(model.data)
model.data<-model.data[,c(13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
rownames(model.data)<-c()


#Calculate lambdas

lambda1 <- qchisq(median(model.data$postSHS_mat_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda2 <- qchisq(median(model.data$mat_smok2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda3 <- qchisq(median(model.data$int_postSHS_mat_mat_smok2_pval,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)

# Export results
write.table(model.data, "PACE_postSHS_lix.EUR_Model7_14112020.txt", na="NA") #change file name 
gzip("PACE_postSHS_lix.EUR_Model7_14112020.txt") #change file name

lambda <- rbind(lambda1, lambda2, lambda3)
write.table(lambda, file=paste0("PACE_postSHS",cohort,"Model7_Lambda_matsmok2",format(Sys.Date(), 
                                                                                      "%d%m%Y"), ".txt"), row.names = TRUE,  append = FALSE)

# QQ plot
jpeg("PACE_postSHS_lix.EUR_Model7_QQPlot1_20082020.jpg")
qq(model.data$postSHS_mat_pval,main="PACE_postSHS_lix.EUR_Model7_QQPlot1_20082020") # please change file name 
dev.off()

jpeg("PACE_postSHS_lix.EUR_Model7_QQPlot2_20082020.jpg")
qq(model.data$mat_smok2_pval,main="PACE_postSHS_lix.EUR_Model7_QQPlot2_20082020") # please change file name 
dev.off()


jpeg("PACE_postSHS_lix.EUR_Model7_QQPlot4_20082020.jpg")
qq(model.data$int_postSHS_mat_mat_smok2_pval,main="PACE_postSHS_lix.EUR_Model7_QQPlot4_20082020") # please change file name 
dev.off()


##########################################################################################


# THANK YOU!!!

