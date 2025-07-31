library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(stringr) 


# loading clinical data
df_clinical_patient = read.delim("prad_tcga_gdc/prad_tcga_gdc/data_clinical_patient.txt",skip=4,header = TRUE) %>% column_to_rownames(., var = "PATIENT_ID")
df_clinical_patient <- df_clinical_patient[
  !is.na(df_clinical_patient$OS_STATUS) & df_clinical_patient$OS_STATUS != "", 
]

# extracting status
status <- tibble(df_clinical_patient$OS_STATUS) %>% separate(`df_clinical_patient$OS_STATUS`,c("Status","Explanation"))

# loading mutation data
library(data.table)
df_mutations <- read.delim("prad_tcga_gdc/prad_tcga_gdc/data_mutations.txt",
                           header = TRUE,
                           sep = "\t",
                           quote = "",
                           fill = TRUE,
                           comment.char = "#",
                           stringsAsFactors = FALSE,
                           check.names = FALSE)

# converting to wide table
df_mutations_smpls <- df_mutations %>% 
  pivot_wider(id_cols = Hugo_Symbol,names_from = Tumor_Sample_Barcode,values_from = Variant_Type,values_fn = length) %>%
  mutate_all(funs(ifelse(is.na(.), 0, .))) %>% 
  column_to_rownames(., var = "Hugo_Symbol")

# uncomment this line if you need only 0 and 1
df_mutations_smpls[df_mutations_smpls != 0] <- 1 

# selecting TCGA samples having mutations data
# we have to cut last characters in mutation barcodes so they coincide with patient IDs
colnames(df_mutations_smpls) <- str_sub(colnames(df_mutations_smpls), end = -5)
# create a new survival object
df_patient_mut <- df_clinical_patient[colnames(df_mutations_smpls),]
status_mut <- tibble(df_patient_mut$OS_STATUS) %>% separate(`df_patient_mut$OS_STATUS`,c("Status","Explanation"))
surv_obj_mut <- Surv(time = df_patient_mut$OS_MONTHS, event = as.numeric(status_mut$Status))

# feature selection
library("glmpath")
library("glmnet")
library("penalized")

# prepring data without NAs and 0 times
df_patient_mut <- df_clinical_patient[colnames(df_mutations_smpls),]
bools <- ((!is.na(df_patient_mut$OS_MONTHS)) & (df_patient_mut$OS_MONTHS!=0)) 
patient_time <- df_patient_mut$OS_MONTHS[bools]
status_mut <- tibble(df_patient_mut$OS_STATUS) %>% separate(`df_patient_mut$OS_STATUS`,c("Status","Explanation"))
status_mut_col <- status_mut$Status[bools]
df_patient_mut_b <- df_mutations_smpls[bools]
# creating survival object without NAs and 0 times
surv_obj_mut <- Surv(time = patient_time, event = as.numeric(status_mut_col))

# running glmnet
cvfit <- cv.glmnet(data.matrix(t(df_patient_mut_b)),surv_obj_mut,family="cox",type.measure = "C")
plot(cvfit)
print(cvfit)

opt_i <- which.max(cvfit$cvm)
opt_l <- cvfit$lambda[opt_i]

fit_glm <- glmnet(as.data.frame(t(df_patient_mut_b)),surv_obj_mut,family="cox") # , alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)
#opt_l <- 0.114700
# analysing results
cfs = coef(cvfit,s=opt_l)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
write.csv(meaning_vals,file='survival_mut_glmnet.csv')
print(paste(meaning_coefs,collapse=" + "))
# if we want to run coxph we have to copy-paste the string
fit.coxph <- coxph(surv_obj_mut ~ CDHR1 + ARMH4 + SLC47A2 + AOC2 + CNN1 + SAV1 + OR3A2 +
                     C17orf78 + CDH2 + SHC2 + C2orf42 + RAC2 + EFNA5 + RPF2 + NSDHL + NRG2 + NPPA +
                     PRAMEF1 + FAM43B + OR52E4 + VPS18 + LRRC28 + CCDC102A + TERF2IP + HOXB8
                   + SOX9 + CBLN2 + ZNF225 + GOLIM4 + KLKB1 , data = as.data.frame(t(df_patient_mut_b)))
ggforest(fit.coxph, data = as.data.frame(t(df_patient_mut_b)))


# now we can plot survival curves for specific mutations
fit_mut <- survfit(surv_obj_mut ~ CDHR1, data = as.data.frame(t(df_patient_mut_b)))
ggsurvplot(fit_mut, data = as.data.frame(t(df_patient_mut_b)), pval = TRUE)