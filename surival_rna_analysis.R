getwd()
setwd("/Users/tharun kota/OneDrive/Desktop/UCD/Trimester_3/AI_Personalised_Medicine/Assignment_1/")
install.packages("survminer")

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library("org.Hs.eg.db")

clin_raw <- read.delim("prad_tcga_gdc/prad_tcga_gdc/data_clinical_patient.txt", sep = '\t',skip = 4)
rownames(clin_raw) <- clin_raw$PATIENT_ID

RNA_raw <- read.delim("prad_tcga_gdc/prad_tcga_gdc/data_mrna_seq_tpm_zscores_ref_all_samples.txt",check.names = FALSE)

#BiocManager::install("org.Hs.eg.db")


# Convert Entrez IDs to character
entrez_ids <- as.character(RNA_raw$Entrez_Gene_Id)
# Map Entrez → HGNC Symbol
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = entrez_ids,
  column = "SYMBOL",   # what you want
  keytype = "ENTREZID",  # what you have
)
RNA_raw$Hugo_symbol <- gene_symbols
RNA_raw <- RNA_raw[, c(1, ncol(RNA_raw), 2:(ncol(RNA_raw)-1))]
RNA_raw <- RNA_raw[!is.na(RNA_raw$Hugo_symbol), ]
#RNA_raw <- RNA_raw[RNA_raw$Hugo_symbol!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Hugo_symbol),]
rownames(RNA_raw) <- RNA_raw$Hugo_symbol
RNA <- as.data.frame(t(RNA_raw[, !colnames(RNA_raw) %in% c("Hugo_symbol", "Entrez_Gene_Id")]))
#RNA <- as.data.frame(t(RNA_raw[-1:-2]))
RNA <- cbind(Trim_ID = stringr::str_sub(rownames(RNA), end = -5), RNA)
#RNA$Trim_ID <- stringr::str_sub(rownames(RNA), end = -5)
# Align clinical data:
clin <- clin_raw[RNA$Trim_ID, ]
#clin <- clin_raw[str_sub(row.names(RNA), end = -4),]
RNA <- RNA[!duplicated(RNA$Trim_ID), ]
clin <- clin[!duplicated(clin$PATIENT_ID), ]
rownames(RNA) <- RNA[[1]]; RNA <- RNA[, -1]
RNA <- RNA[, colSums(is.na(RNA)) == 0]


# create a survival object consisting of times & censoring
surv_obj <- Surv(time = clin$OS_MONTHS, 
                 event = clin$OS_STATUS=="1:DECEASED")

fit <- survfit(surv_obj ~ 1, data = clin)
ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival")


# fit multivariate  COX proportional hazard model for some specific genes
fit.coxph <- coxph(surv_obj ~ KRAS + TP53, 
                   data = RNA)
ggforest(fit.coxph, data = RNA)

library(progeny)
zscores = as.matrix(t(RNA))
pathways <- progeny(zscores, scale=TRUE, organism="Human")  #, top = 100, perm = 1)
path_df = as.data.frame(pathways)

# fit multivariate  COX proportional hazard model for some specific paythways
fit.coxph <- coxph(surv_obj ~ Androgen+EGFR+Estrogen+Hypoxia
                   +MAPK+NFkB+p53+PI3K+TGFb+TNFa+Trail+
                     VEGF+WNT+`JAK-STAT`, data = path_df)
ggforest(fit.coxph, data = path_df)

library(corrplot)
corrplot(cor(path_df), type = "upper", order = "hclust",tl.col = "black", tl.srt = 45)


# analysing on what MAPK activity depend
MAPK_df = cbind(MAPK = path_df$MAPK, PI3K = path_df$PI3K,
                KSR1 = RNA$KSR1, KSR2 = RNA$KSR2, IQGAP1 = RNA$IQGAP1, 
                IQGAP2 = RNA$IQGAP2, IQGAP3 = RNA$IQGAP3, GAB1 = RNA$GAB1, GAB2 = RNA$GAB2,
                KRAS = RNA$KRAS, NRAS = RNA$NRAS, HRAS = RNA$HRAS, BRAF=RNA$BRAF, 
                CRAF=RNA$RAF1, ARAF=RNA$ARAF)
corrplot(cor(MAPK_df), type = "upper", order = "hclust",tl.col = "black", tl.srt = 45)

# now creating object without zero times
clin_filt <- clin[clin$OS_MONTHS > 0,]
RNA_filt <- RNA[clin$OS_MONTHS > 0,]
path_filt <- path_df[clin$OS_MONTHS > 0,]
# create a survival object consisting of times & censoring
surv_filt <- Surv(time = clin_filt$OS_MONTHS, 
                  event = clin_filt$OS_STATUS=="1:DECEASED")
fit <- survfit(surv_filt ~ 1, data = clin_filt)
ggsurvplot(fit, data = clin_filt, xlab = "Month", ylab = "Overall survival")

# glmnet
library("glmpath")
library("glmnet")
library("penalized")
library("glmnetUtils")

# finding optimal alpha
lena = 11
cvafit <- cva.glmnet(as.matrix(RNA_filt),surv_filt,alpha=seq(0, 1, len = lena),nfolds=10,family="cox",type.measure = "C")
# plotting results of alpha scan
cvafit$cvms <- NULL
cvafit$cvss <- NULL
for (i in 1:lena) {
  cvafit$cvms <- c(cvafit$cvms, max(cvafit$modlist[[i]]$cvm))
  cvafit$cvss <- c(cvafit$cvss, max(cvafit$modlist[[i]]$cvs))
}
plot(seq(0, 1, len = lena),cvafit$cvms,type = 'b',xlab='alpha',ylab='C-index',ylim=c(0.4,0.8))
arrows(x0=seq(0, 1, len = lena), y0=cvafit$cvms-cvafit$cvss, x1=seq(0, 1, len = lena), y1=cvafit$cvms+cvafit$cvss, code=3, angle=90, length=0.1)

# chosing optimal alpha
opt_i <- which.max(cvafit$cvms)
opt_a <- cvafit$alpha[opt_i]
# we can select alpha manually instead
opt_a <- 0.3

# CV glmnet
cvfit <- cv.glmnet(as.matrix(RNA_filt),surv_filt,family="cox",type.measure = "C",alpha = opt_a)
plot(cvfit)
# selecting optimal lambda
opt_i <- which.max(cvfit$cvm)
opt_l <- cvfit$lambda[opt_i]

# analysing results
cfs = coef(cvfit,s=opt_l)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
write.csv(meaning_vals,file='survival_coeffs_cvglmnet_RNAseq.csv')
print(paste(meaning_coefs,collapse=" + "))


# fit multivariate  COX proportional hazard model for some specific genes
fit.coxph <- coxph(surv_filt ~ HOXB6 + SIGLEC16 + PRR27 + `HOXB-AS3` + MAGEA5P + `KRTAP5-2` + ELOA2 + LHX5 + GAGE2A + `RAPGEF4-AS1`,
                   data = RNA_filt)
ggforest(fit.coxph, data = RNA)



# glmnet over progeny
fit_glm <- glmnet(path_filt,surv_filt,family="cox") # , alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)
# analysing results
cfs = coef(fit_glm,s=0.0002)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
write.csv(meaning_vals,file='survival_coeffs_progeny_glmnet.csv')
print(paste(meaning_coefs,collapse=" + "))

# plotting Kaplan-Mayer curves
pathway = 'MAPK'
pathway_data = path_filt$MAPK
# sort age 
uni_path = sort(unique(pathway_data))
# store results
results_path = matrix(1,length(uni_path))
# do a for loop for every unique value of age mat
for (i in 2:(length(uni_path)-1)){ # Starting from 2 because the first element would yield no elements higher than.
  path_i = 1*(pathway_data>uni_path[i])
  # survdiff is the function from the survival package 
  logrank = survdiff(surv_filt ~ path_i)
  # store in results_age
  results_path[i] = logrank$pvalue
}

# Plot unique elements of age against p-value
plot(uni_path, results_path, log = 'y')

# Select minimum P-value
min_p_path = which.min(results_path)
# here are 2 good thresholds, -1 and 1
opt_thr = uni_path[min_p_path]
# I recalculated the P-value as a sanity check
pval = survdiff(surv_filt ~ pathway_data>opt_thr)$pvalue
nplus = sum(pathway_data>opt_thr)   # how many patients we have in high group
nminus = sum(pathway_data<opt_thr)   # how many patients we have in low group


# Create MAPK activity groups: "high" if ≥ 2, "low" otherwise
path_filt <- path_filt %>% mutate(MAPK_group = ifelse(MAPK >= 2, "high", "low"))

# Count patients in each group
nhigh <- sum(path_filt$MAPK_group == "high")
nlow <- sum(path_filt$MAPK_group == "low")

# Fit Kaplan-Meier model
KM <- survfit(surv_filt ~ MAPK_group, data = path_filt)

# Plot KM survival curves
p <- ggsurvplot(KM, data = path_filt, pval = TRUE, xlab = "Overall survival time, months",
                legend.labs = c(paste("High MAPK activity,\n", nhigh, " patients", sep = ""),
                                paste("Low MAPK activity,\n", nlow, " patients", sep = "")),
                legend.title = "")

# Enhance plot appearance
ggpar(p, 
      font.main = c(13, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.caption = c(16, "bold"), 
      font.legend = c(13, "bold"), 
      font.tickslab = c(14, "bold"))































