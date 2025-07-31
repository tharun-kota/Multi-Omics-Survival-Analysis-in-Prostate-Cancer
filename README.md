# Prostate Cancer Survival Analysis Using Multi-Omics Data

## Overview

This project investigates molecular heterogeneity in prostate cancer by integrating mRNA expression data, somatic mutation profiles, and clinical survival information. The study applies both unsupervised and supervised machine learning techniques to identify prognostic biomarkers and explore potential molecular subtypes.

## Objectives

- Explore the relationship between gene expression, mutation profiles, and patient survival.
- Identify survival-associated molecular features using regularized Cox regression.
- Investigate whether gene expression patterns reflect somatic mutation signatures.
- Examine pathway-level activity using PROGENy for biological interpretability.

## Data Source

Data was obtained from the [cBioPortal](https://www.cbioportal.org/) platform and includes:

- mRNA-seq expression data (TPM z-scores)
- Somatic mutation profiles
- Clinical survival data (overall survival in months and status)

---

## Methodology

### 1. Data Acquisition and Preprocessing

- Filtered clinical data for valid survival times and event statuses.
- Mapped Entrez IDs to gene symbols and cleaned mRNA expression data.
- Converted mutation data into a binary matrix (mutated vs. wild-type).
- Transformed datasets into patient × gene matrices.

### 2. Survival Analysis

- Constructed `Surv` objects for both omics datasets.
- Applied Elastic Net regularized Cox regression using `glmnet`.
- Used cross-validation to optimize model hyperparameters (`alpha`, `lambda`).
- Selected features with non-zero coefficients as survival-associated genes/mutations.

### 3. Interpretation and Visualization

- Used `coxph` for multivariate survival modeling.
- Visualized Kaplan–Meier survival curves using `ggsurvplot`.
- Inferred pathway activity using PROGENy.
- Conducted correlation and log-rank analysis to evaluate pathway relevance.

---

## Results

### mRNA-Seq Analysis

- **29** genes were selected via Elastic Net, **7** of which were statistically significant for survival.
- MAPK pathway activity was identified as a key stratifier for survival subgroups.
- Other pathways (e.g., NFkB, TNFa) showed strong co-regulation.

### Mutation Analysis

- **57** mutation-associated genes identified; **SOX9** and **CDH2** showed significant survival impact.
- These genes are involved in epithelial–mesenchymal transition (EMT) and tumor aggressiveness.

---

## Key Findings

- MAPK signaling may serve as a prognostic marker in prostate cancer.
- Mutations in SOX9 and CDH2 are potentially clinically relevant.
- No strong alignment between gene expression and mutation layers was observed in survival context.
- Pathways like PI3K/AKT and TP53 did not reach statistical significance, possibly due to cohort limitations.

---

## Conclusion

This project highlights the value of multi-omics integration in personalized oncology. While further validation is needed, the results support MAPK activity and specific gene mutations (e.g., SOX9, CDH2) as potential biomarkers. Future directions may include incorporation of additional omics layers such as proteomics or epigenomics to further refine patient stratification.

---

## Tools & Libraries

- R (`glmnet`, `survival`, `org.Hs.eg.db`, `survminer`, `PROGENy`, `corrplot`)
- Data transformation using `tidyverse`

---


