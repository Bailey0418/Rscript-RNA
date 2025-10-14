##### R #####
#读取数据
library(readr)
data1 <- read_tsv("exp_seq.RECA-EU.tsv")
data2 <- read_tsv("sample.RECA-EU.tsv")
data3 <- read_tsv("donor.RECA-EU.tsv")
head(data1)
head(data2)
head(data3)
#整理数据
library(dplyr)
library(tidyr)
expr <- data1 %>%
  dplyr::select(icgc_sample_id, gene_id, raw_read_count)
expr_wide <- expr %>%
  tidyr::pivot_wider(names_from = icgc_sample_id,
                     values_from = raw_read_count)
expr_mat <- as.data.frame(expr_wide)
rownames(expr_mat) <- expr_mat$gene_id
expr_mat <- expr_mat[, -1] 

meta <- data2 %>%
  dplyr::select(icgc_sample_id, submitted_sample_id, icgc_donor_id, project_code) %>%
  mutate(
    sample_type = ifelse(grepl("-RT-", meta$submitted_sample_id), "Tumor",
                         ifelse(grepl("-RA-", meta$submitted_sample_id), "Normal",
                                ifelse(grepl("-RN-", meta$submitted_sample_id), "Normal",
                                       ifelse(grepl("-ND-", meta$submitted_sample_id), "Germline", NA))))
  )
meta <- data2 %>%
  dplyr::select(icgc_sample_id, submitted_sample_id, icgc_donor_id, project_code) %>%
  mutate(
    sample_type = ifelse(grepl("(N$|RA$)", as.character(submitted_sample_id)), 
                         "Normal", 
                         "Tumor")
  )
head(meta)
table(meta$sample_type)

common_samples <- intersect(colnames(expr_mat), meta$icgc_sample_id)
expr_mat <- expr_mat[, common_samples]
meta1 <- meta[match(common_samples, meta$icgc_sample_id), ]
table(meta1$sample_type)

library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = mart
)
head(genes_annot)

expr_mat$ensembl_gene_id <- rownames(expr_mat)
expr_annot <- merge(genes_annot, expr_mat, by = "ensembl_gene_id")
expr_annot <- expr_annot[expr_annot$hgnc_symbol != "", ]
expr_annot <- expr_annot[!duplicated(expr_annot$hgnc_symbol), ]
rownames(expr_annot) <- expr_annot$hgnc_symbol
expr_annot <- expr_annot[, !(names(expr_annot) %in% c("ensembl_gene_id", "hgnc_symbol"))]
dim(expr_annot)
expr_annot[1:5, 1:5]
write.csv(expr_annot, "RECA_EU_expression_symbol.csv")

library(dplyr)
survival_info <- data3 %>%
  dplyr::select(
    icgc_donor_id,
    donor_vital_status,
    donor_survival_time,
    donor_interval_of_last_followup,
    donor_age_at_diagnosis,
    donor_sex,
    disease_status_last_followup
  ) %>%
  mutate(
    # vital_status: 1=dead, 0=alive
    vital_status = case_when(
      donor_vital_status == "deceased" ~ 1,
      donor_vital_status == "alive" ~ 0,
      TRUE ~ NA_real_
    ),
    # survival_time：如果死亡则用 donor_survival_time，否则用最后随访时间
    survival_time = case_when(
      !is.na(donor_survival_time) ~ as.numeric(donor_survival_time),
      !is.na(donor_interval_of_last_followup) ~ as.numeric(donor_interval_of_last_followup),
      TRUE ~ NA_real_
    )
  )

meta_full <- meta1 %>%
  left_join(survival_info, by = "icgc_donor_id")
head(meta_full)
table(meta_full$sample_type, useNA = "ifany")
table(meta_full$vital_status, useNA = "ifany")
write.csv(meta_full, "RECA_EU_sample_metadata.csv")
