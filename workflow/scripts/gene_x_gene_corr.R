library(corrr)
library(tidyverse)
library(vroom)
library(spqn)
library(SummarizedExperiment)

#expression_fl <- "results/linear_models/female/0/expression.tsv.gz"
expression_fl <- snakemake@input[["expr"]]

dat <- vroom(expression_fl)

dat <- dat %>%
  filter(str_detect(feature,"FBgn"))

corr_df <- dat %>% 
  pivot_longer(-feature,names_to = "strain") %>%
  pivot_wider(names_from = feature, values_from = value) %>%
  correlate(diagonal = 1)

cor_m <- as_matrix(corr_df)
expr_m <- column_to_rownames(dat,"feature")[colnames(cor_m),] %>% as.matrix()

#plot_signal_condition_exp(cor_m, rowMeans(expr_m), signal=0.01)

cor_m_spqn <- normalize_correlation(cor_m, ave_exp=rowMeans(expr_m), ngrp=20, size_grp=300, ref_grp=10)

#plot_signal_condition_exp(cor_m_spqn, rowMeans(expr_m), signal=0.01)

dimnames(cor_m_spqn) <- dimnames(cor_m)

res <- as_cordf(cor_m_spqn)

write_rds(res,snakemake@output[["rds"]])