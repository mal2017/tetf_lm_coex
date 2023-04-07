library(tidyverse)
library(vroom)

subsample <- T
subsample <- snakemake@params$subsample

#dat <- vroom("results/linear_models/male/0/expression.tsv.gz", col_select = "feature")
dat <- vroom(snakemake@input[["expression"]])

#split_invoc <- "split -e -d -a 4 -n r/8 - chunk_"
split_invoc <- snakemake@params[["split_invoc"]]

tes <- dat %>% filter(!str_detect(feature,"FBgn")) %>% pull(feature)
genes <- dat %>% filter(str_detect(feature,"FBgn")) %>% pull(feature)

if (subsample) {
  set.seed(2023)
  genes <- sample(genes,10)
}

res <- expand_grid(tes = tes, genes=genes)

vroom_write(res, pipe(split_invoc), col_names = F)

