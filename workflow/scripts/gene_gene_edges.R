library(tidyverse)
library(corrr)

#current_model <- "female_model_01"
current_model <- snakemake@params[["model_id"]]

cr_fl <- snakemake@input[["cor"]]

# get both triangles of corr matrix in long format --------------------------------

cr <- read_rds(cr_fl)

# make the square matrix long
cr <- cr %>% 
  pivot_longer(-term,names_to = "y",values_to = "r") %>% 
  dplyr::rename(x="term") %>% 
  drop_na()

# the above data is only 1 triangle for space reasons - here we reconstitute
cr <- cr %>% select(x=y, y=x, r) %>% bind_rows(cr)

# establish cutoffs for calling a gene coexpressed with another gene -----------
gg <- read_tsv("http://ftp.flybase.net/releases/FB2022_05/precomputed_files/genes/gene_group_data_fb_2022_05.tsv.gz",skip = 8) %>%
  filter(FB_group_symbol %in% c("RPL")) %>%
  pull(Group_member_FB_gene_id)

# the median in-group correlation of a a plausibly truly coregulated set of genes
cutoff <- cr %>%
  filter(x %in% gg & y %in% gg) %>%
  pull(r) %>%
  abs() %>%
  median()

# get 'valid' edges and a the background set
edges <- cr %>% filter(abs(r) > cutoff & x!=y)
nudges <- cr %>% filter(abs(r) < cutoff & x!=y & x %in% edges$x)

write_rds(edges,snakemake@output[["edges_rds"]])
write_rds(nudges,snakemake@output[["nudges_rds"]])