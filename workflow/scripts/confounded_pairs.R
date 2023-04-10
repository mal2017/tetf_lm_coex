library(tidyverse)

#fixed_overlaps <- read_tsv("../tetf_tidal/results/overlaps/excluded_gene_te_pairs.tsv")
fixed_overlaps <- read_tsv(snakemake@input$fixed)

# we are checking if 'y' in the 'edges' df overlaps
# tes, so we rename name.x as y to do the join
fixed_overlaps <- fixed_overlaps %>% dplyr::select(y=name.x, overlaps = name.y) %>%
  mutate(filter_type = "direct.fixed.overlap")

#edges <- read_rds("results/linear_models/female.gene_2_gene_edges.rds")
edges <- read_rds(snakemake@input$edges)

# x = a gene
# y = a gene correlatd with x
# overlaps = a feature that overlaps y,
# so that gene/te pair x/overlaps 
# can be excluded or monitored downstream
res0 <- edges %>% left_join(dplyr::select(fixed_overlaps,-filter_type)) %>% drop_na()

# combine pairs of genes (x) and tes (y) that
# are excluded due to directly overlap with fixed insertion 
# or fixed insertion overlapping with a highly correlated gene
res <- res0 %>%
  dplyr::select(x, y=overlaps) %>%
  distinct() %>%
  mutate(filter_type = "correlated.gene.fixed.overlap") %>%
  bind_rows(dplyr::rename(fixed_overlaps,x="y",y=overlaps)) %>%
  filter(str_detect(x,"FBgn") & !str_detect(y,"FBgn"))
  
write_tsv(res, snakemake@output$tsv)