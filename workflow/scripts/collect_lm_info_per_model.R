library(tidyverse)

#confounded_pairs_fls <- "results/linear_models/female.confounded_pairs.tsv"
confounded_pairs_fls <- snakemake@input$confounded_pairs

# don't want to duplicate rows in the filter
confounded_pairs <- read_tsv(confounded_pairs_fls) %>%
  group_by(x,y) %>%
  summarise(overlap_filter_type = paste(sort(filter_type),collapse="+"), .groups = "drop")
  
#info_fls <- Sys.glob("results/linear_models/female/0/lm.collected-info.tsv.gz")
info_fls <- snakemake@input[["info"]]

info <- info_fls %>% map_df(read_tsv)

# annotate as valid if it passes breusch pagan and rainbow tests 
info <- info %>%
  mutate(valid = p.value_breuschpagan > 0.05 & p.value_rainbow > 0.05)
  
# note as reproducible if all reps suggest nonsignificance or all reps are
# significant and have same signs for gene expression coefs
info <- info %>%
  group_by(model,feature.x,feature.y) %>%
  mutate(significant_model_0 = all(p.value_ftest_r2 < 0.05)) %>%
  mutate(reproducible = all(valid) & 
           (
             all(!significant_model_0) | 
               (all(significant_model_0) & length(unique(sign(estimate)))==1)
             )
         )

# for each gene/te pair, take the best fitting model overall to be safe with ties broen by highest pvalue for the gene coef, to be conservative
info <- info %>%
  group_by(model,feature.x,feature.y) %>%
  slice_min(p.value_ftest_r2, n=1) %>%
  slice_max(p.value_anova_x,n=1,with_ties = F)
  
# now adjust term p values once we have 1 model per te/gene pair
# note that I also previously adjusted f test p values above,
# but I think it makes more sense to do this after
# finding representative rep for each TE/gene pair such that the adjustment
# is performed on p vals actually in the results table
info <- info %>%  
  group_by(model) %>%
  mutate(across(starts_with("p.value_anova"),
                ~p.adjust(.x,method="BH"),
                .names= "adj_{.col}")) %>%
  mutate(adj_p.value_ftest_r2 = p.adjust(p.value_ftest_r2,method="BH")) %>%
  ungroup() %>%
  dplyr::select(-significant_model_0)

# add info about pairs with fully fixed overlaps
info <- info %>% left_join(confounded_pairs, by=c(feature.x="x",feature.y="y"))

info <- info %>%
  mutate(significant_model = adj_p.value_ftest_r2 < 0.1 & reproducible & valid) %>%
  mutate(significant_x = adj_p.value_anova_x < 0.1 & significant_model & is.na(overlap_filter_type))

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

info <- info %>% 
  group_by(model) %>%
  mutate(coef.quantile = cume_dist(abs(estimate.qnorm))) %>%
  ungroup() %>%
  left_join(lkup, by=c(feature.x = "gene_ID")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol),feature.x,gene_symbol)) %>%
  relocate(gene_symbol, .after = "feature.x")


write_tsv(info,snakemake@output[["tsv"]])