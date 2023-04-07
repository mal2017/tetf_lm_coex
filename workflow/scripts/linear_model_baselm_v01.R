library(tidyverse)
library(lmtest)

# -------- imort pairs for the chunk of lms to be evaluated --------
#pairs_fl <- "results/linear_models/female/0/chunk_0000"
pairs_fl <- snakemake@input[["chunk"]]

pairs <- vroom::vroom(pairs_fl,col_names = c("feature.y","feature.x")) #%>% head(10)

# -------- prep data for the chunk of lms to be evaluated --------

#dat_fl <- "results/linear_models/female/0/expression.tsv.gz"
dat_fl <- snakemake@input[["dat"]]
dat <- vroom::vroom(dat_fl,num_threads = 1)

#cd_fl <- "../tetf_expression_quant/results/meta/metadata.csv"
cd_fl <-  snakemake@input[["cd"]]
cd <- read_csv(cd_fl)

dat <- pivot_longer(dat,-feature,names_to = "sample",values_to = "score")

# ----------- join with metadata and pairs to create working dataframe -----------
df <- left_join(pairs,dat, by=c(feature.y="feature"))

df <- left_join(df,dat, by=c(feature.x="feature",sample="sample"), suffix=c(".y",".x"))

df <- left_join(df,cd, by=c(sample="sample_name"))

# -------------- add the overlap annotation to the dataframe ---------------
message("adding main gene overlap annotation...")
#ol_fl <- "../tetf_tidal/results/overlaps/overlaps.tsv.gz"
ol_fl <- snakemake@input[["ol"]]

ol <- read_tsv(ol_fl) %>%
  dplyr::select(Strain,feature.x=name.x,feature.y=name.y,overlap) %>%
  filter(str_detect(Strain,"DGRP") & feature.x %in% df$feature.x & feature.y %in% df$feature.y)

df <- df %>% left_join(ol, by = c("feature.y", "feature.x", "Strain")) %>%
  mutate(overlap = ifelse(is.na(overlap),F,overlap))

# ------------- add 'coexpressed gene' overlap anno to dataframe -----------
message("adding coexpressed gene overlap annotation...")
#edges_fl <- "results/linear_models/female.gene_2_gene_edges.rds"
edges_fl <- snakemake@input[["edges"]]

edges <- read_rds(edges_fl) %>%
  filter(x %in% df$feature.x)

ol2_coex <- edges %>% 
  dplyr::rename(coex.gene=y) %>%
  left_join(dplyr::select(ol,Strain,x=feature.x,y=feature.y,overlap.coex.gene=overlap),by=c("coex.gene"="x")) %>%
  group_by(Strain,x,y) %>%
  summarize(overlap.coex.gene = any(overlap.coex.gene),.groups="drop") %>%
  dplyr::rename(feature.x=x,feature.y=y)

df <- df %>% 
  left_join(ol2_coex, by = c("feature.y", "feature.x", "Strain")) %>%
  mutate(overlap.coex.gene = ifelse(is.na(overlap.coex.gene),F,overlap.coex.gene))


rm(ol);rm(ol2_coex);rm(edges);gc()

# ----------- add the copies annotation to the dataframe
message("annotating with estimated copy number info...")
#copies_fl <- "../tetf_dgrp_wgs/results/copies/copies.tsv"
copies_fl <- snakemake@input[["copies"]]

copies <- read_tsv(copies_fl) %>%
  group_by(sequence) %>%
  mutate(scaled.copies = scale(est.copies)[,1]) %>%
  ungroup()

#unique(copies$sequence)[!unique(copies$sequence) %in% df$feature.y]
copies <- copies %>% mutate(sequence = str_replace(sequence,"\\\\","_")) #%>% filter(str_detect(sequence,"Ddip"))

stopifnot(all(df$feature.y %in% copies$sequence))

# join estimated copies and set to 1 where no info is available, (0 for scaled)
df <- df %>%
  left_join(copies, by=c(Strain="Strain",feature.x="sequence")) %>%
  left_join(copies, by=c(Strain="Strain",feature.y="sequence")) %>%
  mutate_at(c("est.copies.y","est.copies.x"), replace_na, 1) %>%
  mutate_at(c("scaled.copies.x","scaled.copies.y"),replace_na,0)

# ----------------remove outliers per feature by mads ----------------------
message("filtering outliers...")
#mads_filter <- 3
mads_filter <- snakemake@params[["mads_filter"]]

dat <- dat %>% 
  group_by(feature) %>% 
  mutate(mads = abs(score- median(score))/mad(score)) %>%
  ungroup() %>%
  filter(mads < mads_filter)

# --------begin to evaluate the lms --------
message("evaluating LMs...")
# rename these variables so the config can be more succinct
df <- df %>% dplyr::rename(x="score.x", y="score.y")

# main formula
# TESTING: formula <- "y ~ 0 + x + wolbachia + scaled.copies.y + overlap"
# TESTING: formula <- "y ~ 0 + x + wolbachia + scaled.copies.y + overlap + overlap.coex.gene"
formula <- snakemake@params[["formula"]]

# alternate formula; mainly relevant when using the overlap term, which is frequently all F for any
# given TE/gene pair.
# TESTING: alt_formula <- "y ~ 0 + x + wolbachia + scaled.copies.y"
# TESTING: alt_formula <- "y ~ 0 + x + wolbachia + scaled.copies.y + overlap.coex.gene"
alt_formula <- snakemake@params[["alt_formula"]]


# execute the modeling and tests of signficance -------------------------------
possibly_tidy2 <- possibly(broom::tidy,NULL)
possibly_glance2 <- possibly(broom::glance,NULL)
possibly_anova2 <- possibly(~broom::tidy(drop1(.,test="F")),NULL)

# custom func for running lm
run_lm <- function(.x,.y) {
  # alternate formula is used when overlaps are often all False/True, leading to error. 
  #message(.y);
  frm <- as.formula(if_else(length(unique(.x$overlap))==1,alt_formula, formula))
  lm(frm, data = .x)
}

possibly_run_lm <- possibly(run_lm,NULL)

res <- df %>% 
  nest(-c(feature.x,feature.y)) %>% 
  mutate(fit=map2(data,feature.x,possibly_run_lm)) %>%
  mutate(tidy = map(fit,possibly_tidy2),
         glance = map(fit,possibly_glance2),
         anova = map(fit,possibly_anova2))

# tests of LM assumptions -----------------------------------------------------

possibly_bp2 <- possibly(~broom::tidy(bptest(.)),NULL) # heteroskedasticity
# possibly_dw2 <- possibly(~broom::tidy(dwtest(.)),NULL) # dependence of errors, but not meaningful for non-ordered data
possibly_rain2 <- possibly(function(dat,y) broom::tidy(raintest(y,order.by = ~x, data=dat)),NULL) # linearity


res <- res %>%
  mutate(breusch_pagan = map(fit,possibly_bp2),
         rainbow = map2(data,fit,possibly_rain2))


# ------------------- export ------------------------------------------
tidy_fl <- snakemake@output[["tidy"]]
glance_fl <- snakemake@output[["glance"]]
anova_fl <- snakemake@output[["anova"]]
rainbow_fl <- snakemake@output[["rainbow"]]
breusch_pagan_fl <- snakemake@output[["breusch_pagan"]]

res %>% dplyr::select(feature.x, feature.y, tidy) %>%
  unnest(tidy) %>%
  vroom::vroom_write(tidy_fl)

res %>% dplyr::select(feature.x, feature.y, glance) %>%
  unnest(glance) %>%
  vroom::vroom_write(glance_fl)

res %>% dplyr::select(feature.x, feature.y, anova) %>%
  unnest(anova) %>%
  vroom::vroom_write(anova_fl)

res %>% dplyr::select(feature.x, feature.y, rainbow) %>%
  unnest(rainbow) %>%
  vroom::vroom_write(rainbow_fl)

res %>% dplyr::select(feature.x, feature.y, breusch_pagan) %>%
  unnest(breusch_pagan) %>%
  vroom::vroom_write(breusch_pagan_fl)