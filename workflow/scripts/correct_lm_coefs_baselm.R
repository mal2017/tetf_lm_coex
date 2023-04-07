library(tidyverse)
library(preprocessCore)

#coefs_fl <- "results/linear_models/female_model_01/lm.tidy.tsv.gz"
coefs_fl <- snakemake@input[["coefs"]]

coefs.all <- read_tsv(coefs_fl) %>%
  mutate(relationship = ifelse(sign(estimate) == 1, "pos", "neg")) %>%
  mutate(significant = p.value < 0.05)


#expression_fl <- "results/quantification/vanilla_salmon_tes_transcripts/female.gene.per_feature.vst.tsv.gz"
expression_fl <- snakemake@input[["expression"]]

expression <- read_tsv(expression_fl) %>%
  mutate(type = ifelse(str_detect(feature,"FBgn"),"gene","TE"))

mean.expr <- expression %>%
  mutate(., mean = rowMeans(select(.,!contains("feature") &  !contains("type")))) %>%
  select(feature,type, contains("mean"))

pairs_filter <- coefs.all %>%
  dplyr::select(feature.x,feature.y) %>%
  distinct() %>%
  pivot_longer(c("feature.x","feature.y"),names_to = "feature.side", values_to = "feature") %>%
  dplyr::select(-feature.side) %>%
  distinct()

n_qtiles <- 10

# only retain values actually used in the modeing
mean.expr.filt <- mean.expr %>%
  semi_join(pairs_filter, by=c("feature"))

expr.qtiles <- mean.expr.filt %>%
  group_by(type) %>%
  mutate(across(contains("mean"),ntile,n_qtiles, .names="{.col}_qtile"),.keep="unused") %>%
  ungroup()

coefs.qtiles <- coefs.all %>%
  filter(term == "x") %>%
  left_join(expr.qtiles,by=c("feature.x"="feature")) %>%
  left_join(expr.qtiles,by=c("feature.y"="feature"),suffix=c(".gene",".te"))

#coefs.qtiles %>%
#  ggplot(aes(as.factor(mean_qtile.gene),estimate)) +
#  geom_boxplot(aes(middle=mean(estimate)),outlier.shape = NA) +
#  xlab("gene expression qtile") +
#  ggpubr::stat_compare_means(method = 'anova',label.y = 70) +
#  coord_cartesian(ylim=c(-75,75)) +
#  theme(aspect.ratio = 1)

message("finding reference pairs...")
coefs.qtiles.refs <- coefs.qtiles %>%
  filter(abs(estimate) < 100) %>%
  filter(mean_qtile.gene %in% c(5,6)) %>%
  group_by(relationship) %>%
  summarise(ref = list(estimate)) %>%
  ungroup()

message("setting up data for correction...")
corrected0 <- coefs.qtiles %>%
  left_join(coefs.qtiles.refs) %>%
  group_by(mean_qtile.te, mean_qtile.gene, relationship)

message("performing correction...")
corrected <- mutate(corrected0, estimate.qnorm = normalize.quantiles.use.target(matrix(estimate,ncol = 1),ref[[1]])[,1]) %>% 
  ungroup()

#library(ggridges)
#corrected %>%
#  mutate(across(contains("qtile"),as.factor)) %>%
#  ggplot() +
#  geom_density_ridges2(aes(y=mean_qtile.gene,x=estimate.qnorm,fill=mean_qtile.gene),color=NA) +
#  scale_fill_cyclical(values = c("darkgray", "lightgray")) +
#  coord_cartesian(xlim=c(-50,50))

write_tsv(corrected,snakemake@output[["tsv"]])
