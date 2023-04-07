library(tidyverse)

x <- 0
sex <- "male"
anov_fl <- sprintf("results/linear_models/%s/%s/lm.anova.tsv.gz",sex,x)
rainb_fl <- sprintf("results/linear_models/%s/%s/lm.rainbow.tsv.gz",sex,x)
breusch_fl <- sprintf("results/linear_models/%s/%s/lm.breusch_pagan.tsv.gz",sex,x)
glance_fl <- sprintf("results/linear_models/%s/%s/lm.glance.tsv.gz",sex,x)
tidy_correct_fl <- sprintf("results/linear_models/%s/%s/lm.tidy.corrected.tsv.gz",sex,x)

anov_fl <- snakemake@input[["anova"]]
rainb_fl <- snakemake@input[["rainbow"]]
breusch_fl <- snakemake@input[["breusch_pagan"]]
glance_fl <- snakemake@input[["glance"]]
tidy_correct_fl <- snakemake@input[["tidy_correct"]]

# anova
anov <- read_tsv(anov_fl)

anov <- mutate(anov, rss_full = rss-sumsq) %>%
  dplyr::select(-AIC,-statistic) %>%
  pivot_wider(names_from = term,values_from = c("sumsq","rss","p.value","rss_full","df"),names_sep = "_anova_") %>%
  dplyr::select(feature.x,feature.y,rss_full = `rss_anova_<none>`,contains("sumsq_"),contains("p.value_")) %>%
  dplyr::select(-contains("<none>")) %>%
  mutate(total_variance = rowSums(across(contains("rss") | contains("sumsq")),na.rm = T)) %>%
  mutate(explained_variance = rowSums(across(contains("sumsq")),na.rm = T)/total_variance)

# rainbow
rainb <- read_tsv(rainb_fl)

# this fails when running tests/sanity checks on just a few samples
if (nrow(rainb) >0 ) {
  rainb <- dplyr::select(rainb,feature.x,feature.y,p.value_rainbow=p.value)
} else {
  rainb <- dplyr::select(anov, feature.x,feature.y) %>% mutate(p.value_rainbow = 1)
}

# breusch pagan
breusch <- read_tsv(breusch_fl) %>%
  dplyr::select(feature.x,feature.y,p.value_breuschpagan=p.value)

# glance
glance <- read_tsv(glance_fl) %>%
  dplyr::select(feature.x,feature.y,ftest_r2 = r.squared,ftest_adjr2 = adj.r.squared,p.value_ftest_r2=p.value)

# corrected
tidy.corrected <- read_tsv(tidy_correct_fl) %>%
  dplyr::select(feature.x,feature.y,estimate,estimate.qnorm,p.value_lm_coef=p.value,mean_qtile.gene,mean_qtile.te)

res <- tidy.corrected %>%
  left_join(breusch) %>%
  left_join(rainb) %>%
  left_join(glance) %>%
  left_join(anov)

res <- res %>%
    mutate(model = snakemake@params[["model"]],
        rep = snakemake@params[["rep"]]) %>%
    relocate(model, rep)


write_tsv(res, snakemake@output[["tsv"]])