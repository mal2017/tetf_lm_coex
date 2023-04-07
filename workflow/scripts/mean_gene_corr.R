library(DescTools)
library(tidyverse)
library(corrr)

fls <- c("results/linear_models/female/0/gene_x_gene.corrr.rds")

fls <- snakemake@input

df <- map(fls, read_rds)

is_test <- ifelse(exists("snakemake"),snakemake@params$test,T)

if ( is_test | interactive()) {
  df <- map(df,~{.x[1:1000,1:1001]}) %>% map_df(stretch, .id = "rep")  
} else {
  df <- map_df(df,stretch, .id = "rep")  
}

res <- df %>% mutate(fishz = FisherZ(r))

rm(df);gc()

res <- res %>%
  group_by(x, y) %>%
  summarise(mean_fishz = mean(fishz), .groups = "drop")

res <-  mutate(res,r=FisherZInv(mean_fishz))

# handle unlikely case of perfect correlations between genes.
res <- res %>% mutate(r = if_else(is.nan(r) & is.infinite(mean_fishz), sign(mean_fishz) * 1,r))

res <- res %>% dplyr::select(x,y,r)

# was long tbl, now widen
res <- retract(res,x,y,r) 

# from corrr pkg
res <- as_cordf(res)

# lower tri only
res <- shave(res)

write_rds(res, snakemake@output[["rds"]])