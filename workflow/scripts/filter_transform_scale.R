library(tidyverse)
library(SummarizedExperiment)
library(rlang)

set.seed(1)

#sefile <- "../tetf_expression_quant/results/quantification/vanilla_salmon_tes_transcripts/se.gene.0.rds"
sefile <- snakemake@input[["sefile"]]
se <- read_rds(sefile)

# import counts just for filtering on count scale,
# these aren't used for lms
x <- assay(se,"counts")

# set up filter
#filt <- (((!str_detect(rownames(x),'FBgn')) & rowSums(x > 1) > 10) | (rowSums(x > 10) > 100)) & rowSums(x == 0) < 0.3*ncol(x)
filt <- ifelse(exists("snakemake"), snakemake@params[["filt"]], "T")
print(filt)
filt <- parse_expr(filt)
features_2_use <-  rownames(x[eval(filt),])

rm(x); #rm(se)

#matfile <- "results/quantification/vanilla_salmon_tes_transcripts/female.gene.per_feature.vst.0.tsv.gz"
matfile <- snakemake@input[["mat"]]
x <- vroom::vroom(matfile,num_threads = 1)

# remove NA values - this is most common in sex biased genes
# because i tore these in a big matrix with males and females
x <- x %>% mutate(across(where(is.numeric), replace_na, 0))

# ------------------------- perform filtering -------------------------
x <- column_to_rownames(x, "feature")[features_2_use,]

#any(is.na(x))

# --------------------- perform transforms ---------------------------
#transforms <- "x"
transforms <- snakemake@params[["transforms"]]
transforms <- parse_expr(transforms)

x <- eval(transforms)

# ------------------------- scaling ----------------------------------

print(x[1:5,1:4])

# scale <- T
scale <- snakemake@params[["scale"]]
message(class(x))
message(scale)
if (scale) {
  x <- t(scale(t(x)))
}

print(x[1:5,1:4])

x <- x[!rowAnyNAs(x),]

print(x[1:5,1:4])

# ------------------------- pcs -------------------------------------
#pcs <- 4
# See supplementary materials of:
# Co-expression patterns define epigenetic regulators associated with neurological dysfunction
# Boukas et al. Genome Biology.
pcs <- snakemake@params[["pcs"]]



if (pcs > 0) {
  message(pcs)
  x <- t(WGCNA::removePrincipalComponents(t(x), min(pcs,ncol(x)-1)))
}

print(x[1:5,1:4])

# ------------------------ export ----------------------------------

x <- as_tibble(x, rownames = "feature")

x <- arrange(x, -str_detect(feature,"FBgn"))
message(class(x))

write_tsv(x,file = snakemake@output[["mat"]])
#vroom::vroom(x,snakemake@output[["mat"]])
