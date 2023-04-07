library(SummarizedExperiment)
library(tidyverse)

set.seed(1)

#sex <- "male"
sex <- snakemake@wildcards[["sex"]]
sex <- unlist(ifelse(sex == "both",list(c("male","female")),sex))

#se_fl <- "results/quantification/vanilla_salmon_tes_transcripts/se.gene.rds"
se_fl <- snakemake@input[["se"]]
se <- readRDS(se_fl)
se <- se[,se$sex %in% sex]
se <- se[rownames(se)!="GAL4" & !str_detect(rownames(se),"ERCC")]

#ot <- "normcts" # "normcts" "vst" "rlog" "fpkm" "abundance"
ot <- snakemake@wildcards[["expression_unit"]] # "normcts" "vst" "rlog" "fpkm" "abundance"

stopifnot(ot %in% c("counts","normcts","vst","fpkm","abundance","ppc_qn","edaseq_uq","edaseq_qn"))

# ------------ estimate copy number norm ----------------------------------
# do_copy_norm <- T
#do_copy_norm <- snakemake@params[["copy_adjustment"]]
# copies_fl <- "subworkflows/wgs/results/copies/copies.tsv"
#copies_fl <- snakemake@input[["copies"]]

# if (do_copy_norm) {
#  message("PERMFORMING COPY NUM NORM FOR TES")
#  # construct estimated copy matrix

#  copies <- read_tsv(copies_fl) %>%
#    group_by(Strain,sequence) %>%
#    summarise(est.copies = median(est.copies), .groups = "drop") %>%
#    filter(!sequence %in% c("2L","3L","2R","3R","4","X","Y"))

#  stopifnot(nrow(copies %>% filter(is.na(sequence)))==0)

#  stopifnot(all(copies$sequence %in% rownames(se)))

#  # because the wgs and  rnaseq isn't paired per biosample, join by strain
#  # because not all seqs have detected wgs coverage, first perform crossing
#  copy_lookup <- crossing(sample=colnames(se),sequence=rownames(se)) %>%
#    left_join(copies,by=c("Strain","sequence")) %>%
#    dplyr::select(sample,sequence, est.copies)

#  # when there's no wgs data for a strain, fill in median of other strains
#  # furthermore, when estimated copies is less than 1 we assume that there is
#  # 1 degenerate copy or mismapping, either way doesn't make sense to multiply
#  # expression value by a bajillion when there is a tiny fractional estimated copy
#  copy_lookup <- copy_lookup %>%
#    group_by(sequence) %>%
#    mutate(median.copies = median(est.copies, na.rm = T)) %>%
#    ungroup() %>%
#    mutate(est.copies = ifelse(is.na(sequence),median.copies, est.copies)) %>%
#    mutate(est.copies = ifelse(est.copies < 1, 1, est.copies)) %>%
#    dplyr::select(-median.copies)

#  # widen and make a numeric matrix
#  copy_mat <- pivot_wider(copy_lookup, names_from = "sample",values_from = "est.copies") %>%
#    column_to_rownames("sequence") %>%
#    as.matrix()

#  # any nas - tes without counts or host genes - get replaces with 1
#  copy_mat[is.na(copy_mat)] <- 1

#  stopifnot(all(dim(copy_mat) == dim(se)))

#  ### need to reorder the mat first!!!!!
#  copy_mat <- copy_mat[rownames(se), colnames(se)]

#  stopifnot(all(colnames(se) == colnames(copy_mat)) & all(rownames(se) == rownames(copy_mat)))

#  sanity <- se@assays@data$counts
#  se@assays@data$counts <- se@assays@data$counts/copy_mat
#  se@assays@data$abundance <- se@assays@data$abundance/copy_mat
#}

# normalization options accessible via DESeq2
if (ot %in% c("normcts","vst","rlog","fpkm","cpm")) {
  library(DESeq2)
  dds <- DESeqDataSet(se, ~ 1)

  if (ot == "normcts") {
    dat <- counts(dds,normalized=T)
  } else if (ot == "vst") {
    dat <- assay(vst(dds,blind=T,fitType = snakemake@params[["DESEQ2_FITTYPE"]]))
  } else if (ot == "fpkm") {
    dat <- fpkm(dds)
  } else if (ot == "cpm") {
    dat <- fpm(dds)
  }

}

# other normalization options
if (ot %in% c("counts","abundance","ppc_qn","edaseq_uq","edaseq_qn")) {
  message(paste("extracting data"))
  dat <- assay(se,ifelse(ot=="abundance","abundance","counts"))

  if (ot == "ppc_qn") {
    message(paste("performing",ot))
    dat <- preprocessCore::normalize.quantiles(dat,copy = T)
    dimnames(dat) <- dimnames(dat)
  } else if (ot == "edaseq_uq") {
    message(paste("performing"),ot)
    dat <- EDASeq::betweenLaneNormalization(dat, which="upper", round = F, offset=F)
  } else if (ot == "edaseq_qn") {
    message(paste("performing"),ot)
    dat <- EDASeq::betweenLaneNormalization(dat, which="full", round = F, offset=F)
  }
}

feature <- rownames(dat)

dat <- cbind(as.data.frame(feature),dat)

#saveRDS(dds,snakemake@output[["dds"]])

o_gz <- gzfile(snakemake@output[["txt"]],'w')
write.table(dat,file = o_gz, quote=F, sep = "\t",row.names = F, col.names = T)
close(o_gz)
