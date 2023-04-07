rule se_export_txt:
    """
    Generic rule for exporting raw or normalized expression from a serialized
    SummarizedExperiment object.
    """
    input:
        se=quant_wf("results/quantification/{quant_pipeline}/se.{feature_level}.{quant_rep}.rds"),
    output:
        txt="results/quantification/{quant_pipeline}/{sex}.{feature_level}.{cnnorm}.{expression_unit}.{quant_rep}.tsv.gz"
    resources:
        time=240,
        mem=20000,
        cpus=1
    params:
        DESEQ2_FITTYPE = config.get("DESEQ2_FITTYPE"),
    script:
        "../scripts/se_export_txt.R"

rule filter_transform_scale:
    """
    This rule takes parameters from a given linear model run section in config.yaml
    and filters, transforms, and scales accordingly. Optionally, it will correct the centered and scaled matrix
    by regressing out a specified number of prinicipal components.
    """
    input:
        sefile = lambda wc: quant_wf("results/quantification/{p}/se.{fl}.{{quant_rep}}.rds".format(p=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_PIPELINE"),fl=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_FEATURE_LEVEL"))),
        mat =lambda wc: "results/quantification/{p}/{s}.{fl}.{c}.{u}.{{quant_rep}}.tsv.gz".format(p=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_PIPELINE"),s=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_SEX"),fl=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_FEATURE_LEVEL"), c=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_CNNORM"),u=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_UNITS")),
    output:
        mat = "results/linear_models/{model_id}/{quant_rep}/expression.tsv.gz"
    params:
        filt = config.get("LM_FEATURE_FILT"),
        transforms = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_VARIABLE_TRANSFORM"),
        scale = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_VARIABLE_SCALE"),
        pcs = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_CORRECT_N_PCS") if config.get("RUN_TYPE") != "test" else 0,
    resources:
        time=20,
        mem=24000,
        cpus=1
    script:
        "../scripts/filter_transform_scale.R"

rule gene_x_gene_corr:
    input:
        expr = rules.filter_transform_scale.output.mat
    output:
        rds = "results/linear_models/{model_id}/{quant_rep}/gene_x_gene.corrr.rds",
    conda:
        "../envs/corrr.yaml"
    resources:
        time=20,
        mem=48000,
        cpus=1
    script:
        "../scripts/gene_x_gene_corr.R"

rule mean_gene_corr:
    input:
        expand("results/linear_models/{{model_id}}/{r}/gene_x_gene.corrr.rds",r=REP_LIST)
    output:
        rds = "results/linear_models/{model_id}.expression.corrr.rds"
    conda:
        "../envs/corrr.yaml"
    params:
        test = False if config.get("RUN_TYPE","full") != "test" else True,
    resources:
        time=60,
        mem=96000,
        cpus=1
    script:
        "../scripts/mean_gene_corr.R"

rule gene_gene_edges:
    input:
        cor = rules.mean_gene_corr.output.rds,
    output:
        edges_rds = "results/linear_models/{model_id}.gene_2_gene_edges.rds",
        nudges_rds = "results/linear_models/{model_id}.gene_2_gene_nudges.rds",       
    conda:
        "../envs/corrr.yaml"
    params:
        model_id = "{model_id}"
    resources:
        time=60,
        mem=96000,
    script:
        "../scripts/gene_gene_edges.R"

rule scatter_genes_for_lm:
    """
    Note: Escape brackets on wcs when using in conjunction w/ scattergather.

    https://www.unix.com/shell-programming-and-scripting/247089-creating-all-possible-bi-combinations-list-grep-awk.html
    """
    input:
        expression = rules.filter_transform_scale.output.mat
    output:
        temp(expand("results/linear_models/{{model_id}}/{{quant_rep}}/chunk_{ch}",ch=[str(x).zfill(4) for x in range(0,config.get("LM_CHUNKS"))]))
    params:
        split_invoc = "split -e -d -a 4 -n r/{ch} - results/linear_models/{{model_id}}/{{quant_rep}}/chunk_".format(ch= config.get("LM_CHUNKS"),),
        subsample = config.get("RUN_TYPE") == "test",
    resources:
        time=240,
        mem=48000,
        cpus=1
    script:
        "../scripts/get_all_combos_tes_as_y.R"

rule chunked_linear_model:
    """
    Linear models were fit with R's `lm` function. See configfile and script for other details.
    """
    input:
        chunk = "results/linear_models/{model_id}/{quant_rep}/chunk_{lmchunk}",
        dat = rules.filter_transform_scale.output.mat,
        cd = quant_wf("results/meta/metadata.csv"),
        ol = tidal_wf("results/overlaps/overlaps.tsv.gz"),
        edges = rules.gene_gene_edges.output.edges_rds,
        copies = wgs_wf("results/copies/copies.tsv") if config.get("INCL_COPY_ESTIMATION_IN_EXPORT") else rules.dummy_copies.output.feats
    output:
        tidy = temp("results/linear_models/{model_id}/{quant_rep}/chunk_{lmchunk}.tidy.tsv"),
        glance = temp("results/linear_models/{model_id}/{quant_rep}/chunk_{lmchunk}.glance.tsv"),
        anova = temp("results/linear_models/{model_id}/{quant_rep}/chunk_{lmchunk}.anova.tsv"),
        breusch_pagan = temp("results/linear_models/{model_id}/{quant_rep}/chunk_{lmchunk}.breusch_pagan.tsv"),
        rainbow = temp("results/linear_models/{model_id}/{quant_rep}/chunk_{lmchunk}.rainbow.tsv"),
    params:
        formula = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FORMULA"),
        alt_formula = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_ALT_FORMULA"),
        mads_filter = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_MADS_FILTER"),
    conda:
        "../envs/baselm_v1.yaml"
    resources:
        time=20,
        mem=24000,
        cpus=1
    script:
        "../scripts/linear_model_baselm_v01.R"

rule collect_chunked_linear_models:
    """
    Note: Escape brackets on wcs when using in conjunction w/ scattergather.
    """
    input:
        lambda wc: expand("results/linear_models/{{model_id}}/{{quant_rep}}/chunk_{ch}.{lmr}.tsv",lmr=wc.lmresult,ch = [str(x).zfill(4) for x in range(0,config.get("LM_CHUNKS",80))])
    output:
        "results/linear_models/{model_id}/{quant_rep}/lm.{lmresult}.tsv.gz"
    resources:
        time=60,
        mem=24000,
        cpus=1
    shell:
        "xsv cat rows -d '\t' {input} | tr ',' '\t' | gzip -c > {output}"

rule correct_tls_lm_coefs:
    input:
        coefs = "results/linear_models/{model_id}/{quant_rep}/lm.tidy.tsv.gz",
        expression = lambda wc: "results/quantification/{p}/{s}.{fl}.{c}.{u}.{{quant_rep}}.tsv.gz".format(p=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_PIPELINE"),s=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_SEX"),fl=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_FEATURE_LEVEL"), c=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_CNNORM"),u=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_UNITS")),
    output:
        tsv = "results/linear_models/{model_id}/{quant_rep}/lm.tidy.corrected.tsv.gz"
    resources:
        time=60,
        mem=24000,
        cpus=1
    conda:
        "../envs/correct_coefs.yaml"
    script:
        "../scripts/correct_lm_coefs_baselm.R"


rule collect_lm_info_per_rep:
    input:
        anova= "results/linear_models/{model_id}/{quant_rep}/lm.anova.tsv.gz",
        rainbow= "results/linear_models/{model_id}/{quant_rep}/lm.rainbow.tsv.gz",
        breusch_pagan= "results/linear_models/{model_id}/{quant_rep}/lm.breusch_pagan.tsv.gz",
        glance= "results/linear_models/{model_id}/{quant_rep}/lm.glance.tsv.gz",
        tidy_correct= rules.correct_tls_lm_coefs.output.tsv
    output:
        tsv= "results/linear_models/{model_id}/{quant_rep}/lm.collected-info.tsv.gz"
    params:
        model = "{model_id}",
        rep =  "{quant_rep}"
    conda:
        "../envs/baselm_v1.yaml"
    script:
        "../scripts/collect_lm_info_per_rep.R"

rule collect_lm_info_per_model:
    input:
        info = expand("results/linear_models/{{model_id}}/{r}/lm.collected-info.tsv.gz",r=REP_LIST)
    output:
        tsv= "results/linear_models/{model_id}.collected-info.tsv.gz"
    conda:
        "../envs/baselm_v1.yaml"
    resources:
        time=60,
        mem=96000,
        cpus=1
    script:
        "../scripts/collect_lm_info_per_model.R"


rule export_linear_models:
    input:
        expand("results/linear_models/{m}.collected-info.tsv.gz",m=config.get("LMS_TO_EXPORT"))
    output:
        tsv = temp("results/linear_models/final-models.collected-info.tsv"),
        gz = "results/linear_models/final-models.collected-info.tsv.gz"
    shell:
        """
        zgrep 'anova' {input[0]} > {output.tsv} &&
        zcat {input} | grep -h -v 'anova' >> {output.tsv} &&
        gzip -c {output.tsv} > {output.gz}
        """