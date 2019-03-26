-*- eval: (snakemake-mode); -*-
# Snakemake pipeline for running fine mapping with DAP-G and SuSiE
#
#
# LICENSE: CC0. Do what you want with the code, but it has no guarantees.
#          https://creativecommons.org/share-your-work/public-domain/cc0/
#
# Installation with conda package manager:
#
# conda create -n ptb graphviz imagemagick pandas python=3.6 snakemake xlrd
# source activate ptb 
#

#from snakemake.utils import R

#To run the full pipeline, submit the following line from within the
#same directory as the Snakefile while on the head node (the paths to
#        the data files are relative to the Snakefile):
#    nohup snakemake -s finemap_snakefile.py --latency-wait 10 --keep-going --jobs 96 --cluster "sbatch --output={params.log}_%A.out --error={params.log}_%A.err --cpus-per-task={params.cores} --ntasks=1 --mem-per-cpu={params.mem} --account=pi-xinhe --partition=broadwl --time={params.time} --job-name={params.jobname}" & 
#



from snakemake.utils import report
import pandas as pd 

windows = ["118", "1224", "149", 
          "1574", "15", "1682", "181", 
          "353", "356", "363", "437", 
          "480", "512", "614", "656", 
          "704", "949", "986"]

priors = ["full", "base"]

plink_dir = "plink_output/"
zscore_dir = "zscores2/"
susie_dir = "susie_finemap/"
dap_dir = "dap_finemap/"


####Caviarbf params
caviarbf_dir = "caviarbf_files"
thrs = [0.001]
fit = ["lasso"]
eps = [0.01]
meth = ["cv", "topK"]
lam = [-1, 0] 

rule all:
    input:
        expand(dap_dir + "dapfinemap.{pr}.{b}.result" ,pr = priors, b = windows),


zscores = "ga.zscores.tsv.gz"
data = "/project/mstephens/data/external_public_supp/gwas_summary_statistics/raw_summary_statistics/23_and_me_ptb/23andMe_ga.RDS"
plink_root = "/project2/xinhe/1kg/plink_format/EUR.1kg"

windowstr = " ".join(windows)
rule format1:
    input:
        zscores = zscores,
        prior_base=expand("torus_base1/base/{b}.prior",b=windows)
        prior_full=expand("torus_base1/full/{b}.prior",b=windows)
        data = data
    output:
        data_f=expand(zscore_dir + "data.{b}.RDS", b = windows),
        snp_f=expand(zscore_dir+"snps.{b}.txt", b = windows)
    resources:
        log="f1",
        jobname="f1",
        mem="10G",
    script:
        "R/format1.R"

   
rule frq_plink:
    input:
        snps = zscore_dir + "snps.{b}.txt",
        bed = plink_root + ".bed",
        bim = plink_root + ".bim",
        fam = plink_root + ".fam"
    output:
        frq = plink_dir + "{b}.frq"
    params:
        log="frq",
        jobname="frq",
        mem="10G",
        plink_root = plink_root,
        dir = plink_dir
    shell:
        "plink --bfile {params.plink_root} --extract {input.snps} --freq --out {params.dir}{wildcards.b}"

rule process:
    input:
        data = zscore_dir + "data.{b}.RDS",
        frq = plink_dir + "{b}.frq"
    output:
        data = zscore_dir + "data.filtered.{b}.RDS",
        snps = zscore_dir + "snps.filtered.{b}.txt",
        zscores = zscore_dir + "zscores.filtered.{b}.tsv"
    resources:
        log="f2",
        jobname="f2",
        mem="1G"
    script: "R/strand.R"

rule ld_plink:
    input:
        snps = zscore_dir + "snps.filtered.{b}.txt",
        bed = plink_root + ".bed",
        bim = plink_root + ".bim",
        fam = plink_root + ".fam"
    output:
        frq = plink_dir + "{b}.filtered.frq",
        ld = plink_dir + "{b}.filtered.ld"
    params:
        pdir = plink_dir
    resources:
        log="ld",
        jobname="ld",
        mem="10G",
        plink_root = plink_root,

    shell: "plink --bfile {params.plink_root} --extract {input.snps} --freq \
            --r square --out {params.pdir}{wildcards.b}.filtered"

rule dap:
    input:
        zscores = zscore_dir + "zscores.filtered.{b}.tsv",
        ld = plink_dir + "{b}.filtered.ld",
        prior = "torus_base1/{pr}/{b}.prior"
    output:
        out = dap_dir + "dapfinemap.{pr}.{b}.result"
    resources:
        log="dap",
        jobname="dap",
        mem="40G"
    shell: "~/dap/dap_src/dap-g -d_z {input.zscores} \
                                 -d_ld {input.ld} \
                                 -o {output.out} \
                                 --output_all -p {input.prior}" 

rule susie:
    input: data = zscore_dir + "data.filtered.{b}.RDS", 
           ld = plink_dir + "{b}.filtered.ld"
    output: out = susie_dir + "susie.{b}.RDS"
    params: log="susie", jobname="susie", mem="10G", dir = susie_dir
    shell: "Rscript R/susie.R {wildcards.b} {input.data} {input.ld} {params.dir}"

rule caviarbf:
    input: locus_list = caviarbf_dir + "/trim{t}/ga_loci.txt"
    output: caviarbf_dir + "/thresh{t}/{f}_eps{e}_{m}_lam{l}.loglik"
    params: log = "caviarbf", jobname="cbf", mem="40G", cores=1, time="24:00:00"
    shell: "Rscript R/run_caviarbf.R {wildcards.f} {params.cores} {wildcards.e} {wildcards.m} caviarbf_files/trim{wildcards.t} {wildcards.l}" 

