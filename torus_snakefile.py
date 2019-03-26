-*- eval: (snakemake-mode); -*-
# Snakemake pipeline for running torus multivariate analysis
#
# assumes some files already present form running enrichment_analysis.snakefile.py
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
#    nohup snakemake -s torus_snakefile.py --latency-wait 10 --keep-going --jobs 96 --cluster "sbatch --output={params.log}_%A.out --error={params.log}_%A.err --cpus-per-task=1 --ntasks=1 --mem-per-cpu={params.mem} --account=pi-mstephens --partition=broadwl --time=30:00:00 --job-name={params.jobname}" & 
#



from snakemake.utils import report
import pandas as pd 


dsc_annot1 = ["E8_TCM_T_48h", "E8_TCM_D_48h", 
              "NR2F2", "PGR_Demayo", "DNaseI", 
              "H3K27me3",
              "FAIRE", "PGR-A_Bagchi", "PGR-B_Bagchi", 
              "FOSL2", "FOXO1_DeMayo", "POLII_DeMayo", 
              "1_GATA2::1-Interval-Track", "PLZF",
              "hic_all_interacting_DT1_dTL4_D_48h",
              "eQTL_0.05_FDR"]
marks = ["h3k27ac",  "h3k4me1",  "h3k4me3"]
txs = ["C", "D"]
mark_annots = expand("{mark}-final-{tx}-48h", mark=marks, tx = txs)
dsc_annot1.extend(mark_annots)

base_dir = "base2/"
#Data
#assoc="/project/mstephens/data/external_public_supp/gwas_summary_statistics/raw_summary_statistics/23_and_me_ptb/23andMe_ga.tsv.gz"
#zscores = "23andMe_ga.zscores.tsv.gz"
zscores_ldsc = "../ga.zscores.tsv.gz"


rule all:
    input: expand("torus_" + base_dir + "ga.{annot}.oneplus.rst", annot = dsc_annot1), 
           "torus_" + base_dir + "ga.full.rst", "torus_" + base_dir + "ga.base.rst"
    #input: expand("annotations/{annot}_annotation.oneplus.ldsc.tsv.gz", annot = dsc_annot1)

#Z score files already made
#####

#Annotation files are in one_plus_annotations. Made with select_one_plus_annot.R

#Annotation files already made

#Multi
rule qtl_ldscannot:
    input:
        zscores = zscores_ldsc,
        annot='../../../one_plus_annot_torus/' + base_dir + '{annot}.oneplus.annot.tsv.gz'
    output:
        qtl = 'torus_' + base_dir + 'ga.{annot}.oneplus.rst'
    params:
        log="qtl", jobname="qtl", mem="15G"
    shell:
        "~/dap/torus_src/torus -d {input.zscores} -annot {input.annot} --load_zval -est -qtl > {output.qtl}"

rule qtl_full:
    input: zscores = zscores_ldsc, annot='../../../one_plus_annot_torus/' + base_dir + 'full.annot.tsv.gz'
    output: qtl = 'torus_' + base_dir + 'ga.full.rst'
    params: log="qtl", jobname="qtl", mem="15G"
    shell: "~/dap/torus_src/torus -d {input.zscores} -annot {input.annot} --load_zval -est -qtl -dump_prior full > {output.qtl}"

rule qtl_base:
    input: zscores = zscores_ldsc, annot='../../../one_plus_annot_torus/' + base_dir + 'base.annot.tsv.gz'
    output: qtl = 'torus_' + base_dir + 'ga.base.rst'
    params: log="qtl", jobname="qtl", mem="15G"
    shell: "~/dap/torus_src/torus -d {input.zscores} -annot {input.annot} --load_zval -est -qtl -dump_prior base > {output.qtl}"

