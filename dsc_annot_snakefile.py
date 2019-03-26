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
#    nohup snakemake -s dsc_annot_snakefile.py --latency-wait 10 --keep-going --jobs 96 --cluster "sbatch --output={params.log}_%A.out --error={params.log}_%A.err --cpus-per-task=1 --ntasks=1 --mem-per-cpu={params.mem} --account=pi-mstephens --partition=broadwl --time=1:00:00 --job-name={params.jobname}" & 
#



from snakemake.utils import report
import pandas as pd 


dsc_annotations = ["E8_TCM_T_48h", "E8_TCM_D_48h", 
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
dsc_annotations.extend(mark_annots)


#Data

rule all:
    input:  expand("dsc.{chrom}.annot.gz", chr=range(1, 23))

## Make z-scores file
rule var_bed:
    input:
        bim="../1000G_ldscores/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}.bim"
    output:
        var_bed =  "../1000G_ldscores/1000G_EUR_Plink_Variants/1000G.EUR.QC.{chrom}.variants.bed"
    params:
        log="var", jobname="var", mem="5G",
        temp = "../1000G_ldscores/1000G_EUR_Plink_Variants/1000G.EUR.QC.{chrom}.variants.bed.temp"
    shell:  
        """
        awk '{{print "chr"$1 "\t"$4"\t" ($4 + 1) "\t" $2}}' {input.bim} > {params.temp}
        ~/bedops/bin/sort-bed {params.temp} > {output.var_bed}
        rm {params.temp}
        """ 
        

rule bed_overlap:
    input:
        var_bed = "../1000G_ldscores/1000G_EUR_Plink_Variants/1000G.EUR.QC.{chrom}.variants.bed",
        dsc_bed = "../dsc_annotation/{annot}.bed"
    output:
        var_annot = "{annot}.{chrom}.annot"
    resources:
        log="procvar2",
        jobname="procvar2",
        mem="2G"
    shell: 
        "~/bedops/bin/bedops -e {input.var_bed} {input.dsc_bed} > {output.var_annot}"


rule annot:
    input:
        dsc_beds = expand("{annot}.{{chrom}}.annot", annot = dsc_annotations),
        bim="../1000G_ldscores/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}.bim"
    output:
        annot = "dsc.{chrom}.annot.gz"
    params:
        chrom="{chrom}"
    resources:
        log="annot2",
        jobname="annot2",
        mem="5G"
    script:
        "R/write_annot.R"
