A summary of what is in this directory

# Annotation files:

`dsc_annotation`: This directory contains bed files and annotation files (one line per SNP) for DSC annotations.
    dsc.*.annot.gz are annotation files for 1k genomes variants for only dsc annotations. One file per chromosome. Made using Snakemake file dsc_annot_snakefile.py
    ldscores: contains ldscores for each DSC annotation. One file per annotation per chromosome. 

`1000G_ldscores`: Downloaded ldscores for 1k genomes variants. From https://data.broadinstitute.org/alkesgroup/LDSCORE/

`baseline_bedfiles`: .bed files for baseline annotations. I don't think I have used these.

`eQTL`, `Hi-C` both contain original data supplied by Carole and Marcelo respectively. These were used to create the relevant bed files in `dsc_annotation`

`annot_name_dict.txt`: Text file giving a transaltion of ugly file name roots into readable annotation labels. 

`one_plus_annot_torus`: Torus (as far as I can tell) can't use a subset of annotations, so if you want to change the annotation set you run with you need to write a whole new annotation file. I experimented with a couple different "base" sets labeled as 
    `base1`
    `base2`
    `big_base` (everything)
I have only successfully run anything with base1. So these other directories could be deleted potentially. Files created with `base1.R`. 
    `base1`: Contains "one plus" annotation files for 1k genomes variants. Each of these files is the "base1" annotations (there are 38 of these) plus one of the DSC annotations. For example FOXO1.oneplus.annot.tsv.gz is base1 plus FOXO1. 
        full.annot.tsv.gz: base1 + all DSC annotations
        base.annot.tsv.gz: base1 annotations only
    full_annotations.torus.RDS: A data frame of all 98 DSC and non-DSC annotation for 1k genomes SNPs. 

`one_plus_annot_ldsc`: Contains files used for "one plus" annotation enrichment using stratified LDSC.


# Results files

`gestational_age`: Most results for Zhang et al gestational age data are here. Generally, results are produced by a snakemake file.
    `R`: R scripts live here
    `torus_base1`: Torus results for "base1" and one-plus annotation sets.  See torus_snakefile.py
        `base`: prior files created using only base1
        `full`: prior files created using base1 + DSC
        .`rst` files are output by torus.
        `.RDS` files assorted tables of results formated with R
    `zscores2`: Data for just the 18 regions we have been focusing on. See finemap_snakefile.py
        `data.*.RDS` Full GWAS data in these regions
        `data.filtered.*.RDS`: After filtering out ambiguous variants and flipping strands to match reference.
        `snps.filtered.*.txt` Just snplist
        `fullannot.filtered.*.tsv` Annotations for variants in each region
    `susie_finemap`: Susie results. See finemap_snakefile.py
    `dap_finemap`: DAP finemapping results finemap_snakefile.py
    `caviarbf_files`: Files and results for caviarbf. No Snakemake files for these. 

strat_ldscore: Old stratified ld score regression results for ptb, and two ukbb traits.

Where do GWAS results live?
    Zhang et al GWAS:
    `/project/mstephens/data/external_public_supp/gwas_summary_statistics/summary_statistics/23andme_ga_summary_statistics.tsv.gz`: Formatted results.
    `/project/mstephens/data/external_public_supp/gwas_summary_statistics/raw_summary_statistics/23_and_me_ptb/`: Everything except ga_gwas_181510 is raw data from Zhang et al. The directories `Zhang_2017*` are what I was given plus my processing scripts.

    new GWAS: 
    `/project/mstephens/data/external_public_supp/gwas_summary_statistics/raw_summary_statistics/23_and_me_ptb/ga_gwas_181510`:  Raw data plus assorted formatting attempts to get it to run through the pipeline. Formatted results coming soon. 


    
