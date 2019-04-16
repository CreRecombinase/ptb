configfile:
    "ptb_config.yml"

include: "finemap_snakefile.py"
include: "torus_snakefile.py"
include: "dsc_annot_snakefile.py"
include: "zscores_snakefile"

rule all:
    input:
        config["ptb_scratch"]+"0_ga_gwas_z.txt.gz",
        config["ptb_scratch"]+"0_localga_gwas_z.txt.gz"
