#!/bin/bash
snakemake --verbose -j 10 --use-singularity --singularity-args="--bind /project2/xinhe:/project2/xinhe,/scratch/midway2/nwknoblauch/ptb/" --singularity-prefix=/project2/xinhe/software/singularity --cluster-config rcc.yaml --cluster "sbatch  --nodes={cluster.nodes}  --time={cluster.time} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpuspertask} --partition={cluster.partition} --mem={cluster.mem}"