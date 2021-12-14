# MAG Snakemake Workflow


## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Running pipeline](#running-pipeline)
- [License](./LICENSE)
- [Citation](#citation)

# Overview

This pipeline has the scripts and modules used to generate, quality assess, and taxonomically classify skin MAGs which the pipeline will generate using different per sample and co-assembly approaches. The prok.Snakefile and euk.Snakefile detail the prokaryotic and eukaryotic analyses respectively. Note that for the eukaryotic analyses, EukRep and EukCC need to already have been run. To use the pipeline a file called coassembly_runs.txt and runs.txt must be produced which detail the co-assembly approaches and the per sample approaches used in our pipeline. The metagenomes must be put in the directory data/singlerun and data/coassembly for the single run and co-assembly approaches respectively. A sample coassembly_runs.txt and runs.txt file has been provided. This code has been adapted from https://github.com/Finn-Lab/MAG_Snakemake_wf/ accompanying the paper (https://www.nature.com/articles/s41596-021-00508-2). The viral analysis was done via the VIRify pipeline (https://github.com/EBI-Metagenomics/emg-viral-pipeline). 


# System Requirements

## Hardware Requirements
HPC with at least 500 gigabytes of memory

## Software Requirements
- parallel-fastq-dump v0.6.6
- KneadData v0.7.4 (Bowtie2  v2.3.5.1, Trimmommatic v0.39)
- SPAdes v.3.13.0
- metaWRAP v1.2.1
- INFERNAL v1.1.2
- tRNAScan-SE v2.0
- BBMap v.37.62
- GUNC v1.0.1
- CAT v5.2.1
- CheckM  v1.1.2
- Mash v.2.0
- MUMmer v3.23
- QUASTv5.0.2
- dRep v2.3.2 
- GTDB-Tkv1.0.2
- ncbi-genome-downloadv0.2.12
- EukRep v.0.6.7
- EukCC v0.2


# Running pipeline 

### Submitting jobs

To run pipeline on the runs specified in runs.txt and coassembly_runs.txt, submit jobs with SLURM scheduler:
```
snakemake -s prok.Snakefile --use-singularity --restart-times 3 -k -j 50 --cluster-config clusterconfig.yaml --cluster "sbatch -n {cluster.nCPU} --mem {cluster.mem} -e {cluster.error} -o {cluster.output} -t {cluster.time}"
```

Submit jobs with LSF scheduler:
```
snakemake -s prok.Snakefile --use-singularity --restart-times 3 -k --jobs 50 --cluster-config clusterconfig.yaml --cluster "bsub -n {cluster.nCPU} -M {cluster.mem} -e {cluster.error} -o {cluster.output}"
```

