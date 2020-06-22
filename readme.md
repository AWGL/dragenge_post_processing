# Dragen Post Processing

## Introduction

Perform post processing of the raw variant calls from the Dragen server. Works with the DragenGE pipeline.

Performs the following processes:

- Calculates Contamination
- Calculates Depth Statistics
- Calculates Sensitivity
- Filters Variant Calls
- Prepares VCF for Database
- Annotates VCF with VEP

## Install

```
git clone https://github.com/AWGL/dragenge_post_processing.git
cd dragenge_post_processing
conda env create -f env/dragenge_post_processing.yaml

source activate dragenge_post_processing

wget https://console.cloud.google.com/storage/browser/_details/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
gatk3-register GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

```

Configuration files are located in the config directory. Locations of files such as reference genome will need to be edited to get pipeline to work.

## Run

The pipeline uses nextflow as a workflow manager to make it easy to deploy the pipeline in local and cluster environments.

To run on cluster:

```
nextflow -C config/IlluminaTruSightOne/IlluminaTruSightOne_pbs.config \
run \
-E \
dragenge_post_processing.nf \
--bams /data/results/dragen_results/191010_D00501_0366_BH5JWHBCX3/IlluminaTruSightOne/\*/\*\{.bam,.bam.bai\} \
--vcf /data/results/dragen_results/191010_D00501_0366_BH5JWHBCX3/IlluminaTruSightOne/191010_D00501_0366_BH5JWHBCX3\{.vcf.gz,.vcf.gz.tbi\} \
--variables /data/results/dragen_results/191010_D00501_0366_BH5JWHBCX3/IlluminaTruSightOne/\*/\*.variables \
--publish_dir /share/data/results/dragen_temp/191010_D00501_0366_BH5JWHBCX3/IlluminaTruSightOne/results \
--sequencing_run 191010_D00501_0366_BH5JWHBCX3 \
-work-dir /share/data/results/dragen_temp/191010_D00501_0366_BH5JWHBCX3/IlluminaTruSightOne/work

```
