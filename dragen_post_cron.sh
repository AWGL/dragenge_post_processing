#!/bin/bash

set -euo pipefail

dragen_results_dir=/share/data/results/dragen_results/
dragen_temp_dir=/share/data/results/dragen_temp/

# loop through each folder and find runs which have finished the dragen side of stuff

for path in $(find "$dragen_results_dir" -maxdepth 3 -mindepth 3 -type f -name "dragen_complete.txt" -exec dirname '{}' \;); do

  # for each of those runs find the post processing pipeline we need
  echo $path

  if [ ! -f "$path"/post_processing_started.txt ]; then

  touch "$path"/post_processing_started.txt
 


  runid=$(basename $(dirname "$path"))
  panel=$(basename "$path")
  echo $runid
  echo $panel

  # load the variables file which says which pipeline we need

  . "$path"/"$panel".variables

  echo $post_processing_pipeline
  echo $post_processing_pipeline_version
  
  source /home/transfer/miniconda3/bin/activate $post_processing_pipeline

  echo nextflow -C /data/diagnostics/pipelines/"$post_processing_pipeline"/"$post_processing_pipeline"-"$post_processing_pipeline_version"/"$panel"/"$panel"_pbs.config run -E /data/diagnostics/pipelines/"$post_processing_pipeline"/"$post_processing_pipeline"-"$post_processing_pipeline_version"/"$post_processing_pipeline".nf \
  --bams '"$dragen_results_dir"/"$runid"/"$panel"/*/*{.bam,.bam.bai}' \
  --vcf '"$dragen_results_dir"/"$runid"/"$panel"/"$runid"{.vcf.gz,.vcf.gz.tbi}' \
  --variables '"$dragen_results_dir"/"$runid"/"$panel"/*/*.variables' \
  --alignment_metrics '"$dragen_results_dir"/"$runid"/"$panel"/*/*.mapping_metrics.csv' \
  --publish_dir "$dragen_temp_dir"/"$runid"/"$panel"/results \
  --sequencing_run "$runid" \
  -work-dir "$dragen_temp_dir"/"$runid"/"$panel"/work 

  fi

done






