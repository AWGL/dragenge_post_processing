#!/bin/bash
set -euo pipefail

# Author: Christopher Medway
# Use: calculate_sensitivity.sh <seqId sampleId>
# i.e. bash /data/diagnostics/scripts/calculate_sensitivity.sh /data/results/180615_D00501_0201_AHK2TVBCX2/IlluminaTruSightOne 18M01315

seqId=$1
sampleId=$2

echo $seqId
echo $sampleId

# generate bed of regions with > 20x
zcat "$seqId"/"$sampleId/"*_DepthOfCoverage.gz \
     | awk '$3>=20{print $1"\t"$2-1"\t"$2}' \
     | /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools sort \
         -faidx /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai \
     | /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge \
         > "$seqId"/"$sampleId"/"$sampleId"_gt_eq_20x.bed  


bash rtg.sh \
    $sampleId \
    "$seqId"/*_filtered_annotated_padded.vcf.gz \
    "$seqId"/"$sampleId"/"$sampleId"_gt_eq_20x.bed

# R code to calculate sensitivity and 95% CI
Rscript -e 'df <- read.table("./RTG/summary.txt", skip=2); x <- binom.test(x = as.numeric(df[2,3]), n = (as.numeric(df[2,3]) + as.numeric(df[2,5]))); paste0("sensitivity: ",round(x$estimate, 3)," ",round(x$conf.int[1], 3),"-",round(x$conf.int[2], 3))'

rm -r ./RTG
