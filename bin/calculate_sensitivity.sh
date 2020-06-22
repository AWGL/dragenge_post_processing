#!/bin/bash
set -euo pipefail

depth_of_coverage_file=$1
sample_id=$2
call_vcf=$3
giab_baseline=$4
reference_faidx=$5
giab_high_confidence_bed=$6
reference_genome_sdf=$7

# generate bed of regions with > 20x
zcat $depth_of_coverage_file \
     | awk '$3>=20{print $1"\t"$2-1"\t"$2}' | grep -v "\-1" \
     | bedtools sort -faidx $reference_faidx \
     | bedtools merge > "$sample_id"_gt_eq_20x.bed  


bgzip $call_vcf
tabix "$call_vcf".gz

rtg vcfeval \
-b $giab_baseline \
--bed-regions "$sample_id"_gt_eq_20x.bed   \
-c "$call_vcf".gz \
-e $giab_high_confidence_bed \
-o "RTG" \
-t $reference_genome_sdf \
--sample "HG001,$sample_id"


# R code to calculate sensitivity and 95% CI
Rscript -e 'df <- read.table("./RTG/summary.txt", skip=2); x <- binom.test(x = as.numeric(df[2,3]), n = (as.numeric(df[2,3]) + as.numeric(df[2,5]))); paste0("sensitivity: ",round(x$estimate, 3)," ",round(x$conf.int[1], 3),"-",round(x$conf.int[2], 3))'
