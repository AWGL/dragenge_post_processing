#!/bin/bash
set -euo pipefail

bam=${1}
sample_id=${2}
call_vcf=${3}
giab_baseline=${4}
reference_faidx=${5}
giab_high_confidence_bed=${6}
reference_genome_sdf=${7}
roi_bed=${8}
min_mapping_quality_coverage=${9} 
min_base_quality_coverage=${10}
reference_genome=${11}


# recalcualte depth of coverage using roi bed
gatk3 -T DepthOfCoverage \
-R $reference_genome \
-o "$sample_id"_depth_of_coverage \
-I $bam \
-L $roi_bed \
--countType COUNT_FRAGMENTS \
--minMappingQuality $min_mapping_quality_coverage \
--minBaseQuality $min_base_quality_coverage \
-ct 20 \
--omitLocusTable \
-rf MappingQualityUnavailable \
-dt NONE

sed 's/:/\t/g' "$sample_id"_depth_of_coverage \
| grep -v 'Locus' | sort -k1,1 -k2,2n | bgzip > "$sample_id"_depth_of_coverage.csv.gz

tabix -b 2 -e 2 -s 1 "$sample_id"_depth_of_coverage.csv.gz


# generate bed of regions with > 20x
zcat "$sample_id"_depth_of_coverage.csv.gz \
     | awk '$3>=20{print $1"\t"$2-1"\t"$2}' | grep -v "\-1" \
     | bedtools sort -faidx $reference_faidx \
     | bedtools merge > "$sample_id"_gt_eq_20x.bed

bedtools intersect -a "$sample_id"_gt_eq_20x.bed -b $roi_bed > "$sample_id"_final_roi.bed

rtg vcfeval \
-b $giab_baseline \
--bed-regions "$sample_id"_final_roi.bed   \
-c "$call_vcf" \
-e $giab_high_confidence_bed \
-o "RTG" \
-t $reference_genome_sdf \
--sample "HG001,$sample_id"


# R code to calculate sensitivity and 95% CI
Rscript -e 'df <- read.table("./RTG/summary.txt", skip=2); x <- binom.test(x = as.numeric(df[2,3]), n = (as.numeric(df[2,3]) + as.numeric(df[2,5]))); paste0("sensitivity: ",round(x$estimate, 3)," ",round(x$conf.int[1], 3),"-",round(x$conf.int[2], 3))'
