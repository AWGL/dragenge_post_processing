#!/bin/bash
set -euo pipefail

#Legacy script from the old GermlineEnrichment Pipline

ROI_BED=$1
HOTSPOTS=$2
DEPTH=$3
ID=$4
MINIMUM_COVERAGE=$5
REF_SEQ=$6


#Make BED file of all genes overlapping ROI
bedtools intersect -wa \
-a $REF_SEQ \
-b $ROI_BED | \
awk -F "\t" '$3 == "gene" { print $1"\t"$4-1"\t"$5 }' | \
sort -k1,1V -k2,2n -k3,3n | \
bedtools merge > "$ID"_TargetGenes.bed

#make target bed
bedtools intersect \
-a $REF_SEQ \
-b "$ID"_TargetGenes.bed | \
grep "NM_[0-9]*\.[0-9]*" | \
awk -F "\t" '$3 == "exon" { print $1"\t"$4-1"\t"$5 }' | \
sort -k1,1V -k2,2n -k3,3n | \
bedtools merge > "$ID"_Targets.bed

awk '$1 ~ /^MT/' $ROI_BED >> "$ID"_Targets.bed

#Intersect CDS for all genes, pad by p=n and merge coordinates by gene
bedtools intersect \
-a $REF_SEQ \
-b "$ID"_Targets.bed | \
awk -F'[\t|;|=]' -v p=5 '$3 == "CDS" { gene="null"; for (i=9;i<NF;i++) if ($i=="gene"){gene=$(i+1); break}; genes[gene] = genes[gene]$1"\t"($4-1)-p"\t"$5+p"\t"gene";" } END { for (gene in genes) print genes[gene] }' | \
while read line; do
    echo "$line" | \
    tr ';' '\n'| \
    sort -k1,1V -k2,2n -k3,3n | \
    bedtools merge -c 4 -o distinct;
done | \
sort -k1,1V -k2,2n -k3,3n > "$ID"_ClinicalCoverageTargets.bed

cat $HOTSPOTS "$ID"_ClinicalCoverageTargets.bed | \
sort -k1,1V -k2,2n -k3,3n > "$ID"_ClinicalCoverageTargetsHotspots.bed

#Make PASS BED
tabix -R "$ID"_ClinicalCoverageTargetsHotspots.bed \
$DEPTH | \
awk -v minimumCoverage="$MINIMUM_COVERAGE" '$3 >= minimumCoverage { print $1"\t"$2-1"\t"$2 }' | \
sort -k1,1V -k2,2n -k3,3n | \
bedtools merge > "$ID"_PASS.bed

#Calculate overlap between PASS BED and ClinicalCoverageTargets
bedtools coverage \
-a "$ID"_ClinicalCoverageTargetsHotspots.bed \
-b "$ID"_PASS.bed | \
tee "$ID"_clinical_coverage_target_metrics.txt | \
awk '{pass[$4]+=$6; len[$4]+=$7} END { for(i in pass) printf "%s\t %.2f%%\n", i, (pass[i]/len[i]) * 100 }' | \
sort -k1,1 > "$ID"_clinical_coverage_gene_coverage.txt

#Make GAP BED
bedtools subtract \
-a "$ID"_ClinicalCoverageTargetsHotspots.bed \
-b "$ID"_PASS.bed | \
sort -k1,1V -k2,2n -k3,3n \
> "$ID"_gaps.bed

