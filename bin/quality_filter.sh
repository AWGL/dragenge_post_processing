 #!/bin/bash
set -euo pipefail

# Add Quality filter flags to the VCF FILTER field

vcf=${1}
sequencing_run=${2}
snp_qual=${3}
indel_qual=${4} 
reference_genome=${5}
java_options=${6}
min_dp_variant_reporting=${7}


#Select SNPs
gatk3 $java_options \
-T SelectVariants \
-R $reference_genome \
-V $vcf \
-selectType SNP \
-o "$sequencing_run"_snps.vcf \
-dt NONE

#Filter SNPs
gatk3 $java_options \
-T VariantFiltration \
-R $reference_genome \
-V "$sequencing_run"_snps.vcf \
--filterExpression "QUAL < $snp_qual" \
--filterName "LowQual" \
-o "$sequencing_run"_snps_filtered.vcf \
-dt NONE

#Select non-snps (INDEL, MIXED, MNP, SYMBOLIC, NO_VARIATION)
gatk3 $java_options \
-T SelectVariants \
-R $reference_genome \
-V $vcf \
--selectTypeToExclude SNP \
-o "$sequencing_run"_non_snps.vcf \
-dt NONE

#Filter non-snps (INDEL, MIXED, MNP, SYMBOLIC, NO_VARIATION)
gatk3 $java_options \
-T VariantFiltration \
-R $reference_genome \
-V "$sequencing_run"_non_snps.vcf \
--filterExpression "QUAL < $indel_qual" \
--filterName "LowQual" \
-o "$sequencing_run"_non_snps_filtered.vcf \
-dt NONE

#Combine filtered VCF files
gatk3 $java_options \
-T CombineVariants \
-R $reference_genome \
--variant "$sequencing_run"_snps_filtered.vcf \
--variant "$sequencing_run"_non_snps_filtered.vcf \
-o "$sequencing_run"_roi_qual1.vcf \
-genotypeMergeOptions UNSORTED \
-dt NONE

#filter genotypes
gatk3 $java_options \
-T VariantFiltration \
-R $reference_genome \
-V "$sequencing_run"_roi_qual1.vcf \
--genotypeFilterExpression "DP < $min_dp_variant_reporting" \
--genotypeFilterName "LowDP" \
-o "$sequencing_run"_roi_qual.vcf \
-dt NONE


bgzip "$sequencing_run"_roi_qual.vcf
tabix "$sequencing_run"_roi_qual.vcf.gz