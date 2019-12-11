#!/bin/bash
set -euo pipefail

#script to compare data with NA12878 goldstandard
#args: <test-sample-name> <test-sample-vcf> <test-sample-bed>


#compare test NA12878 with gold standard
/share/apps/rtg-distros/rtg-tools-3.7.1/rtg vcfeval \
-b /data/db/human/giab/3.3.2/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
--bed-regions "$3" \
-c "$2" \
-e /data/db/human/giab/3.3.2/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed \
-o "RTG" \
-t /data/db/human/gatk/2.8/b37/human_g1k_v37-SDF \
--sample "HG001,$1"
