#!/usr/bin/env nextflow

/*
========================================================================================
Nextflow pipeline for germline variant calling
========================================================================================

Author: Joseph Halstead

*/


/*
========================================================================================
Set up Parameters
========================================================================================
*/

params.sequencing_run = '191010_D00501_0366_BH5JWHBCX3'
params.bams = '/media/joseph/Storage/test_data/IlluminaTruSightOne/*/*{.bam,.bam.bai}'
params.vcf =  "/media/joseph/Storage/test_data/IlluminaTruSightOne/${params.sequencing_run}{.vcf.gz,.vcf.gz.tbi}"
params.variables = '/media/joseph/Storage/test_data/IlluminaTruSightOne/*/*.variables'
params.alignment_metrics =  '/media/joseph/Storage/test_data/IlluminaTruSightOne/*/*.mapping_metrics.csv'
params.reference_genome =  '/media/joseph/Storage/genomic_resources/bwa/human_g1k_v37.fasta'
params.sequence_dict = '/media/joseph/Storage/genomic_resources/bwa/human_g1k_v37.dict'
params.capture_bed = 'config/IlluminaTruSightOne_ROI_b37.bed'
params.coverage_bed = 'config/IlluminaTruSightOne_ROI_b37_coverage.bed'
params.coverage_groups = 'config/IlluminaTruSightOne_ROI_b37_coverage_groups.txt'
params.gnomad_exomes = '/media/joseph/Storage/genomic_resources/gnomad/gnomad.exomes.r2.0.1.sites.vcf.gz'
params.gnomad_genomes = '/media/joseph/Storage/genomic_resources/gnomad/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz'
params.vep_threads = 1
params.vep_cache = '/media/joseph/Storage/genomic_resources/vep_caches/vep'
params.coverage_calculator_path = '/home/joseph/Documents/dragen/coverage/CoverageCalculatorPy/CoverageCalculatorPy.py'
params.min_coverage = 20
params.pipeline_version = '1'
params.high_confidence_snps = '/media/joseph/Storage/genomic_resources/1000g/1000G_phase1.snps.high_confidence.b37.vcf'
params.publish_dir = 'results/'

params.snps_min_qual = 10.0
params.snps_max_fs = 60.0
params.snps_min_mqranksum = -15.0
params.snps_min_readposranksum = -8.0

params.indels_min_qual = 10.0
params.indels_max_fs = 200.0
params.indels_max_sor = 10.0
params.indels_min_readposranksum = -20.0

params.min_genotype_depth = 8


/*
========================================================================================
Define initial files
========================================================================================
*/

reference_genome = file(params.reference_genome)
capture_bed = file(params.capture_bed)
sequence_dict = file(params.sequence_dict)
gnomad_exomes = file(params.gnomad_exomes)
gnomad_genomes = file(params.gnomad_genomes)
vep_cache = file(params.vep_cache)
coverage_bed = file(params.coverage_bed)
coverage_groups = file(params.coverage_groups)
high_confidence_snps = file(params.high_confidence_snps)

/*
========================================================================================
Define initial channels
========================================================================================
*/

Channel
  .fromPath(params.variables)
  .set{ variables_channel }

Channel
  .fromFilePairs(params.vcf, flat: true) 
  .set { raw_vcf }

Channel
  .fromFilePairs(params.bams, flat: true) 
  .set { original_bams }

Channel
  .fromPath(params.alignment_metrics)
  .set{ alignment_metrics }


original_bams.into{
    coverage_bams
    contamination_bams
}


variables_channel.into{
    variables_meta
    variables_ped
}

/*
========================================================================================
Main pipeline
========================================================================================
*/


// Use bcftools to restrict vcf to region of interest and remove chr from vcf chromosomes
process restrict_vcf_to_roi{

	input:
    set val(id), file(vcf), file(vcf_index) from raw_vcf 

    output:
    set file("${params.sequencing_run}.roi.vcf.gz"), file("${params.sequencing_run}.roi.vcf.gz.tbi") into roi_vcf_channel

	"""
	bcftools view -R $capture_bed $vcf > ${params.sequencing_run}.roi.vcf

	bgzip ${params.sequencing_run}.roi.vcf
	tabix ${params.sequencing_run}.roi.vcf.gz
	"""
}


// split channel into two so we can process snps and indels seperately
roi_vcf_channel.into {
  roi_vcf_for_snp_filtering
  roi_vcf_for_indel_filtering
}

// filter snps
process select_and_filter_snps{

    input:
    set file(vcf), file(vcf_index) from roi_vcf_for_snp_filtering 

    output:
    set file("${params.sequencing_run}.roi.snps.filtered.vcf.gz"), file("${params.sequencing_run}.roi.snps.filtered.vcf.gz.tbi") into filtered_snps_channel 

    """
    gatk SelectVariants \
    -V $vcf \
    -O ${params.sequencing_run}.roi.snps.vcf \
    --select-type-to-include SNP 

    gatk VariantFiltration \
    -V ${params.sequencing_run}.roi.snps.vcf \
    -O ${params.sequencing_run}.roi.snps.filtered.vcf \
    --filter-expression "QUAL < $params.snps_min_qual" \
    --filter-name "SNPLowQual" \
    --filter-expression "FS > $params.snps_max_fs" \
    --filter-name "SNPHighFS" \
    --filter-expression "MQRankSum < $params.snps_min_mqranksum" \
    --filter-name "SNPLowMQRankSum" \
    --filter-expression "ReadPosRankSum < $params.snps_min_readposranksum" \
    --filter-name "SNPLowReadPosRankSum"

    bgzip ${params.sequencing_run}.roi.snps.filtered.vcf
    tabix ${params.sequencing_run}.roi.snps.filtered.vcf.gz

    """

}

// filter non snps
process select_and_filter_non_snps{

    input:
    set file(vcf), file(vcf_index) from roi_vcf_for_indel_filtering 

    output:
    set file("${params.sequencing_run}.roi.indels.filtered.vcf.gz"), file("${params.sequencing_run}.roi.indels.filtered.vcf.gz.tbi") into filtered_indels_channel 

    """
    gatk SelectVariants \
    -V $vcf \
    -O ${params.sequencing_run}.roi.indels.vcf \
    --select-type-to-exclude SNP 

    gatk VariantFiltration \
    -V ${params.sequencing_run}.roi.indels.vcf \
    -O ${params.sequencing_run}.roi.indels.filtered.vcf \
    --filter-expression "QUAL < $params.snps_min_qual" \
    --filter-name "IndelLowQual" \
    --filter-expression "FS > $params.indels_max_fs" \
    --filter-name "IndelHighFS" \
    --filter-expression "SOR > $params.indels_max_sor" \
    --filter-name "IndelHighSOR" \
    --filter-expression "ReadPosRankSum < $params.indels_min_readposranksum" \
    --filter-name "IndelReadPosRankSum" 

    bgzip ${params.sequencing_run}.roi.indels.filtered.vcf
    tabix ${params.sequencing_run}.roi.indels.filtered.vcf.gz

    """

}

// merge filtered variants
process merge_indels_and_snps{

    input:
    set file(snp_vcf), file(snp_vcf_index) from filtered_snps_channel
    set file(indel_vcf), file(indel_vcf_index) from filtered_indels_channel

    output:
    set file("${params.sequencing_run}.roi.filtered.vcf.gz"), file("${params.sequencing_run}.roi.filtered.vcf.gz.tbi") into filtered_vcf_channel

    """
    gatk MergeVcfs \
    -I $snp_vcf \
    -I $indel_vcf \
    -O ${params.sequencing_run}.roi.filtered.vcf

    bgzip ${params.sequencing_run}.roi.filtered.vcf
    tabix ${params.sequencing_run}.roi.filtered.vcf.gz

    """
}

// mark genotypes with low depth with a filter
process filter_genotypes_with_low_depth{

    input:
    set file(vcf), file(vcf_index) from filtered_vcf_channel

    output:
    set file("${params.sequencing_run}.roi.filtered.gts.vcf.gz"), file("${params.sequencing_run}.roi.filtered.gts.vcf.gz.tbi") into filtered_gts_vcf_channel

    """
    gatk VariantFiltration \
    -V $vcf \
    -O ${params.sequencing_run}.roi.filtered.gts.vcf \
    --genotype-filter-expression "DP < $params.min_genotype_depth" \
    --genotype-filter-name "LowDP"

    bgzip ${params.sequencing_run}.roi.filtered.gts.vcf
    tabix ${params.sequencing_run}.roi.filtered.gts.vcf.gz
    """

}

// split filtered vcf channel into four
filtered_gts_vcf_channel.into{
    vcf_annotation_channel
    vcf_database_channel
    vcf_relatedness_channel
    vcf_contamination_channel
}


// split multiallelics and normalise
process split_multiallelics_and_normalise{

	input:
	set file(vcf), file(vcf_index) from vcf_annotation_channel 

	output:
	set file("${params.sequencing_run}.hard-filtered.roi.norm.vcf.gz"), file("${params.sequencing_run}.hard-filtered.roi.norm.vcf.gz.tbi") into normalised_vcf_channel

	"""
	zcat $vcf | vt decompose -s - | vt normalize -r $reference_genome - > ${params.sequencing_run}.hard-filtered.roi.norm.vcf
	bgzip ${params.sequencing_run}.hard-filtered.roi.norm.vcf
	tabix ${params.sequencing_run}.hard-filtered.roi.norm.vcf.gz

	"""

}

// Annotate using VEP
process annotate_with_vep{

    publishDir "${params.publish_dir}/annotated_vcf/", mode: 'copy'

    input:
    file(normalised_vcf) from normalised_vcf_channel

    output:
    set file("${params.sequencing_run}.hard-filtered.roi.norm.anno.vcf.gz"), file("${params.sequencing_run}.hard-filtered.roi.norm.anno.vcf.gz.tbi")  into annotated_vcf

    """
    vep \
    --verbose \
    --format vcf \
    --everything \
    --fork $params.vep_threads \
    --species homo_sapiens \
    --assembly GRCh37 \
    --input_file $normalised_vcf \
    --output_file ${params.sequencing_run}.hard-filtered.roi.norm.anno.vcf \
    --force_overwrite \
    --cache \
    --dir  $vep_cache \
    --fasta $reference_genome \
    --offline \
    --cache_version 94 \
    --no_escape \
    --shift_hgvs 1 \
    --vcf \
    --refseq \
    --flag_pick \
    --pick_order biotype,canonical,appris,tsl,ccds,rank,length \
    --exclude_predicted \
    --custom ${gnomad_genomes},gnomADg,vcf,exact,0,AF_POPMAX \
    --custom ${gnomad_exomes},gnomADe,vcf,exact,0,AF_POPMAX 

    bgzip ${params.sequencing_run}.hard-filtered.roi.norm.anno.vcf
    tabix ${params.sequencing_run}.hard-filtered.roi.norm.anno.vcf.gz

    """
}

// Calculate relatedness between samples
process calculate_relatedness {

    publishDir "${params.publish_dir}/relatedness/", mode: 'copy'

    input:
    set file(vcf), file(vcf_index) from vcf_relatedness_channel

    output:
    file("${params.sequencing_run}.relatedness2") into relatedness_channel

    """
    vcftools --relatedness2 \
    --out $params.sequencing_run \
    --gzvcf $vcf

    """

}

// calculate intersample contamination
process calculate_contamination {

    publishDir "${params.publish_dir}/contamination/", mode: 'copy'

    input:
    set val(id), file(bam), file(bam_index) from contamination_bams 

    output:
    file("${id}_contamination.selfSM") into contamination_channel

    """
    gatk SelectVariants \
    -R $reference_genome \
    -V $high_confidence_snps \
    -O high_confidence_snps.vcf \
    --select-type-to-include SNP \
    --restrict-alleles-to BIALLELIC \
    -L $capture_bed \
    --exclude-intervals X \
    --exclude-intervals Y \
    --exclude-intervals MT \
    --exclude-non-variants \
    --exclude-filtered \

    verifyBamID \
    --vcf high_confidence_snps.vcf \
    --bam $bam \
    --out ${id}_contamination \
    --verbose \
    --ignoreRG \
    --chip-none \
    --minMapQ 20 \
    --maxDepth 1000 \
    --precise
    """
}

// use the alignment metrics file to calculate the sex
process calculate_sex{

    publishDir "${params.publish_dir}/calculated_sex/", mode: 'copy'

    input:
    file(metrics_file) from alignment_metrics

    output:

    file("${metrics_file.simpleName}_calculated_sex.txt") into sex_channel

    """
    calculate_sex.py --file $metrics_file > ${metrics_file.simpleName}_calculated_sex.txt

    """

}

// create ped file
process create_ped {

    publishDir "${params.publish_dir}/ped/", mode: 'copy'

    input:
    file(variables) from variables_ped.collect()

    output:
    file("${params.sequencing_run}.ped") into ped_channel

    """
    create_ped.py --variables "*.variables" > ${params.sequencing_run}.ped

    """

}


// Variant database needs extra metadata from variables files
process collect_metadata_for_vcf{

    input:
    file(variables) from variables_meta.collect()

    output:
    file("${params.sequencing_run}_meta.txt") into meta_data_channel

    """
    for i in ${variables.collect { "$it" }.join(" ")}; do

        # source variables
        . \$i
        
        echo \\#\\#SAMPLE\\=\\<ID\\=\$sampleId,Tissue\\=Germline,WorklistId\\="\$worklistId",SeqId\\="\$seqId",Assay\\="\$panel",PipelineName\\=DragenGE,PipelineVersion\\="$params.pipeline_version"\\,RawSequenceQuality=0,PercentMapped=0,ATDropout=0,GCDropout=0,MeanInsertSize=0,SDInsertSize=0,DuplicationRate=0,TotalReads=0,PctSelectedBases=0,MeanOnTargetCoverage=0,PctTargetBasesCt=0,EstimatedContamination=0,GenotypicGender=None,TotalTargetedUsableBases=0,RemoteVcfFilePath=None,RemoteBamFilePath=None\\> >> ${params.sequencing_run}_meta.txt
    done

    """

}

// prepare vcf for upload to existing database
process prepare_vcf_for_database {

    publishDir "${params.publish_dir}/database_vcf/", mode: 'copy'

    input:
    set file(vcf), file(vcf_index) from vcf_database_channel 
    file(metadata) from meta_data_channel

    output:
    file("${params.sequencing_run}.hard-filtered.roi.database.vcf") into database_vcf

    """
    # get header 

    zcat $vcf | grep "^##" > header.txt

    cat header.txt $metadata > ${params.sequencing_run}.hard-filtered.roi.header.vcf

    zcat $vcf | grep -v "^##" >> ${params.sequencing_run}.hard-filtered.roi.header.vcf

    # remove mt variants

    awk '\$1 !~ /^MT/ { print \$0 }' ${params.sequencing_run}.hard-filtered.roi.header.vcf > ${params.sequencing_run}.hard-filtered.roi.database.vcf

    """


}


// Use sambamba to generate the per base coverage
process generate_coverage_file{

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

	input:
    set val(id), file(bam), file(bam_index) from coverage_bams 

    output:
    set val(id), file("${id}.per_base_coverage.gz"), file("${id}.per_base_coverage.gz.tbi") into per_base_coverage_channel

	"""
	sambamba depth base \
	-L $capture_bed \
	--min-coverage=0 \
	--min-base-quality=10 \
	--filter='mapping_quality > 19 and not duplicate and not failed_quality_control' \
	$bam | \
	grep -v REF | \
	awk '{OFS="\t"} {print \$1,\$2,\$3,\$3,\$3 }' | \
	sort -k1,1 -k2,2n | \
	bgzip > ${id}.per_base_coverage.gz

	tabix -b 2 -e 2 -s 1  ${id}.per_base_coverage.gz

	"""
}


// Use coverage calculator to get gaps etc
process generate_gaps_files{

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

	input:
    set val(id), file(depth_file), file(depth_file_index) from per_base_coverage_channel 

    output:
    set file("${id}.coverage"), file("${id}.gaps"),  file("${id}.missing"), file("${id}.totalCoverage") into gaps_channel

	"""
	python $params.coverage_calculator_path \
	-D $depth_file \
	-B $coverage_bed \
	-d $params.min_coverage \
	-o ${id} \
	-O ./ \
	-g $coverage_groups

	"""
}
