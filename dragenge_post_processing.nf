#!/usr/bin/env nextflow

/*
========================================================================================
Nextflow pipeline for germline variant calling
========================================================================================

Author: Joseph Halstead

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
hotspots = file(params.hotspots)
refseq = file(params.refseq)

reference_genome_faidx = file(params.reference_genome_faidx)
giab_baseline = file(params.giab_baseline)
giab_high_confidence_bed = file(params.giab_high_confidence_bed)
giab_reference_genome_sdf = file(params.giab_reference_genome_sdf)

/*
========================================================================================
Define initial channels
========================================================================================
*/

Channel
  .fromPath(params.variables)
  .ifEmpty { exit 1, "Cannot find any variables files matching: ${params.variables}" }
  .set{ variables_channel }

Channel
  .fromFilePairs(params.vcf, flat: true)
  .ifEmpty { exit 1, "Cannot find a vcf file files matching: ${params.vcf}" }
  .set { raw_vcf }

Channel
  .fromFilePairs(params.bams, flat: true)
  .ifEmpty { exit 1, "Cannot find any bam files files matching: ${params.bams}" }
  .set { original_bams }


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


// Use bcftools to restrict vcf to region of interest
process restrict_vcf_to_roi{

    cpus params.vcf_processing_cpus

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


// filter snps
process select_and_filter_snps{

    cpus params.vcf_processing_cpus 

    input:
    set file(vcf), file(vcf_index) from roi_vcf_channel 

    output:
    set file("${params.sequencing_run}.roi.filtered.vcf.gz"), file("${params.sequencing_run}.roi.filtered.vcf.gz.tbi") into filtered_vcf_channel 

    """

    gatk VariantFiltration \
    -V $vcf \
    -O ${params.sequencing_run}.roi.filtered.vcf \
    --filter-expression "QUAL < $params.snps_min_qual" \
    --filter-name "LowQual" \
    --filter-expression "FS > $params.snps_max_fs" \
    --filter-name "HighFS" \
    --filter-expression "MQRankSum < $params.snps_min_mqranksum" \
    --filter-name "LowMQRankSum" \
    --filter-expression "ReadPosRankSum < $params.snps_min_readposranksum" \
    --filter-name "LowReadPosRankSum"

    bgzip ${params.sequencing_run}.roi.filtered.vcf
    tabix ${params.sequencing_run}.roi.filtered.vcf.gz

    """

}


// mark genotypes with low depth with a filter
process filter_genotypes_with_low_depth{

    cpus params.vcf_processing_cpus

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

    cpus params.vcf_processing_cpus

	input:
	set file(vcf), file(vcf_index) from vcf_annotation_channel 

	output:
	set file("${params.sequencing_run}.roi.filtered.norm.vcf.gz"), file("${params.sequencing_run}.roi.filtered.norm.vcf.gz.tbi") into normalised_vcf_channel

	"""
	zcat $vcf | vt decompose -s - | vt normalize -r $reference_genome - > ${params.sequencing_run}.roi.filtered.norm.vcf
	bgzip ${params.sequencing_run}.roi.filtered.norm.vcf
	tabix ${params.sequencing_run}.roi.filtered.norm.vcf.gz

	"""

}

// Annotate using VEP
process annotate_with_vep{

    cpus params.vep_cpus

    publishDir "${params.publish_dir}/annotated_vcf/", mode: 'copy'

    input:
    file(normalised_vcf) from normalised_vcf_channel

    output:
    set file("${params.sequencing_run}.roi.filtered.norm.anno.vcf.gz"), file("${params.sequencing_run}.roi.filtered.norm.anno.vcf.gz.tbi")  into annotated_vcf

    """
    vep \
    --verbose \
    --format vcf \
    --everything \
    --fork $params.vep_cpus \
    --species homo_sapiens \
    --assembly GRCh37 \
    --input_file $normalised_vcf \
    --output_file ${params.sequencing_run}.roi.filtered.norm.anno.vcf \
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

    bgzip ${params.sequencing_run}.roi.filtered.norm.anno.vcf
    tabix ${params.sequencing_run}.roi.filtered.norm.anno.vcf.gz

    """
}

// Calculate relatedness between samples
process calculate_relatedness {

    cpus params.relatedness_cpus

    publishDir "${params.publish_dir}/relatedness/", mode: 'copy'

    input:
    set file(vcf), file(vcf_index) from vcf_relatedness_channel

    output:
    file("${params.sequencing_run}.relatedness2")

    """
    vcftools --relatedness2 \
    --out $params.sequencing_run \
    --gzvcf $vcf

    """

}

// calculate intersample contamination
process calculate_contamination {

    cpus params.vcf_processing_cpus

    publishDir "${params.publish_dir}/contamination/", mode: 'copy'

    input:
    set val(id), file(bam), file(bam_index) from contamination_bams 

    output:
    file("${id}_contamination.selfSM") 

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


// create ped file
process create_ped {

    cpus params.small_task_cpus

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

    cpus params.small_task_cpus

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

    cpus params.vcf_processing_cpus

    publishDir "${params.publish_dir}/database_vcf/", mode: 'copy'

    input:
    set file(vcf), file(vcf_index) from vcf_database_channel 
    file(metadata) from meta_data_channel

    output:
    file("${params.sequencing_run}.roi.filtered.database.vcf") into sensitivity_calculation_vcf

    """
    # get header 

    zcat $vcf | grep "^##" > header.txt

    cat header.txt $metadata > ${params.sequencing_run}.roi.filtered.header.vcf

    zcat $vcf | grep -v "^##" >> ${params.sequencing_run}.roi.filtered.header.vcf

    # remove mt variants

    awk '\$1 !~ /^MT/ { print \$0 }' ${params.sequencing_run}.roi.filtered.header.vcf > ${params.sequencing_run}.roi.filtered.database.vcf

    """


}


// Use sambamba to generate the per base coverage
process generate_coverage_file{

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

	input:
    set val(id), file(bam), file(bam_index) from coverage_bams 

    output:
    set val(id), file("${id}_depth_of_coverage.gz"), file("${id}_depth_of_coverage.gz.tbi") into per_base_coverage_channel
    file("${id}_depth_of_coverage.sample_summary") into depth_summary

	"""
    gatk3 -T DepthOfCoverage \
    -R $reference_genome \
    -o ${id}_depth_of_coverage \
    -I $bam \
    -L $capture_bed \
    --countType COUNT_FRAGMENTS \
    --minMappingQuality 20 \
    --minBaseQuality 10 \
    -ct $params.min_coverage \
    --omitLocusTable \
    -rf MappingQualityUnavailable \
    -dt NONE

    sed 's/:/\t/g' ${id}_depth_of_coverage \
    | grep -v 'Locus' | sort -k1,1 -k2,2n | bgzip > ${id}_depth_of_coverage.gz

    tabix -b 2 -e 2 -s 1 ${id}_depth_of_coverage.gz

	"""
}

// split filtered vcf channel into four
per_base_coverage_channel.into{
    per_base_coverage_channel_gaps
    per_base_coverage_channel_sex
    per_base_coverage_channel_sensitivity
}


// Use coverage calculator to get gaps etc
process generate_gaps_files{

    params.small_task_cpus

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

	input:
    set val(id), file(depth_file), file(depth_file_index) from per_base_coverage_channel_gaps 

    output:
    set file("${id}_gaps.bed"), file("${id}_clinical_coverage_target_metrics.txt"),  file("${id}_clinical_coverage_gene_coverage.txt")

	"""
    get_coverage.sh \
    $capture_bed \
    $hotspots \
    $depth_file \
    $id \
    $params.min_coverage \
    $refseq

	"""
}


// use the alignment metrics file to calculate the sex
process calculate_sex{

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/calculated_sex/", mode: 'copy'

    input:
    set val(id), file(depth_file), file(depth_file_index) from per_base_coverage_channel_sex

    output:
    file("${id}_calculated_sex.txt")

    """
    calculate_sex.py --file $depth_file > ${id}_calculated_sex.txt

    """

}

// filter so we only get the coverage for the giab sample
per_base_coverage_channel_sensitivity.filter( {it[0] =~ /$params.giab_sample.*/ } ).set{ giab_per_base_coverage_channel }

// calculate sensitivity using giab
process calculate_sensitivity{

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/sensitivity/", mode: 'copy'

    input:
    set val(id), file(depth_file), file(depth_file_index) from giab_per_base_coverage_channel 
    file(vcf) from sensitivity_calculation_vcf

    output:
    file("${id}_sensitivity.txt")

    when:
    params.calculate_sensitivity == true 

    """
    calculate_sensitivity.sh $depth_file \
    $params.giab_sample \
    $vcf \
    $giab_baseline \
    $reference_genome_faidx \
    $giab_high_confidence_bed \
    $giab_reference_genome_sdf > ${id}_sensitivity.txt
    """


}

workflow.onComplete{

    if (workflow.success){

        ran_ok = "${params.sequencing_run} success!.\n"
    }
    else{

        ran_ok = "${params.sequencing_run} fail!.\n"

    }

    def newFile = new File("${params.publish_dir}/pipeline_complete.txt")

    newFile.createNewFile()
    newFile.append(ran_ok)

}