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

capture_bed = file(params.capture_bed)
coverage_bed = file(params.coverage_bed)
reference_genome = file(params.reference_genome)
reference_genome_index = file(params.reference_genome_index)
vep_cache = file(params.vep_cache)
clinvar = file(params.clinvar)
vep_cache_mt = file(params.vep_cache_mt)
mitomap = file(params.mitomap_vcf)
gnotate_gnomad = file(params.gnotate_gnomad)
gnotate_spliceai = file(params.gnotate_spliceai)

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
    gene_coverage_bam
    contamination_bams
}


variables_channel.into{
    variables_ped
    variables_config
}

raw_vcf.into{
    raw_vcf_relatedness
    raw_vcf_roi
}



// Chromosomes for when we do VEP in parallel
chromosome_ch = Channel.from('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT' )



/*
========================================================================================
Main pipeline
========================================================================================
*/

// Create PED file
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

ped_channel.into{
    ped_channel_config
    ped_channel_upd
}


// Create json for qiagen upload
process create_config_for_database {

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/database_config/", mode: 'copy'

    input:
    file(ped) from ped_channel_config
    file(variables) from variables_config.collect()

    output:
    file("${params.sequencing_run}_config.csv")

    """
    make_database_config.py \
    --ped_file $ped \
    --variables "*.variables" \
    --pipeline_name dragenwgs_post_processing_$params.version \
    --output_name ${params.sequencing_run}_config.csv
    """

}

// Use bcftools to restrict vcf to region of interest
process restrict_vcf_to_roi{

    cpus params.vcf_processing_cpus

	input:
    set val(id), file(vcf), file(vcf_index) from raw_vcf_roi 

    output:
    set val(id), file("${params.sequencing_run}.roi.vcf.gz"), file("${params.sequencing_run}.roi.vcf.gz.tbi") into roi_vcf_channel

	"""
	bcftools view -R $capture_bed $vcf > ${params.sequencing_run}.roi.vcf

	bgzip ${params.sequencing_run}.roi.vcf
	tabix ${params.sequencing_run}.roi.vcf.gz
	"""
}


// Split multiallelics and normalise also annotate with spliceai + gnomad
process split_multiallelics_and_normalise_and_annotate{

    cpus params.vcf_processing_cpus

    input:
    set val(id), file(vcf), file(vcf_index) from roi_vcf_channel 
    file(reference_genome_index)
    file(gnotate_spliceai)
    file(gnotate_gnomad)

    output:
    set val(id), file("${params.sequencing_run}.norm.vcf.gz"), file("${params.sequencing_run}.norm.vcf.gz.tbi") into normalised_vcf_channel

    """
    zcat  < $vcf | vt decompose -s - | vt normalize -r $reference_genome - > ${params.sequencing_run}.norm.inter.vcf
    bgzip ${params.sequencing_run}.norm.inter.vcf
    tabix ${params.sequencing_run}.norm.inter.vcf.gz

    slivar expr --gnotate $gnotate_spliceai -o ${params.sequencing_run}.norm.inter2.vcf -v ${params.sequencing_run}.norm.inter.vcf.gz

    bgzip ${params.sequencing_run}.norm.inter2.vcf
    tabix ${params.sequencing_run}.norm.inter2.vcf.gz

    slivar expr --gnotate $gnotate_gnomad -o ${params.sequencing_run}.norm.vcf -v ${params.sequencing_run}.norm.inter2.vcf.gz

    bgzip ${params.sequencing_run}.norm.vcf
    tabix ${params.sequencing_run}.norm.vcf.gz

    """

}


normalised_vcf_channel.into{
    for_splitting_channel
    for_mito_vcf_channel
}

// Split the normalised vcf up per chromosome for parallel processing
process split_vcf_per_chromosome{

    cpus params.vcf_processing_cpus

    input:
    set val(id), file(normalised_vcf), file(normalised_vcf_index) from for_splitting_channel
    each chromosome from chromosome_ch

    output:
    set val(chromosome), file("${params.sequencing_run}.norm.${chromosome}.vcf.gz"), file("${params.sequencing_run}.norm.${chromosome}.vcf.gz.tbi")  into for_vep_channel

    """
    bcftools view -r $chromosome $normalised_vcf > ${params.sequencing_run}.norm.${chromosome}.vcf

    bgzip ${params.sequencing_run}.norm.${chromosome}.vcf
    tabix ${params.sequencing_run}.norm.${chromosome}.vcf.gz

    """
}

// Annotate using VEP 
process annotate_with_vep{

    cpus params.vep_cpus

    input:
    set val(chromosome), file(normalised_vcf), file(normalised_vcf_index) from for_vep_channel
    file reference_genome_index
    file vep_cache

    output:
    file("${params.sequencing_run}.norm.${chromosome}.anno.vcf") into annotated_vcf_per_chromosome

    """
    vep \
    --verbose \
    --format vcf \
    --hgvs \
    --symbol \
    --numbers \
    --domains \
    --regulatory \
    --canonical \
    --protein \
    --biotype \
    --uniprot \
    --tsl \
    --appris \
    --gene \
    --variant_class \
    --clin_sig_allele 0 \
    --check_existing \
    --fork $params.vep_cpus \
    --species homo_sapiens \
    --assembly GRCh37 \
    --input_file $normalised_vcf \
    --output_file ${params.sequencing_run}.norm.${chromosome}.anno.vcf \
    --force_overwrite \
    --cache \
    --dir  $vep_cache \
    --fasta $reference_genome \
    --offline \
    --cache_version $params.vepversion \
    --no_escape \
    --shift_hgvs 1 \
    --vcf \
    --refseq \
    --flag_pick \
    --pick_order biotype,canonical,appris,tsl,ccds,rank,length \
    --exclude_predicted \
    --custom ${clinvar},clinvar,vcf,exact,0,CLNSIG,CLNSIGCONF

    """
}

// Merge per chromosome vcfs back into one vcf
process merge_annotated_vcfs{

    publishDir "${params.publish_dir}/annotated_vcf/", mode: 'copy'

    input:
    file(vcfs) from annotated_vcf_per_chromosome.collect()

    output:
    set file("${params.sequencing_run}.norm.anno.vcf.gz"), file("${params.sequencing_run}.norm.anno.vcf.gz.tbi") into annotated_vcf

    """
    bcftools concat ${vcfs.collect { "$it " }.join()} > ${params.sequencing_run}.norm.anno.unsorted.vcf

    bcftools sort -T $params.tmp_dir ${params.sequencing_run}.norm.anno.unsorted.vcf > ${params.sequencing_run}.norm.anno.vcf

    bgzip ${params.sequencing_run}.norm.anno.vcf
    tabix ${params.sequencing_run}.norm.anno.vcf.gz

    """
}


// Calculate relatedness between samples
process calculate_relatedness {

    cpus params.relatedness_cpus

    publishDir "${params.publish_dir}/relatedness/", mode: 'copy'

    input:
    set val(id), file(vcf), file(vcf_index) from raw_vcf_relatedness

    output:
    file("${params.sequencing_run}.relatedness2")

    """
    vcftools --relatedness2 \
    --out $params.sequencing_run \
    --gzvcf $vcf
    """
}

// Extract mitochrondrial variants and annotate with VEP - use ensembl transcripts
process get_mitochondrial_variant_and_annotate{

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/annotated_mitochrondrial_vcf/", mode: 'copy'

    input:
    set val(id), file(normalised_vcf), file(normalised_vcf_index) from for_mito_vcf_channel
    file reference_genome_index
    file vep_cache_mt

    output:
    set val(id), file("${params.sequencing_run}.mito.anno.vcf.gz"), file("${params.sequencing_run}.mito.anno.vcf.gz.tbi") into mito_vcf_channel

    """
    bcftools view -r MT $normalised_vcf > ${params.sequencing_run}.mito.vcf
    bgzip ${params.sequencing_run}.mito.vcf
    tabix ${params.sequencing_run}.mito.vcf.gz

    vep \
    --verbose \
    --format vcf \
    --everything \
    --fork $params.vep_cpus \
    --species homo_sapiens \
    --assembly GRCh37 \
    --input_file ${params.sequencing_run}.mito.vcf.gz \
    --output_file ${params.sequencing_run}.mito.anno.vcf \
    --force_overwrite \
    --cache \
    --dir  $vep_cache_mt \
    --fasta $reference_genome \
    --offline \
    --cache_version $params.vepversion_mt \
    --no_escape \
    --shift_hgvs 1 \
    --vcf \
    --merged \
    --flag_pick \
    --pick_order biotype,canonical,appris,tsl,ccds,rank,length \
    --exclude_predicted \
    --custom ${mitomap},mitomap,vcf,exact,0,AF,AC \
    --custom ${clinvar},clinvar,vcf,exact,0,CLNSIG,CLNSIGCONF
    bgzip ${params.sequencing_run}.mito.anno.vcf
    tabix ${params.sequencing_run}.mito.anno.vcf.gz
    """
}

// Use sambamba to generate the per base coverage
process generate_coverage_file{

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

    input:
    set val(id), file(bam), file(bam_index) from coverage_bams 

    output:
    set val(id), file("${id}_depth_of_coverage.csv.gz"), file("${id}_depth_of_coverage.csv.gz.tbi") into per_base_coverage_channel

    """
    sambamba depth base \
    -L $capture_bed \
    --min-coverage=0 \
    --min-base-quality=10 \
    --filter 'mapping_quality > 20 and not duplicate and not failed_quality_control' $bam > ${id}_depth_of_coverage.csv

    grep -v REF ${id}_depth_of_coverage.csv | bgzip > ${id}_depth_of_coverage.csv.gz

    tabix -b 2 -e 2 -s 1 ${id}_depth_of_coverage.csv.gz


    """
}

// Use sambamba to generate the per base coverage
process get_region_coverage_data{

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

    input:
    set val(id), file(coverage_file), file(coverage_file_index) from  per_base_coverage_channel

    output:
    set val(id), file("${id}_region_coverage_data.csv") into region_coverage_data

    """
    get_coverage_region_stats.py \
    --input_file $coverage_file \
    --bed $capture_bed \
    --sample_id ${id} \
    --output_name ${id}_region_coverage_data.csv \
    20 160


    """
}

// Merge per chromosome vcfs back into one vcf
process collect_region_coverage_data_and_annotate{

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

    input:
    file(vcfs) from region_coverage_data.collect()

    output:
    set file("${params.sequencing_run}_merged_coverage_data.csv") into merged_coverage_data

    """
    bcftools concat ${vcfs.collect { "$it " }.join()} > ${params.sequencing_run}.norm.anno.unsorted.vcf

    bcftools sort -T $params.tmp_dir ${params.sequencing_run}.norm.anno.unsorted.vcf > ${params.sequencing_run}.norm.anno.vcf

    bgzip ${params.sequencing_run}.norm.anno.vcf
    tabix ${params.sequencing_run}.norm.anno.vcf.gz

    """
}