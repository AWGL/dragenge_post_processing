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

params.pipeline_version = "300"

capture_bed = file(params.capture_bed)
coverage_per_base_bed = file(params.coverage_per_base_bed)
coverage_report_bed = file(params.coverage_report_bed)
reference_genome = file(params.reference_genome)
reference_genome_index = file(params.reference_genome_index)
vep_cache = file(params.vep_cache)
clinvar = file(params.clinvar)
gnotate_gnomad = file(params.gnotate_gnomad)
gnotate_spliceai = file(params.gnotate_spliceai)
high_confidence_snps = file(params.high_confidence_snps)
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
  .set{ variables_ch }

Channel
  .fromFilePairs(params.vcf, flat: true)
  .ifEmpty { exit 1, "Cannot find a vcf file files matching: ${params.vcf}" }
  .set { raw_vcf_ch }

Channel
  .fromFilePairs(params.bams, flat: true)
  .ifEmpty { exit 1, "Cannot find any bam files files matching: ${params.bams}" }
  .set { original_bams_ch }


Channel
  .fromFilePairs(params.sv_vcf, flat: true)
  .set { original_sv_vcf_ch }

variables_ch.into{
    variables_ped_ch
    variables_config_ch
    variables_reporting_ch
    variables_old_variant_database_ch
}

original_bams_ch.into{
    original_bams_coverage_ch
    original_bams_contamination_ch
    original_bams_sensitivity_ch
}


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
    file(variables) from variables_ped_ch.collect()

    output:
    file("${params.sequencing_run}.ped") into (ped_config_ch, ped_reporting_ch)

    """
    create_ped.py --variables "*.variables" > ${params.sequencing_run}.ped
    """
}



// Create json for qiagen upload
process create_config_for_database {

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/database_config/", mode: 'copy'

    input:
    file(ped) from ped_config_ch
    file(variables) from variables_config_ch.collect()

    output:
    file("${params.sequencing_run}_config.csv")

    """
    make_database_config.py \
    --ped_file $ped \
    --variables "*.variables" \
    --pipeline_name dragenge_post_processing_$params.pipeline_version \
    --output_name ${params.sequencing_run}_config.csv
    """

}


// Merge the joint called per Chromosome VCFs and restrict to original ROI
process restrict_vcf_to_roi{

    cpus params.gatk_cpus
    scratch params.scratch 

    input:
    set val(id), file(vcf), file(vcf_index) from raw_vcf_ch 

    output:
    set file("${params.sequencing_run}_roi.vcf.gz"), file("${params.sequencing_run}_roi.vcf.gz.tbi") into roi_vcf_ch

    """
    bcftools view -R $capture_bed $vcf | bcftools sort -O v - > "${params.sequencing_run}_roi.vcf"

    bgzip ${params.sequencing_run}_roi.vcf
    tabix ${params.sequencing_run}_roi.vcf.gz
    """
}


// Apply quality filters to VCF
process quality_filter_vcf{

    cpus params.gatk_cpus
    scratch params.scratch 

    input:
    set file(vcf), file(vcf_idx) from roi_vcf_ch

    output:
    set file("${params.sequencing_run}_roi_qual.vcf.gz"), file("${params.sequencing_run}_roi_qual.vcf.gz.tbi") into (quality_filtered_vcf_annotation_ch,
                                                                                                                    quality_filtered_vcf_old_variant_database_ch, 
                                                                                                                    quality_filtered_vcf_sensitivity_ch )

    """
    quality_filter.sh \
    $vcf \
    $params.sequencing_run \
    $params.snp_qual \
    $params.indel_qual \
    $params.reference_genome \
    "$params.medium_java_options" \
    $params.min_dp_variant_reporting
    """

}

// Annotate using VEP and gnomad and SpliceAI
process normalise_annotate_with_vep_and_gnomad{

    cpus params.big_task_cpus
    scratch params.scratch 

    publishDir "${params.publish_dir}/annotated_vcf/", mode: 'copy'

    input:
    set file(quality_vcf), file(quality_vcf_idx) from quality_filtered_vcf_annotation_ch

    output:
    set file("${params.sequencing_run}_anno.vcf.gz"), file("${params.sequencing_run}_anno.vcf.gz.tbi") into (annotated_vcf_relatedness_ch, annotated_vcf_samples_ch, annotated_vcf_reporting_ch)
    file("${params.sequencing_run}_anno.vcf.gz.md5")

    """
    zcat $quality_vcf | vt decompose -s - | vt normalize -r $reference_genome - > ${params.sequencing_run}_norm.vcf

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
    --assembly $params.genome_build \
    --input_file ${params.sequencing_run}_norm.vcf \
    --output_file ${params.sequencing_run}_vep.vcf \
    --force_overwrite \
    --cache \
    --dir $vep_cache \
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

    bgzip ${params.sequencing_run}_vep.vcf
    tabix ${params.sequencing_run}_vep.vcf.gz

    slivar expr --gnotate $gnotate_spliceai -o ${params.sequencing_run}_vep.spliceai.vcf -v ${params.sequencing_run}_vep.vcf.gz
    bgzip ${params.sequencing_run}_vep.spliceai.vcf
    tabix ${params.sequencing_run}_vep.spliceai.vcf.gz

    slivar expr --gnotate $gnotate_gnomad -o ${params.sequencing_run}_anno.vcf -v ${params.sequencing_run}_vep.spliceai.vcf.gz
    bgzip ${params.sequencing_run}_anno.vcf
    tabix ${params.sequencing_run}_anno.vcf.gz

    md5sum ${params.sequencing_run}_anno.vcf.gz > ${params.sequencing_run}_anno.vcf.gz.md5
    """
}

// Calculate relatedness between samples
process calculate_relatedness {

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/relatedness/", mode: 'copy'

    input:
    set file(vcf), file(vcf_index) from annotated_vcf_relatedness_ch

    output:
    file("${params.sequencing_run}.relatedness2") into relatedness_ch

    """
    vcftools --relatedness2 \
    --out $params.sequencing_run \
    --gzvcf $vcf
    """
}


/*
========================================================================================
Variant Reporting
========================================================================================
*/

// Extract sample names from vcf
process get_sample_names_from_vcf{

    cpus params.small_task_cpus

    input:
    set file(vcf), file(vcf_index) from annotated_vcf_samples_ch

    output:
    file("${params.sequencing_run}_sample_names.txt") into sample_names_from_vcf_ch

    """
    bcftools query -l $vcf > ${params.sequencing_run}_sample_names.txt
    """

}

// Create a channel with sample names as the values
sample_names_from_vcf_ch.splitCsv(header:['col1']).map{ row-> tuple(row.col1)}.set { samples_ch }

// Create variant reports csvs for non mitochondrial variants
process create_variant_reports {

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/variant_reports/", mode: 'copy'

    input:
    set file(vcf), file(vcf_index) from annotated_vcf_reporting_ch
    each sample_names from samples_ch
    file(ped) from ped_reporting_ch
    file(variables) from variables_reporting_ch.collect()

    output:
    file("${params.sequencing_run}_${sample_names[0]}_variant_report.csv") into variant_report_channel

    """
    . ${sample_names[0]}.variables
    germline_variant_reporter.py \
    --vcf $vcf \
    --proband_id ${sample_names[0]} \
    --ped $ped \
    --gene_list $params.gene_panel \
    --gnomad_ad $params.gnomad_ad \
    --gnomad_r $params.gnomad_r \
    --min_dp $params.min_dp_variant_reporting \
    --min_gq $params.min_gq_variant_reporting \
    --max_parental_alt_ref_ratio $params.max_parental_alt_ref_ratio \
    --output ${params.sequencing_run}_${sample_names[0]}_variant_report.csv \
    --worklist \$worklistId \
    --splice_ai $params.splice_ai_cutoff \
    """
}


// make a VCF for the old variant database
process create_vcf_for_old_variant_database{

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/old_variant_database_vcf/", mode: 'copy'

    input:
    set file(vcf), file(vcf_idx) from quality_filtered_vcf_old_variant_database_ch
    file(variables) from variables_old_variant_database_ch.collect()

    output:
    file("${params.sequencing_run}_old_variant_database.vcf")

    when:
    params.create_old_database_vcf

    """
    # loop through variables and make into a metadata file
    for i in ${variables.collect { "$it" }.join(" ")}; do
        # source variables
        . \$i
        
        echo \\#\\#SAMPLE\\=\\<ID\\=\$sampleId,Tissue\\=Germline,WorklistId\\="\$worklistId",SeqId\\="\$seqId",Assay\\="\$panel",PipelineName\\=dragenge_post_processing,PipelineVersion\\="$params.pipeline_version"\\,RawSequenceQuality=0,PercentMapped=0,ATDropout=0,GCDropout=0,MeanInsertSize=0,SDInsertSize=0,DuplicationRate=0,TotalReads=0,PctSelectedBases=0,MeanOnTargetCoverage=0,PctTargetBasesCt=0,EstimatedContamination=0,GenotypicGender=None,TotalTargetedUsableBases=0,RemoteVcfFilePath=None,RemoteBamFilePath=None\\> >> ${params.sequencing_run}_meta.txt
    done

    # get header 
    zcat $vcf | grep "^##" > header.txt
    cat header.txt ${params.sequencing_run}_meta.txt > ${params.sequencing_run}.roi.filtered.header.vcf

    zcat $vcf | grep -v "^##" >> ${params.sequencing_run}.roi.filtered.header.vcf


    fix_ploidy.py --input_vcf ${params.sequencing_run}.roi.filtered.header.vcf --output_vcf ${params.sequencing_run}.roi.filtered.header.fixed.vcf

    # remove mt variants
    awk '\$1 !~ /^MT/ { print \$0 }' ${params.sequencing_run}.roi.filtered.header.fixed.vcf > ${params.sequencing_run}_old_variant_database.vcf

    """
}

/*
========================================================================================
Depth & QC - Sex Calculations, Coverage Reports and Contamination
========================================================================================
*/

// calculate intersample contamination
process calculate_contamination {

    cpus params.medium_task_cpus

    publishDir "${params.publish_dir}/contamination/", mode: 'copy'

    input:
    set val(id), file(bam), file(bam_index) from original_bams_contamination_ch 

    output:
    file("${id}_contamination.selfSM") 

    script:
    if (params.genome_build == 'GRCh37')

    """
    gatk3 -T SelectVariants \
    -R $reference_genome \
    --variant $high_confidence_snps \
    -o high_confidence_snps.vcf \
    -selectType SNP \
    -restrictAllelesTo BIALLELIC \
    -L $capture_bed \
    -XL X \
    -XL Y \
    -XL MT \
    -env \
    -ef \
    -dt NONE

    verifyBamID \
    --vcf high_confidence_snps.vcf \
    --bam $bam \
    --out ${id}_contamination \
    --verbose \
    --ignoreRG \
    --chip-none \
    --minMapQ $params.min_mapq_contamination \
    --maxDepth $params.max_depth_contamination \
    --precise
    """

    else if (params.genome_build == 'GRCh38')

    """
    gatk3 -T SelectVariants \
    -R $reference_genome \
    --variant $high_confidence_snps \
    -o high_confidence_snps.vcf \
    -selectType SNP \
    -restrictAllelesTo BIALLELIC \
    -L $capture_bed \
    -XL chrX \
    -XL chrY \
    -XL chrM \
    -env \
    -ef \
    -dt NONE

    verifyBamID \
    --vcf high_confidence_snps.vcf \
    --bam $bam \
    --out ${id}_contamination \
    --verbose \
    --ignoreRG \
    --chip-none \
    --minMapQ $params.min_mapq_contamination \
    --maxDepth $params.max_depth_contamination \
    --precise
    """

    else
    error "Invalid genome_build mode: ${genome_build}"
}


// Get per base coverage using GATK3
process get_per_base_coverage{

    cpus params.medium_task_cpus

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

    input:
    set val(id), file(bam), file(bam_index) from original_bams_coverage_ch

    output:
    set val(id), file("${id}_depth_of_coverage.csv.gz"), file("${id}_depth_of_coverage.csv.gz.tbi") into (sex_reporting_ch, coverage_reporting_ch, sensitivity_coverage_ch,  per_base_coverage_summary_ch)

    """
    gatk3 $params.medium_java_options -T DepthOfCoverage \
    -R $reference_genome \
    -o ${id}_depth_of_coverage \
    -I $bam \
    -L $coverage_per_base_bed \
    --countType COUNT_FRAGMENTS \
    --minMappingQuality $params.min_mapping_quality_coverage \
    --minBaseQuality $params.min_base_quality_coverage\
    -ct 20 \
    --omitLocusTable \
    -rf MappingQualityUnavailable \
    -dt NONE

    sed 's/:/\t/g' ${id}_depth_of_coverage \
    | grep -v 'Locus' | sort -k1,1 -k2,2n | bgzip > ${id}_depth_of_coverage.csv.gz

    tabix -b 2 -e 2 -s 1 ${id}_depth_of_coverage.csv.gz

    """

}

// Generate a depth/coverage summary report e.g. mean depth etc
process create_depth_summary {

    cpus params.medium_task_cpus

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

    input:
    set val(key), file(coverage), file(coverage_idx) from per_base_coverage_summary_ch

    output:
    file("${key}.depth_summary") into depth_summary_ch

    """
    zcat $coverage | awk '{print \$1"\t"\$2-1"\t"\$2"\t"\$3}' > "$key"_coverage_bed.bed

    bedtools intersect -a "$key"_coverage_bed.bed -b $capture_bed > "$key"_coverage_bed_final.bed

    awk '{print \$1"\t"\$2"\t"\$4"\t"\$4}' "$key"_coverage_bed_final.bed > "$key"_per_base_final.csv

    create_depth_summary.py \
    --input_file ${key}_per_base_final.csv \
    --output_file ${key}.depth_summary \
    ${params.depth_thresholds.collect { "$it" }.join(" ")}
    """
}


// Calculate relatedness between samples
process calculate_sex {

    cpus params.medium_task_cpus

    publishDir "${params.publish_dir}/sex/", mode: 'copy'

    input:
    set val(id), file(coverage), file(coverage_idx) from sex_reporting_ch

    output:
    file("${id}_calculated_sex.txt")

    """
    calculate_sex.py \
    --file $coverage \
    --female_theshold $params.female_threshold \
    --male_theshold $params.male_threshold \
    --genome_build $params.genome_build > ${id}_calculated_sex.txt
    """
}


// Generate a report of the depth/coverage over the regions
process get_region_coverage_data {

    cpus params.medium_task_cpus

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

    input:
    set val(id), file(coverage_file), file(coverage_file_index) from coverage_reporting_ch

    output:
    set val(id), file("${id}_region_coverage_data.csv") into region_coverage_data
    set val(id), file("${id}_region_coverage_summary.csv") into region_coverage_summary_data
    set val(id), file("${id}_gene_coverage_summary.csv") into gene_summary_data
    set val(id), file("${id}_gaps_*x.csv") into gaps_data

    """
    get_coverage_region_stats.py \
    --input_file $coverage_file \
    --bed $coverage_report_bed \
    --sample_id ${id} \
    --output_name ${id} \
    --calculate_gene_stats_gaps true \
    ${params.depth_thresholds.collect { "$it" }.join(" ")}
    """
}

// merge coverage data - experimental for variant database
process merge_coverage_data {

    cpus params.medium_task_cpus

    publishDir "${params.publish_dir}/coverage/", mode: 'copy'

    input:
    file(coverage_files) from gene_summary_data.collect()

    output:
    file("${params.sequencing_run}_coverage.tar")

    """
    for i in ./*_gene_coverage_summary.csv;  do cp \$i \$i.copy && gzip \$i.copy; done

    tar -cvf ${params.sequencing_run}_coverage.tar *_gene_coverage_summary.csv.copy.gz
    """
}



// filter so we only get the coverage for the giab sample
original_bams_sensitivity_ch.filter( {it[0] =~ /$params.giab_sample.*/ } ).set{ giab_original_bams_sensitivity_ch }

// calculate sensitivity using giab
process calculate_sensitivity{

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/sensitivity/", mode: 'copy'

    input:
    set val(id), file(bam), file(bam_index) from giab_original_bams_sensitivity_ch 
    set file(vcf), file(vcf_idx) from quality_filtered_vcf_sensitivity_ch

    output:
    file("${id}_sensitivity.txt")

    when:
    params.calculate_sensitivity == true 

    """
    calculate_sensitivity.sh $bam \
    $params.giab_sample \
    $vcf \
    $giab_baseline \
    $reference_genome_index \
    $giab_high_confidence_bed \
    $giab_reference_genome_sdf \
    $capture_bed \
    $params.min_mapping_quality_coverage \
    $params.min_mapping_quality_coverage \
    $reference_genome \
    "$params.medium_java_options" \
    > ${id}_sensitivity.txt
    """
}


process annotate_cnvs_svs {

    cpus params.sv_cpus

    publishDir "${params.publish_dir}/annotated_sv_vcf/", mode: 'copy'

    when:
    params.annotate_svs

    input:
    set val(id), file(vcf), file(vcf_index) from original_sv_vcf_ch

    output:
    set file("${params.sequencing_run}.sv.vep.vcf.gz"), file("${params.sequencing_run}.sv.vep.vcf.gz.tbi")
    file("${params.sequencing_run}.sv.vep.vcf.gz.md5")


    when:
    params.annotate_svs

    """
    bcftools norm -m - $vcf > ${params.sequencing_run}.sv.norm.vcf

    vep --verbose \
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
    --variant_class \
    --check_existing \
    --species homo_sapiens \
    --assembly $params.genome_build \
    --input_file ${params.sequencing_run}.sv.norm.vcf \
    --output_file ${params.sequencing_run}.sv.vep.vcf \
    --cache \
    --dir $vep_cache \
    --fasta $reference_genome \
    --offline \
    --cache_version $params.vepversion \
    --no_escape \
    --shift_hgvs 1 \
    --vcf \
    --refseq \
    --flag_pick \
    --exclude_predicted \
    --max_sv_size $params.sv_max_size \
    --buffer_size $params.sv_vep_buffer_size

    bgzip ${params.sequencing_run}.sv.vep.vcf
    tabix ${params.sequencing_run}.sv.vep.vcf.gz

    md5sum ${params.sequencing_run}.sv.vep.vcf.gz > ${params.sequencing_run}.sv.vep.vcf.gz.md5

    """
}

// Create marker once complete
workflow.onComplete{

	if (workflow.success){

		ran_ok = "${params.sequencing_run} success!.\n"
	}
	else{

		ran_ok = "${params.sequencing_run} fail!.\n"

	}

	def newFile = new File("${params.publish_dir}/post_processing_finished.txt")

	newFile.createNewFile()
	newFile.append(ran_ok)
}

