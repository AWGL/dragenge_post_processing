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

params.bams = '/media/joseph/Storage/test_data/IlluminaTruSightOne/*/*{.bam,.bam.bai}'
params.vcf =  '/media/joseph/Storage/test_data/IlluminaTruSightOne/*.hard-filtered{.vcf.gz,.vcf.gz.tbi}'
params.variables = '/media/joseph/Storage/test_data/IlluminaTruSightOne/*/*.variables'
params.reference_genome =  '/media/joseph/Storage/genomic_resources/bwa/human_g1k_v37.fasta'
params.sequence_dict = '/media/joseph/Storage/genomic_resources/bwa/human_g1k_v37.dict'
params.capture_bed = 'config/IlluminaTruSightOne_ROI_b37.bed'
params.coverage_bed = 'config/IlluminaTruSightOne_ROI_b37_coverage.bed'
params.coverage_groups = 'config/IlluminaTruSightOne_ROI_b37_coverage_groups.txt'
params.gnomad_exomes = '/media/joseph/Storage/genomic_resources/gnomad/gnomad.exomes.r2.0.1.sites.vcf.gz'
params.gnomad_genomes = '/media/joseph/Storage/genomic_resources/gnomad/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz'
params.vep_threads = 1
params.vep_cache = '/media/joseph/Storage/genomic_resources/vep_caches/vep'
params.sequencing_run = '191010_D00501_0366_BH5JWHBCX3'
params.coverage_calculator_path = '/home/joseph/Documents/dragen/coverage/CoverageCalculatorPy/CoverageCalculatorPy.py'
params.min_coverage = 20
params.pipeline_version = "1"
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
  .set { filtered_vcf }

Channel
  .fromFilePairs(params.bams, flat: true) 
  .set { original_bams }

filtered_vcf.into {
  filtered_vcf_for_relatedness
  filtered_vcf_for_annotation
}


/*
========================================================================================
Main pipeline
========================================================================================
*/


// Use bcftools to restrict vcf to region of interest and remove chr from vcf chromosomes
process restrict_vcf_to_roi{

	input:
    set val(id), file(vcf), file(vcf_index) from filtered_vcf_for_annotation 

    output:
    set file("${params.sequencing_run}.hard-filtered.roi.vcf.gz"), file("${params.sequencing_run}.hard-filtered.roi.vcf.gz.tbi") into roi_vcf_channel

	"""
	bcftools view -R $capture_bed ${params.sequencing_run}.hard-filtered.vcf.gz > ${params.sequencing_run}.hard-filtered.roi.vcf

	bgzip ${params.sequencing_run}.hard-filtered.roi.vcf
	tabix ${params.sequencing_run}.hard-filtered.roi.vcf.gz
	"""
}

roi_vcf_channel.into {
  roi_vcf_for_annotation_channel
  roi_vcf_for_database_channel
}



// split multiallelics and normalise
process split_multiallelics_and_normalise{

	input:
	set file(vcf), file(vcf_index) from roi_vcf_for_annotation_channel 

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

    input:
    set val(id), file(vcf), file(vcf_index) from filtered_vcf_for_relatedness

    output:
    file("${id}.relatedness2") into relatedness_channel

    """
    vcftools --relatedness2 \
    --out ${id} \
    --gzvcf $vcf

    """

}


// Variant database needs extra metadata from variables files
process collect_metadata_for_vcf{

    input:
    file(variables) from variables_channel.collect()

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


process prepare_vcf_for_database {

    input:
    set file(vcf), file(vcf_index) from roi_vcf_for_database_channel 
    file(metadata) from meta_data_channel

    output:
    set file("${params.sequencing_run}.hard-filtered.roi.database.vcf") into database_vcf

    """
    # add meta data

    # remove mt variants



    """


}


// Use sambamba to generate the per base coverage
process generate_coverage_file{

	input:
    set val(id), file(bam), file(bam_index) from original_bams 

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