#!/usr/bin/env python

import argparse
import pandas as pd
from pyvariantfilter.family import Family
from pyvariantfilter.family_member import FamilyMember
from pyvariantfilter.variant_set import VariantSet

"""
Script using pyvariantfilter to filter a multisample VCF by criteria such as:

1) Inheritance
2) Consequence
3) Gene Panel
4) Presence in whitelist

Usage:
    
python germline_variant_reporter.py \
    --vcf $vcf \
    --proband_id ${sample_names[0]} \
    --ped $ped \
    --gene_list $gene_panel \
    --gnomad_ad $gnomad_ad \
    --gnomad_r $gnomad_r \
    --minqual_snps $snp_qual \
    --minqual_indels $indel_qual \
    --min_dp $min_dp \
    --min_gq $min_gq \
    --output $output_variant_report.csv \
    --worklist $worklist_id \
    --whitelist $whitelist \
    --apply_panel 

"""

parser = argparse.ArgumentParser(description='Create a text report with priority')
parser.add_argument('--vcf', type=str, nargs=1, required=True,
				help='The path to the VCF file. Must be bgzipped and tabixed.')
parser.add_argument('--proband_id', type=str, nargs=1, required=True,
				help='The Sample ID of the proband in the VCF.')
parser.add_argument('--ped', type=str, nargs=1, required=True,
				help='A ped file describing the family relationships.')
parser.add_argument('--gene_list', type=str, nargs=1, required=True,
				help='A list of genes to use a virtual panel. Should be a CSV (comma seperated) with the gene names as the first column with a header line.')
parser.add_argument('--gnomad_ad', type=float, nargs=1, required=True,
				help='The maximum gnomad frequency to filter on for dominant variants.')
parser.add_argument('--gnomad_r', type=float, nargs=1, required=True,
				help='The maximum gnomad frequency to filter on for reccessive variants.')
parser.add_argument('--min_dp', type=int, nargs=1, required=True,
				help='The minimum depth for variants.')
parser.add_argument('--min_gq', type=int, nargs=1, required=True,
				help='The minimum GQ for variants.')
parser.add_argument('--max_parental_alt_ref_ratio', type=float, nargs=1, required=True,
				help='The maximum alt / ref ratio the parent can have and still be called de novo.')
parser.add_argument('--output', type=str, nargs=1, required=True, help='The output name')
parser.add_argument('--worklist', type=str, nargs=1, required=True, help='Worklist ID for Audit purposes. Can be used as a comment.')
parser.add_argument('--apply_panel', action='store_true', help='Whether to apply a virtual pannel?')
parser.add_argument('--splice_ai', type=float, nargs=1, required=True,
				help='Variants with spliceai value above this are kept in first round filter.')

args = parser.parse_args()

vcf = args.vcf[0]
proband_id = args.proband_id[0]
ped = args.ped[0]
gene_list = args.gene_list[0]
gnomad_ad = args.gnomad_ad[0]
gnomad_r = args.gnomad_r[0]
min_dp = args.min_dp[0]
min_gq = args.min_gq[0]
max_parental_alt_ref_ratio = args.max_parental_alt_ref_ratio[0]
output_name = args.output[0]
worklist = args.worklist[0]
apply_panel = args.apply_panel
splice_ai_cutoff = args.splice_ai[0]

# get maximum of both population frequencies to initially filter on. We have to do it in two stages to avoid filtering compound hets out.
initial_af = max(gnomad_ad, gnomad_r)

# define a few functions to help us out later on

def fix_genotype(df, column_key):
	"""
	Convert genotype in format G/A to HET etc
	"""

	genotype = df[column_key]
	variant  = df['variant_id']

	alt = variant.split('>')[1]

	if '/' in genotype:

		genotype = genotype.split('/')

	elif '|' in genotype:

		genotype = genotype.split('|')

	else:

		raise Exception('weird genotype')

	if genotype.count('.') == 2:

		return 'MISSING'

	elif genotype.count(alt) == 2:

		return 'HOM_ALT'

	elif genotype.count(alt) == 1:

		return 'HET'

	elif genotype.count(alt) == 0:

		return 'HOM_REF'

	else:

		return 'UNKNOWN'

def pathogenic_in_clinvar(clinvar_vep, clinvar_custom, clinvar_conflicting):
	"""
	Check multiple fields to see if pathogenic in clinvar.

	We check:

	1) CLIN_SIG - added by VEP
	2) CLIN_SIG - added by our own clinvar annnotation.
	3) CLINVAR_CONFLICTING - added by our own clinvar annotation.

	We check all three for presense of string pathogenic.


	"""
	for anno in clinvar_vep:
			

		if 'pathogenic' in anno.lower():

			return True

	for anno in clinvar_custom:

		# ignore this as we will pick it up via clinvar_conflicting
		if anno == 'Conflicting_interpretations_of_pathogenicity':

			pass

		elif 'pathogenic' in anno.lower():

			return True

	for anno in clinvar_conflicting:

		if 'pathogenic' in anno.lower():

			return True


	return False

def passes_initial_filter(variant, proband_id, gene_dict, min_gq, min_dp):
	"""
	An initial filter to 
	
	We import if the variant passes quality filtering and is below 1% in gnomad exomes and gnomad genomes AND
	
	a) Is listed as pathogenic in clinvar OR
	b) Has a a relevant consequence

	OR is in the whitelist
	
	"""

	# If the proband has the variant and we pass the genotype and variant level filters
	if variant.has_alt(proband_id) and variant.passes_gt_filter(proband_id, min_gq=min_gq, min_dp=min_dp):

		# ignore mt variants for now
		if variant.chrom == 'MT':

			return False

		else:

			if not variant.passes_filter():

				return False


		# The filter_on_numerical_transcript_annotation_lte() function allows us to filter on numerical values 
		# we can set different cutoffs for different variant types. For example ad_het is variants in which the 
		# proband is heterozygous on an autosome. In this case we get two boolean values describing whether the 
		# variant is below x% in the gnomad genomes and gnomad exomes datasets.
		freq_gnomad = variant.filter_on_numerical_info_annotation_lte(annotation_key='gnomad_popmax_af',
																						  ad_het=initial_af,
																						  ad_hom_alt=initial_af,
																						  x_male =initial_af,
																						  x_female_het=initial_af,
																						  x_female_hom=initial_af,
																						  compound_het=initial_af,
																						  y=initial_af,
																						  mt=initial_af,
																						  zero_values=['.', '', None, -1]
																						  )

		# if we haven't bothered with a panel
		if gene_dict == None:

			variant.info_annotations['in_panel'] = False
			in_panel = False

		else:

			genes = variant.get_genes(feature_key='SYMBOL')

			in_panel = False

			for gene in genes:

				if gene in gene_dict:

					# its in the panel so set annotation in_panel to True
					variant.info_annotations['in_panel'] = True
					in_panel = True
					break

		if freq_gnomad:
			
			return True
		
	return False

def passes_final_filter_trio(variant, compound_het_dict , inheritance, min_gq, min_dp ,max_parental_alt_ref_ratio):
	"""
	We do the filtering in 2 steps so we can work out which are comp hets.

	This function allows us to keep variants matching an inheritance pattern.

	"""

	freq_gnomad = variant.filter_on_numerical_info_annotation_lte(annotation_key='gnomad_popmax_af',
																  ad_het=gnomad_ad,
																  ad_hom_alt=initial_af,
																  x_male =initial_af,
																  x_female_het=initial_af,
																  x_female_hom=initial_af,
																  compound_het=initial_af,
																  y=initial_af,
																  mt=initial_af,
																  compound_het_dict=compound_het_dict,
																  zero_values=['.', '', None, -1]
																  )

	# Coopt the get_genes() function to get the clinvar annotation VEP field.
	clinvar_vep = variant.get_genes(feature_key='CLIN_SIG')
	clinvar_custom = variant.get_genes('clinvar_CLNSIG')
	clinvar_conflicting = variant.get_genes('clinvar_CLNSIGCONF')

	is_path_in_clinvar = pathogenic_in_clinvar(clinvar_vep, clinvar_custom, clinvar_conflicting)
			
	# If the variant is below 1% and pathogenic in clinvar then import
	if freq_gnomad and is_path_in_clinvar:
		
		return True

	
	# Get variants which match certain inheritance models
	if variant.matches_inheritance_model(inheritance,
										   compound_het_dict,
										   min_parental_gq_dn = min_gq,
										   min_parental_depth_dn = min_dp,
										   min_parental_gq_upi = min_gq,
										   min_parental_depth_upi = min_dp,
										   max_parental_alt_ref_ratio_dn = max_parental_alt_ref_ratio
										   ) and freq_gnomad:
		
		return True
		
	return False


# read ped into df
ped_df = pd.read_csv(ped, sep='\t', names=['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'affected'])

if apply_panel == True:

	# read gene list to dict
	gene_df = pd.read_csv(gene_list, sep=',')

	gene_dict = {}

	for row in gene_df.itertuples():
		
		gene_dict[row[1]] = row[1]

else:

	gene_dict = None


# have we got a singleton or a family?
filtered_ped = ped_df[ped_df['sample_id']==proband_id]

family_id = filtered_ped['family_id'].iloc[0]

affected = filtered_ped['affected'].iloc[0]

sex = filtered_ped['sex'].iloc[0]

# if any family members don't have sex set in ped
family_ped = ped_df[ped_df['family_id'] == family_id]

print(sex)

if sex == 0 or sex == '0':
	print('Sex cannot be zero - not running program.')
	f = open(output_name, 'w')
	f.write(f'Sample {proband_id} has no sex. Program not run.')
	f.close()
	exit()

has_family = True

if family_id == '0' or family_id == 0:
	
	has_family = False
	
is_affected = False
	
if affected == 2:
	
	is_affected = True

# create family object

# if we have a family and sample is affected
if has_family == True and is_affected == True:
	
	my_family = Family(family_id)
	
	my_family.read_from_ped_file(ped, family_id, proband_id)

elif has_family == False and is_affected == True:
	
	proband = FamilyMember(proband_id, proband_id, sex, is_affected)

	my_family = Family(proband_id)

	my_family.add_family_member(proband)

	my_family.set_proband(proband.get_id())

else:
	print('not running as not affected.')
	f = open(output_name, 'w')
	f.write(f'Sample {proband_id} is not affected. Program not run.')
	f.close()
	exit()


# Create a new VariantSet object
my_variant_set = VariantSet()

# Associate the my_family object with my_variant_set
my_variant_set.add_family(my_family)

# read from vcf
my_variant_set.read_variants_from_vcf(vcf,
									proband_variants_only=True,
									filter_func=passes_initial_filter,
									args=(proband_id, gene_dict, min_gq, min_dp))

# see whether we can phase comp hets by inheritance
if my_variant_set.family.proband_has_both_parents() == True:

	# Create an attribute my_variant_set.candidate_compound_het_dict where each transcript is a key the variants 
	# within that transcript are the values
	my_variant_set.get_candidate_compound_hets()

	# As we have both parents we can filter the compound hets.
	my_variant_set.filter_compound_hets()

	# Flatten the filtered compound hets (my_variant_set.filtered_compound_het_dict) into a a dictionary with each 
	# genuine compound het as a key.
	my_variant_set.get_filtered_compound_hets_as_dict()

	inheritance = ['uniparental_isodisomy','autosomal_reccessive','x_reccessive','de_novo', 'compound_het']
	
else:
	
	# Create an attribute my_variant_set.candidate_compound_het_dict where each transcript is a key the variants 
	# within that transcript are the values
	my_variant_set.get_candidate_compound_hets()

	# As we are pretending we do not have any parents we cannot phase the compound hets
	my_variant_set.get_unfiltered_compound_hets_as_dict()

	inheritance = ['uniparental_isodisomy', 'autosomal_dominant', 'autosomal_reccessive','x_reccessive','x_dominant','de_novo', 'compound_het', 'na']

# now apply second filtering function
my_variant_set.filter_variants(passes_final_filter_trio, args=(my_variant_set.final_compound_hets, inheritance, min_gq, min_dp, max_parental_alt_ref_ratio ))


# convert to pandas dataframe
variant_df = my_variant_set.to_df(min_parental_gq_dn= min_gq, min_parental_depth_dn=min_dp, min_parental_gq_upi=min_gq, min_parental_depth_upi= min_dp , max_parental_alt_ref_ratio_dn=max_parental_alt_ref_ratio)

# catch NTC
if variant_df.shape[0] ==0:

	f = open(output_name, 'w')
	f.write(f'Sample {proband_id} has no variants. Program not run.')
	f.close()
	exit()


# convert genotype fields
for fm in my_family.get_all_family_member_ids():

	variant_df[f'{fm}_Genotype'] = variant_df.apply(fix_genotype, axis=1, args=( f'{fm}_GT',))

# format pandas dataframe
gt_fields = []

proband_id = my_family.get_proband_id()
variant_df['#SampleId'] = proband_id
variant_df['WorklistId'] = worklist
variant_df['Variant'] = variant_df['variant_id']
variant_df['Genotype'] = variant_df[f'{proband_id}_Genotype']
variant_df['SYMBOL'] = variant_df['csq_SYMBOL']
variant_df['Feature'] = variant_df['csq_Feature']
variant_df['Consequence'] = variant_df['csq_Consequence']
variant_df['HGVSc'] = variant_df['csq_HGVSc']
variant_df['HGVSp'] = variant_df['csq_HGVSp']
variant_df['CLIN_SIG'] = variant_df['csq_CLIN_SIG']
variant_df['Existing_variation'] = variant_df['csq_Existing_variation']
variant_df['AutoPick'] = variant_df['csq_PICK']
variant_df['gnomad_popmax_af'] = variant_df['info_gnomad_popmax_af']
variant_df['splice_ai'] = variant_df['info_MaxSpliceAI']

# add sample specific columns such as genotype and depth
for fm in my_family.get_all_family_member_ids():

	# ignore proband as we did it earlier
	if fm != proband_id:

		for field in ['_Genotype', '_DP', '_GQ', '_AD']:

			gt_fields.append(fm + field)


csv_fields = ['#SampleId', 'WorklistId', 'Variant', 'Genotype', f'{proband_id}_DP', f'{proband_id}_GQ',  f'{proband_id}_AD', 'Feature', 'SYMBOL', 'worst_consequence', 'Consequence', 'inheritance_models',  'HGVSc', 'HGVSp', 'CLIN_SIG', 'Existing_variation', 'AutoPick', 'info_in_panel', 'gnomad_popmax_af', 'splice_ai'] + gt_fields

# save to file
variant_df[csv_fields].to_csv(output_name, index=False, sep='\t')










