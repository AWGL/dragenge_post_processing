#!/usr/bin/env python

"""
Script for calculating the coverage over various regions given a Per base Depth of Coverage file

Produces

1) _region_coverage_data.csv : Raw Data file
2) _region_coverage_summary.csv: Summary of coverage over each region
3) _gene_coverage_summary.csv: Summary of coverage over each gene
4) _gaps_20x.csv: Gaps at a certian depth

Input:

min_depths: int: list of thresholds for calculating metrics e.g. pct >x
input_file: str: path to the gatk depth of coverage file (tabixed and bgzipped)
bed: str: bed file of regions to report coverage over
sample_id: str : the sample id
output_name: str: root name of the output
calculate_gene_stats_gaps: bool: whether to output gaps etc files
"""

import argparse
import csv
from itertools import groupby
from operator import itemgetter

import tabix
import numpy as np
import pandas as pd

def get_gene(df):
	"""
	Extract the custom annotation from the gene bed file
	"""
	
	gene = df['custom_annotation']
	
	if ':' in gene:
		
		return gene.split(':')[0]
	
	elif '.' in gene:
		
		return gene.split('.')[0]
	
	else:
		
		return gene

def normalise_value_on_region_length(df, gene_length_dict, value):
	"""
	For certain values such as mean coverage we need to normalise \
	as regions have different lengths. For example one exon might be very short but have \
	high coverage throwing off mean coverage for a gene.
	
	gene_length_dict: dict of gene lengths i.e. sum of all length of all exons in gene.
	value: which value to normalise

	"""
	
	value = df[value]
	region_length = df['region_length']
	gene = df['gene']
	
	
	gene_length = gene_length_dict[gene]
	
	return value * (region_length / gene_length)



def merge_gaps(gaps_list):
	"""
	Given a list of base positions which are below a threshold e.g.

	[1,2,3,10,11,12]

	merge these into contigous regions e.g.

	[[1-3], [10-12]]

	"""
	
	gap_list = []
	
	for k, g in groupby(enumerate(gaps_list), lambda ix : ix[0] - ix[1]):
		gaps = list(map(itemgetter(1), g))
		gap_list.append([gaps[0], gaps[-1]])
		
	return gap_list

def get_region_stats(tabix_query, min_depths, bed_region_length):
	"""
	given a tabix query of the depth of coverage file calculate various stats on it.

	min_depths: list: list of depth thresholds e.g. [20,160]
	bed_region_length: int: the length of the bed region of interest. 



	"""
	# init some variables
	mean_cov_array = []
	pct_gtr_dict = {}
	region_length = 0
	gaps_dict = {}


	# create initial dictionary keys
	for min_depth in min_depths:
		
		gaps_dict[min_depth] = []
	
	# now loop through the tabix query
	for base in tabix_query:
		
		# get position and depth from query								
		position = base[1]
		depth = int(base[2])
		
		# increment length
		region_length = region_length +1

		# add to mean coverage array (for calculation of mean depth over region)
		mean_cov_array.append(depth)
		
		# get pct > than for min depths
		for min_depth in min_depths:
			
			if depth >= min_depth:
		
				if min_depth not in pct_gtr_dict:
				
					pct_gtr_dict[min_depth] = 1
					
				else:
					
					pct_gtr_dict[min_depth] = pct_gtr_dict[min_depth] +1
					
			if depth < min_depth:
					
				gaps_dict[min_depth].append(int(position))
					
	
	for min_depth in min_depths:
		
		gaps_dict[min_depth] = merge_gaps(gaps_dict[min_depth])

		try:
			pct_gtr_dict[min_depth] = pct_gtr_dict[min_depth] / region_length
		
		# in case region length is 0
		except KeyError:
			
			pct_gtr_dict[min_depth] = 0
	
	# get mean using numpy
	mean_coverage = np.mean(mean_cov_array)
	
	# qc check to ensure that errors don't occur
	if region_length == bed_region_length +1:
		
		error = 'none'

	# qc check to ensure that errors don't occur
	elif region_length == bed_region_length +2:
		
		error = 'none'
		
	else:
		print (region_length,bed_region_length )
		error = 'region length error'

	# create results dict for each region
	results_dict = {'mean_coverage': mean_coverage, 'pct_gtr': pct_gtr_dict, 'gaps': gaps_dict, 'region_length': region_length, 'error': error}  
	
	return results_dict


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Calculate statistics across coverage reporting regions.')

	parser.add_argument('min_depths', metavar='N', type=int, nargs='+',
						help='Coverage cutoffs')
	parser.add_argument('--input_file', type=str, nargs=1, required=True,
					help='The per base coverage file')
	parser.add_argument('--bed', type=str, nargs=1, required=True,
					help='Bed file of regions to report coverage over.')
	parser.add_argument('--sample_id', type=str, nargs=1, required=True,
					help='sample_id for adding to csv.')
	parser.add_argument('--output_name', type=str, nargs=1, required=True,
					help='Where to store and name the output file.')
	parser.add_argument('--calculate_gene_stats_gaps', type=bool, nargs=1, required=True,
					help='Create files with gene level gaps and stats.')

	args = parser.parse_args()


	input_file = args.input_file[0]
	bed = args.bed[0]
	sample_id = args.sample_id[0]
	output_name = args.output_name[0]
	calculate_gene_stats_gaps = args.calculate_gene_stats_gaps[0]
	min_depths = args.min_depths

	# read in bed
	bed_list = []
	with open(bed) as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			bed_list.append(row)

	# get tabix object from input
	tb = tabix.open(input_file)

	bed_region_dict = {}

	# for each region get the stats
	for region in bed_list:
		
		chrom = region[0]
		start = int(region[1])
		end = int(region[2])

		try:

			anno = region[3]

		except:

			anno = 'none'

		region_key = f'{chrom}:{start}-{end}-{anno}'
		try:
			query = tb.query(chrom, start-1, end)
		except:
			raise Exception (f'Could not process region {region_key}')
		
		# get bed region length
		bed_region_length = end - start
		bed_region_length = bed_region_length -1
		
		# get region stats
		region_stats = get_region_stats(query, min_depths, bed_region_length)

		# add additional info
		region_stats['chrom'] = chrom
		region_stats['start'] = start
		region_stats['end'] = end
		region_stats['custom_annotation'] = anno

		# add to bed_region_dict - cached to ensure
		if region_key not in bed_region_dict:

			bed_region_dict[region_key] = region_stats


	# create a dataframe
	df_row_list = []

	columns = ['mean_coverage',
				  'region_length',
				  'error',
				  'chrom',
				  'start',
				  'end',
				  'gaps',
				  'custom_annotation']

	for cov in min_depths:
			
		columns.append(f'pct_gtr_{cov}x')


	for key in bed_region_dict:
		
		mean_coverage = round(bed_region_dict[key]['mean_coverage'],4)
		region_length = bed_region_dict[key]['region_length']
		error = bed_region_dict[key]['error']
		chrom = str(bed_region_dict[key]['chrom'])
		start = bed_region_dict[key]['start']
		end = bed_region_dict[key]['end']
		gaps = bed_region_dict[key]['gaps']
		custom_annotation = bed_region_dict[key]['custom_annotation']
		
		new_row = [mean_coverage,region_length, error, chrom, start, end, gaps, custom_annotation]
		
		for cov in min_depths:
					
			new_row.append(round(bed_region_dict[key]['pct_gtr'][cov],4))
			
		df_row_list.append(new_row)


	df = pd.DataFrame(df_row_list, columns=columns)

	df['sample_id'] = sample_id

	# write raw data
	df.to_csv(f'{output_name}_region_coverage_data.csv', index=False)


	if calculate_gene_stats_gaps:

		# create region summary

		df['gene'] = df.apply(get_gene,axis=1)

		columns = ['chrom', 'start', 'end', 'mean_coverage']

		for depth in min_depths:
		
			columns.append(f'pct_gtr_{depth}x')
		
		columns.append('custom_annotation')
		columns.append('gene')


		df[columns].to_csv(f'{output_name}_region_coverage_summary.csv', index=False, sep='\t')


		# create gene summary

		# get lengths of gene roi
		gene_length_dict = df[['region_length', 'gene']].groupby('gene').sum().to_dict()['region_length']

		# normalise mean coverage by gene
		df['normalised_mean_coverage'] = df.apply(normalise_value_on_region_length, axis=1, args=(gene_length_dict, 'mean_coverage',))

		# normalise 
		for depth in min_depths:
		
			df[f'normalised_pct_gtr_{depth}x'] = df.apply(normalise_value_on_region_length, axis=1, args=(gene_length_dict, f'pct_gtr_{depth}x',))


		# get normalised mean coverage by gene
		mean_coverage_dict = df[['normalised_mean_coverage', 'gene']].groupby('gene').sum().to_dict()['normalised_mean_coverage']

		# get normalised pct >20x for each gene
		pct_gt_dict = {}

		for depth in min_depths:

			pct_gt_dict[f'pct_gtr_{depth}x'] = df[[f'normalised_pct_gtr_{depth}x', 'gene']].groupby('gene').sum().to_dict()[f'normalised_pct_gtr_{depth}x']

		# write gene summary to file
		with open(f'{output_name}_gene_coverage_summary.csv', 'w') as csvfile:
			spamwriter = csv.writer(csvfile, delimiter=',')
			
			headers = ['gene']

			for depth in min_depths:
					
				headers.append(f'pct_gtr_{depth}x')

			headers.append('mean_depth')
				
			spamwriter.writerow(headers)
			
			for key in gene_length_dict:
				
				row = [key]
				
				for depth in min_depths:
					
					row.append(round(pct_gt_dict[f'pct_gtr_{depth}x'][key],4))

				row.append(round(mean_coverage_dict[key],4))
					
					
				spamwriter.writerow(row)

		# process gaps
		gaps_dict = {}

		for depth in min_depths:
			
			gaps_dict[depth] = []


		for row in df.itertuples():
			
			chromosome = row.chrom
			start = row.start
			end = row.end
			gaps = row.gaps
			custom_annotation = row.custom_annotation
			gene = row.gene
					
			for key in gaps:
						
				for specific_gap in gaps[key]:
						
					gaps_dict[key].append([chromosome, specific_gap[0], specific_gap[1], custom_annotation, gene])


		# write gaps to file
		for key in gaps_dict:
				
			if len(gaps_dict[key]) > 0:
			
				with open(f'{output_name}_gaps_{key}x.csv', 'w') as csvfile:
					spamwriter = csv.writer(csvfile, delimiter='\t')
					
					spamwriter.writerow(['chromosome', 'gap_start', 'gap_end', 'annotation', 'gene'])


					for gap in gaps_dict[key]:

						spamwriter.writerow(gap)
			
			else:
				
				with open(f'{output_name}_gaps_{key}x.csv', 'w') as csvfile:
					
					spamwriter = csv.writer(csvfile, delimiter='\t')
					
					spamwriter.writerow(['chromosome', 'gap_start', 'gap_end', 'annotation', 'gene'])
