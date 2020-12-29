#!/usr/bin/env python


import argparse
import csv
from itertools import groupby
from operator import itemgetter

import tabix
import numpy as np
import pandas as pd



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

args = parser.parse_args()


input_file = args.input_file[0]
bed = args.bed[0]
sample_id = args.sample_id[0]
output_name = args.output_name[0]
min_depths = args.min_depths

bed_list = []
with open(bed) as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	for row in spamreader:
		bed_list.append(row)


tb = tabix.open(input_file)

def merge_gaps(gaps_list):
	
	gap_list = []
	
	for k, g in groupby(enumerate(gaps_list), lambda ix : ix[0] - ix[1]):
		gaps = list(map(itemgetter(1), g))
		gap_list.append([gaps[0], gaps[-1]])
		
	return gap_list

def get_region_stats(tabix_query, min_depths, bed_region_length, annotate=True):
	
	mean_cov_array = []
	pct_gtr_dict = {}
	
	region_length = 0
	
	gaps_dict = {}
	
	for min_depth in min_depths:
		
		gaps_dict[min_depth] = []
	
	for base in query:
										
		position = base[1]
		
		region_length = region_length +1
		
		depth = int(base[2])

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
			
		except KeyError:
			
			pct_gtr_dict[min_depth] = 0
					
	mean_coverage = np.mean(mean_cov_array)
	
	if region_length != bed_region_length +1:
		
		error = 'no results for this region'
		
	else:
		
		error = 'none'
	
	results_dict = {'mean_coverage': mean_coverage, 'pct_gtr': pct_gtr_dict, 'gaps': gaps_dict, 'region_length': region_length, 'error': error}  
	
	return results_dict


bed_region_dict = {}


for region in bed_list:
	
	chrom = region[0]
	start = int(region[1])
	end = int(region[2])

	try:

		anno = region[3]

	except:

		anno = 'none'

	region_key = f'{chrom}:{start}-{end}'
	
	query = tb.query(chrom, start-1, end)
	
	bed_region_length = end - start
	
	bed_region_length = bed_region_length -1
	
	region_stats = get_region_stats(query, min_depths, bed_region_length)

	region_stats['chrom'] = chrom
	region_stats['start'] = start
	region_stats['end'] = end
	region_stats['custom_annotation'] = anno

	if region_key not in bed_region_dict:

		bed_region_dict[region_key] = region_stats

df_row_list = []

columns = ['mean_coverage',
              'region_length',
              'error',
              'chrom',
              'start',
              'end',
              'gaps']

for cov in min_depths:
        
    columns.append(f'pct_gtr_{cov}x')


for key in bed_region_dict:
    
    mean_coverage = bed_region_dict[key]['mean_coverage']
    region_length = bed_region_dict[key]['region_length']
    error = bed_region_dict[key]['error']
    chrom = bed_region_dict[key]['chrom']
    start = bed_region_dict[key]['start']
    end = bed_region_dict[key]['end']
    gaps = bed_region_dict[key]['gaps']
    custom_annotation = bed_region_dict[key]['custom_annotation']
    
    new_row = [mean_coverage,region_length, error, chrom, start, end, gaps, custom_annotation]
    
    for cov in min_depths:
                
        new_row.append(bed_region_dict[key]['pct_gtr'][cov])
        
    df_row_list.append(new_row)


df = pd.DataFrame(df_row_list, columns=columns)

df['sample_id'] = sample_id

df.to_csv(output_name, index=False)