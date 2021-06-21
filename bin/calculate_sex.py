#!/usr/bin/env python

"""
Calculate sex using X and Y coverage.
"""

import argparse
from pathlib import Path

import pandas as pd


parser = argparse.ArgumentParser(description='Calculate the sex from the dragen alignment metrics file')
parser.add_argument('--file', type=str, nargs=1, required=True,
				help='The input metrics file')
parser.add_argument('--female_theshold', type=int, nargs=1, required=True,
				help='female_theshold')
parser.add_argument('--male_theshold', type=int, nargs=1, required=True,
				help='male_theshold')
parser.add_argument('--genome_build', type=str, nargs=1, required=True,
				help='genome build - either GRCh37 or GRCh38')

args = parser.parse_args()

file_name = Path(args.file[0]).stem

female_theshold = args.female_theshold[0]
male_theshold = args.male_theshold[0]
genome_build = args.genome_build[0]

print (genome_build)

df = pd.read_csv(args.file[0], sep='\t', header=None, names=['chr', 'pos', 'cov1', 'cov2', 'cov3'], dtype={'chr': str} )

df_grouped = df.groupby('chr').mean()
df_grouped['chr'] = df_grouped.index

if genome_build == 'GRCh37':

	x_cov = df_grouped[df_grouped['chr']=='X'].iloc[0]['cov1']
	y_cov = df_grouped[df_grouped['chr']=='Y'].iloc[0]['cov1']

elif genome_build == 'GRCh38':

	x_cov = df_grouped[df_grouped['chr']=='chrX'].iloc[0]['cov1']
	y_cov = df_grouped[df_grouped['chr']=='chrY'].iloc[0]['cov1']

else:

	raise Exception('Invalid genome_build')

x_y_ratio = x_cov / y_cov

if x_y_ratio < male_theshold:

	print('sample_name,calculated_sex')
	print (f'{file_name},Male')

elif x_y_ratio > female_theshold:

	print('sample_name,calculated_sex')
	print (f'{file_name},Female')

else:

	print('sample_name,calculated_sex')
	print (f'{file_name},Unknown')
