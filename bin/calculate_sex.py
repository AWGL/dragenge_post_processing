#!/usr/bin/env python

import argparse
import pandas as pd
from pathlib import Path


parser = argparse.ArgumentParser(description='Calculate the sex from the dragen alignment metrics file')
parser.add_argument('--file', type=str, nargs=1, required=True,
				help='The input metrics file')

args = parser.parse_args()

file_name = Path(args.file[0]).stem

female_theshold = 150
male_theshold = 50

df = pd.read_csv(args.file[0], sep='\t', header=None, names=['chr', 'pos', 'cov1', 'cov2', 'cov3'])

df_grouped = df.groupby('chr').mean()
df_grouped['chr'] = df_grouped.index

x_cov = df_grouped[df_grouped['chr']=='X'].iloc[0]['cov1']
y_cov = df_grouped[df_grouped['chr']=='Y'].iloc[0]['cov1']

x_y_ratio = x_cov/y_cov

if x_y_ratio < male_theshold:

	print('sample_name,calculated_sex')
	print (f'{file_name},Male')

elif x_y_ratio > female_theshold:

	print('sample_name,calculated_sex')
	print (f'{file_name},Female')

else:

	print('sample_name,calculated_sex')
	print (f'{file_name},Unknown')