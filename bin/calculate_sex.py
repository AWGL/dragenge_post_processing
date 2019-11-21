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

df = pd.read_csv(args.file[0], names =['desc', 'summary', 'name', 'value1', 'value_pct'])

filtered_df = df[(df['desc'] == 'MAPPING/ALIGNING SUMMARY') & (df['name'] =='XAvgCov/YAvgCov ratio over target region')]

x_y_ratio = filtered_df.iloc[0]['value1']

if x_y_ratio < male_theshold:

	print('sample_name,calculated_sex')
	print (f'{file_name},Male')

elif x_y_ratio > female_theshold:

	print('sample_name,calculated_sex')
	print (f'{file_name},Female')

else:

	print('sample_name,calculated_sex')
	print (f'{file_name},Unknown')