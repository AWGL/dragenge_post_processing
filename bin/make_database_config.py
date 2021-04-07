#!/usr/bin/env python

import argparse
import glob
import csv
from pathlib import Path

"""
Get the ped file and add relevant variables and then output csv

"""


parser = argparse.ArgumentParser(description='Create a config file for the variant database.')
parser.add_argument('--ped_file', type=str, nargs=1, required=True,
				help='file path for ped file.')
parser.add_argument('--variables', type=str, nargs=1, required=True,
				help='file path to folder containing variables files.')
parser.add_argument('--pipeline_name', type=str, nargs=1, required=True,
				help='pipeline name')
parser.add_argument('--output_name', type=str, nargs=1, required=True,
				help='output name')

args = parser.parse_args()

ped_file = args.ped_file[0]
pipeline_name = args.pipeline_name[0]
output_name = args.output_name[0]
variable_files = glob.glob(args.variables[0])

# read variables files
ped_dict = {}

for file in variable_files:

	sample_name = Path(file).stem

	f = open(file, 'r')

	variables_dict = {}

	for line in f:

		new_line = line.strip()

		if len(new_line) == 0:

			pass

		elif new_line[0] == '#':

			pass

		else:

			split_line = new_line.split('=')

			variables_dict[split_line[0]] = split_line[1]

	ped_dict[sample_name] = variables_dict

# read ped file
config_rows = []

with open(ped_file) as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	for row in spamreader:
		config_rows.append(row + [ped_dict[row[1]]['panel']] + [pipeline_name] + [ped_dict[row[1]]['seqId']] + [ped_dict[row[1]]['worklistId']])

# write config
with open(output_name, 'w') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    for row in config_rows:
    	spamwriter.writerow(row)