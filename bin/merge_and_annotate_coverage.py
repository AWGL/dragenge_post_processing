#!/usr/bin/env python

import argparse
import glob

import pyranges as pr
import pandas as pd


def map_chromosome(chrom):


	chrom_dict = {
		'1':'NC_000001.10',
		'2':'NC_000002.11',
		'3':'NC_000003.11',
		'4':'NC_000004.11',
		'5':'NC_000005.9',
		'6':'NC_000006.11',
		'7':'NC_000007.13',
		'8':'NC_000008.10',
		'9':'NC_000009.11',
		'10':'NC_000010.10',
		'11':'NC_000011.9',
		'12':'NC_000012.11',
		'13':'NC_000013.10',
		'14':'NC_000014.8',
		'15':'NC_000015.9',
		'16':'NC_000016.9',
		'17':'NC_000017.10',
		'18':'NC_000018.9',
		'19':'NC_000019.9',
		'20':'NC_000020.10',
		'21':'NC_000021.8',
		'22':'NC_000022.10',
		'X':'NC_000023.10',
		'Y':'NC_000024.9',
		'MT':'NC_012920.1'
	}

	return chrom_dict[chrom]



def get_overlapping_genome_features(chrom, pos, end, annotations_pr):

	chrom = map_chromosome(chrom)

	variant = {'Chromosome': [chrom], 'Start': [pos], 'End': [end]}

	variant_pr = pr.from_dict(variant)

	overlapping_genome_features = variant_pr.join(annotations_pr).as_df()
	
	return overlapping_genome_features


def parse_genome_features(overlapping_genome_features):

	gene_features = {'gene': None, 'pseudogene': None,'miRNA': None,'tRNA': None  }  
	transcript_features = {'transcript': None, 'mRNA': None, 'lnc_RNA': None,'primary_transcript': None, 'snoRNA': None, 'antisense_RNA': None, 'snRNA': None, 'guide_RNA': None, 'C_gene_segment': None, 'V_gene_segment': None, 'J_gene_segment': None, 'D_gene_segment': None, 'scRNA': None }

	genes = {}
	transcripts = {}
	transcript_exons = {}

	for row in overlapping_genome_features.itertuples():

		if row.Feature in gene_features:

			if row.Feature == 'miRNA':

				biotype = 'miRNA'

			elif row.Feature == 'tRNA':

				biotype = 'tRNA'

			else:

				biotype = row.gene_biotype

			genes[row.ID] = {'name': row.Name, 'start': row.Start_b, 'end': row.End_b, 'biotype': biotype}


		if row.Feature in transcript_features:


			if pd.isnull(row.Name):

				name = row.ID

			else:

				name = row.Name

			transcripts[row.ID] = {'name': row.Name, 'start': row.Start_b, 'end': row.End_b, 'gene': row.Parent}


	for row in overlapping_genome_features.itertuples():

		if row.Feature == 'exon':

			parent = row.Parent
			start = row.Start_b
			end = row.End_b

			if parent in transcripts:

				exon_n = row.ID.split('-')[-1]

				if row.transcript_id not in transcript_exons:

					gene_id = genes[transcripts[parent]['gene']]

					transcript_exons[row.transcript_id] = {'exons': [[start, end, exon_n]], 'gene':gene_id }

				else:

					transcript_exons[row.transcript_id]['exons'].append([start, end, exon_n])

			elif parent in genes:

				if parent not in transcript_exons:

					gene_id = genes[parent]

					transcript_exons[parent] = {'exons': [[start, end, 'na']], 'gene':gene_id }

				else:

					transcript_exons[parent]['exons'].append([start, end, 'na'])  


			else:

				print('oh no')

	parsed_overlapping_genome_features = {'transcripts': transcripts, 'genes': genes, 'transcript_exons': transcript_exons}
	
	return parsed_overlapping_genome_features


files_to_merge = 