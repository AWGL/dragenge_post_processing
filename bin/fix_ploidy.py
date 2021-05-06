#!/usr/bin/env python

"""
Fix ploidy for old variant database import in dragen 3.7+i

Dragen 3.7 reports sex chromosomes as ploidy 1 instead of 2 which is fine except old variant database doesn't like

Also dragen 3.7 phases variants which causes occasional fail with old variant database.


"""
import argparse
from pysam import VariantFile

parser = argparse.ArgumentParser(description='Fix ploidy for old variant database import in dragen 3.7+')
parser.add_argument('--input_vcf', type=str, nargs=1, required=True,
				help='The input vcf')
parser.add_argument('--output_vcf', type=str, nargs=1, required=True,
				help='output vcf')


args = parser.parse_args()

input_vcf = args.input_vcf[0]
output_vcf = args.output_vcf[0]

bcf_in = VariantFile(input_vcf)  

new_file = open(output_vcf, 'w')


# write header
new_file.write(str(bcf_in.header))


# loop through each variant
for rec in bcf_in.fetch():
	
	# loop through each genotype
	for sample in rec.samples:
	
		sample_gt = rec.samples[sample]['GT']
		
		# unphase
		rec.samples[sample].phased = False
		
		# if the genotype is missing and ploidy 2
		if rec.samples[sample]['GT'] == (None, None):
			
			pass
		
		else:
			
			# sort genotype - e.g. 1/0 to 0/1 - old variant database is stupid and reads these wrong
			rec.samples[sample]['GT'] = sorted((rec.samples[sample]['GT']))
		
		# if ploidy 1
		if len(sample_gt) ==1:
			
			# get ploidy one genotype
			gt_val = rec.samples[sample]['GT'][0]
			
			# exception for missing ploidy 1 values
			if rec.samples[sample]['GT'] == (None,):
				
				rec.samples[sample]['GT'] = (None, None)
				
			else:
				# make ploidy 2 e.g. 1 > 1/1 and 0 to 0/0
				rec.samples[sample]['GT'] = sorted((gt_val, gt_val))
			
	new_file.write(str(rec))

new_file.close()


