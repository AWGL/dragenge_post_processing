#!/usr/bin/env python

"""
Fix ploidy for old variant database import in dragen 3.7+
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

new_file.write(str(bcf_in.header))

for rec in bcf_in.fetch():
	
	for sample in rec.samples:
				
		sample_gt = rec.samples[sample]['GT']
		
		# if ploidy is 1 e.g. X and Y then may ploidy 2 using whatever genotype is availible
		if len(sample_gt) ==1:
			
			gt_val = rec.samples[sample]['GT'][0]
			
			rec.samples[sample]['GT'] = (gt_val, gt_val)
			
	new_file.write(str(rec))

new_file.close()


