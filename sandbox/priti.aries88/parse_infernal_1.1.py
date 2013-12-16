#!/usr/bin/env python3

#import biothings
import argparse
import os
import re

def interface():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-i','--input_file', type=str, required=True, help='Path to infernal 1.1 cmscan output file.')
    parser.add_argument('-r','--rfam_file', type=str, required=False, help='Path to Rfam database file.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser = parser.parse_args()
    return parser


def read_rfam(args) :
	argsdict = vars(args)
	
	if (argsdict['rfam_file'] == None) :
		argsdict['rfam_file'] = os.path.dirname(os.path.realpath(__file__)) + "/rfam.db"
	
	rfam = {}
	for line in open(argsdict['rfam_file']) :
		if (line.startswith("Rfam")) :
			continue
		else :
			line= line.rstrip()
			val = line.split('\t')
			rfam[val[0]] = [val[2], val[3]]
	
	return rfam 
	


def parse_infernal(args, rfam) :
	fout = open(args.output_file,'w')
	fout.write ("##gff-version 3\n")
	start = 0
	c = 0
	repeat_count = 0
	for line in open(args.input_file) :

		if (line.startswith('Query:')) :
			asm_id = line[6:].split('[')[0].strip()
			continue
		
		if (line.startswith('Hit scores:') ):
			start = 1
			continue
		
		if (line.startswith('Hit alignments:')) :
			start = 0
			continue

		if start:
			pattern = re.compile('^\s+\(.\)')
			if pattern.match(line) :
				
				val = re.split('\s+',line)
				strand = val[9]
				if strand == "+" :
					start_position = val[7]
					stop_position = val[8]
				else :
					start_position = val[8]
					stop_position = val[7]
					
				feature = val[6]
				score = val[4]

				if (rfam.get(feature, 1) != 1) :
					product = rfam[feature][0]
					id = rfam[feature][1]
				else :
					print ('error')
					
			
				repeat_pattern = re.compile('repeat')
			
				if ((repeat_pattern.search(product) != None) or (id == 'regulatory_region')) :
					
					repeat_count = repeat_count + 1
					rna_id = id + "." + str(repeat_count)
					fout.write (asm_id + "\t" + "Infernal_1.1" + "\t" + id + "\t" + str(start_position) + "\t" + str(stop_position) + "\t" + str(score) \
					+ "\t" + strand + "\t.\t" + "ID="+ rna_id + ";description=" + product + ";Rfam_ID=" + feature +"\n")

				
				else :
					c = c + 1
					gene_id = "gene." + str(c)
					rna_id = gene_id + "." + id + ".1"
					exon_id = rna_id + ".exon.1"
			
					fout.write (asm_id + "\t" + "Infernal_1.1" + "\t" + "gene" + "\t" + str(start_position) + "\t" + str(stop_position) + "\t" + str(score) \
					+ "\t" + strand + "\t.\t" + "ID="+ gene_id + "\n")  

					fout.write (asm_id + "\t" + "Infernal_1.1" + "\t" + id + "\t" + str(start_position) + "\t" + str(stop_position) + "\t" + str(score) \
					+ "\t" + strand + "\t.\t" + "ID="+ rna_id + ";Parent=" + gene_id + ";product_name=" + product + ";Rfam_ID=" + feature + "\n")
				
					fout.write (asm_id + "\t" + "Infernal_1.1" + "\t" + "exon" + "\t" + str(start_position) + "\t" + str(stop_position) + "\t" + str(score) \
					+ "\t" + strand + "\t.\t" + "ID="+ exon_id + ";Parent=" + rna_id + "\n")



if __name__ == '__main__':
    args = interface()
    rfam = read_rfam(args)
    parse_infernal(args, rfam)







