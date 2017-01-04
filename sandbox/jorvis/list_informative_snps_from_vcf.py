#!/usr/bin/env python3.2

'''
Reads the SNPs from a reference VCF file and and then any number of
subject VCF files.  Every SNP in the reference is compared with those
in the same position in the subjects.  It is then marked as
informative if any of the following are true:

  - The subject doesn't have a SNP in the same position
  - The subject has a SNP, but a different allele, indicating no consensus among all subjects

Currently only generates counts

Skips entries in a VCF file which don't pass the following conditions:

    - The 5th column must not be '.'
    - The 7th column must be 'PASS'

Comment lines are retained.

Will NOT currently work with merged sample files.

Warning: currently only checks exact allele field matches, so the following would all
be considered different and flagged as 'informative'

   ref:T vs. qry:AT
   ref:CATATATAT,C vs. qry:CATAT,CATATAT
   ref:T,TTA vs. qry:T
'''


import argparse
import gzip


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    parser.add_argument('-r', '--reference_vcf_file', type=str, required=True, help='Input VCF file' )
    parser.add_argument('-q', '--query_vcf_files', type=str, required=True, help='Comma-separated list of VCF files to which the reference should be compared.')
    args = parser.parse_args()

    ref_mols = parse_snps_from_file(args.reference_vcf_file)

    qry_files = args.query_vcf_files.split(",")
    
    for path in qry_files:
        print("INFO: processing query VCF: {0}".format(path))
        qry_mols = parse_snps_from_file(path)
        compare_ref_with_qry(ref_mols, qry_mols)

    total_ref_snp_count   = 0
    informative_snp_count = 0

    for mol_id in ref_mols:
        for coord in ref_mols[mol_id]:
            total_ref_snp_count += 1

            if ref_mols[mol_id][coord]['inf'] is True:
                informative_snp_count += 1
    
    print("RESULTS: Total SNP count in reference: {0}".format(total_ref_snp_count) )
    print("RESULTS: Informative SNP count: {0}".format(informative_snp_count) )

    
def compare_ref_with_qry( refs, qrys ):
    for ref_mol_id in refs:
        ## problem if this doesn't exist at all
        if ref_mol_id not in qrys:
            raise Exception("ERROR: reference molecule ID ({0}) not found in qry VCF".format(ref_mol_id))
        
        ## check each SNP position
        for ref_snp_coord in refs[ref_mol_id]:
            ref_snp = refs[ref_mol_id][ref_snp_coord]
            ## informative if this SNP position isn't in the qry set at all
            if ref_snp_coord not in qrys[ref_mol_id]:
                ref_snp['inf'] = True
                #print("DEBUG: {0}:{1} SNP informative because not found in qry".format(ref_mol_id, ref_snp_coord) )

            ## informative if both have a SNP in the same coordinate but with different alleles
            elif not have_shared_allelic_variation(ref_snp['alt'], qrys[ref_mol_id][ref_snp_coord]['alt']):
                ref_snp['inf'] = True
                #print("DEBUG: {0}:{1} SNP informative because variant allele in qry (ref:{2} vs. qry:{3})".format( \
                #        ref_mol_id, ref_snp_coord, ref_snp['alt'], qrys[ref_mol_id][ref_snp_coord]['alt']) )
                

def have_shared_allelic_variation( ref_string, qry_string ):
    """
    Because the VCF files have ellele columns with comma-separated values, this function
    takes two of them and looks for any single shared allele among them, returning True
    if one is found.
    """
    refs = ref_string.split(",")
    qrys = qry_string.split(",")

    for ref in refs:
        for qry in qrys:
            if ref == qry:
                #print("FOUND: ref:{0} == qry:{1}".format( ref_string, qry_string ) )
                return True

    #print("NOT FOUND: ref:{0} == qry:{1}".format( ref_string, qry_string ) )
    return False
        

def parse_snps_from_file( file_path ):
    molecules = dict()
    snp_count = 0

    fh = None
    if file_path.endswith(".gz"):
        fh = gzip.GzipFile(file_path, "r")
    else:
        fh = open(file_path)

    total_seq_length = 0
    sample_id = None
        
    for line in fh:
        line = line.rstrip()
        cols = line.split("\t")

        if line.startswith("#"):
            continue
        
        cols = line.split("\t")

        mol = cols[0]
        pos = cols[1]
        alt = cols[4]
        
        if len(cols) < 7 or alt == '.' or cols[6] != 'PASS':
            continue

        ## store this SNP
        if mol not in molecules:
            molecules[mol] = dict()

        molecules[mol][pos] = { 'alt': alt, 'inf': False }
        snp_count += 1

    print("INFO: parsed {0} SNPs in file: {1}".format(snp_count, file_path) )

    return molecules



if __name__ == '__main__':
    main()












    
