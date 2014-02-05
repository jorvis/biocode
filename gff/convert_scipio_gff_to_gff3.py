#!/usr/bin/env python3

"""
The default Scipio output, while technically GFF3, has no parent/child relationships between
matches and match parts, so any displays reading the GFF show very segmented matches.  This
script corrects this:

Example input:

    ##query	jgi|Conco1|10118|gm1.8358_g 1 477
    # gap 1 41
    jcf7180000788868	Scipio	protein_match	29004	29234	0.644	+	0	ID=4604;Query=jgi|Conco1|10118|gm1.8358_g 42 122;Mismatches=54 56 57 60 63 66 77 78 82 83 89 104 115 118;Gap=M74 I1 M2 I3 M1
    jcf7180000788868	Scipio	protein_match	29288	29420	0.644	+	0	ID=4604;Query=jgi|Conco1|10118|gm1.8358_g 123 166;Mismatches=144 146 152 160 163 165
    jcf7180000788868	Scipio	protein_match	29475	29784	0.644	+	2	ID=4604;Query=jgi|Conco1|10118|gm1.8358_g 167 270;Mismatches=177 180 181 184 191 192 224 246 256 268
    jcf7180000788868	Scipio	protein_match	29839	30075	0.644	+	1	ID=4604;Query=jgi|Conco1|10118|gm1.8358_g 271 349;Mismatches=281 284 310 314 315 331
    jcf7180000788868	Scipio	protein_match	30136	30182	0.644	+	1	ID=4604;Query=jgi|Conco1|10118|gm1.8358_g 350 364;Mismatches=359
    jcf7180000788868	Scipio	protein_match	30312	30616	0.644	+	2	ID=4604;Query=jgi|Conco1|10118|gm1.8358_g 365 466;Mismatches=365 377 383 384 385 388 403 406 407 410 414 416 417 419 420 426 459 461 463 465
    # gap 467 477
    # protein sequence = [.........................................GGLSDEDRIFTNxYxxHDxKLxGAxKRGDWYKTKExx
    # LKGxxWIINExKESGLRGRGGAGFPxGLKWSFMNKPx.Dx...RPRYLVINADEGEPGTCKDREIxRxDPHKLxEGCLVAGxAMxAxAAYIYIRGEFY
    # xEAxxLQxAINEAYxxGLIGKNACGSGYDFDVYIHRGAGAYICGEETxLIESIEGKQGKPRLKPPFPADxGLFGCPTTVxNVETVAVAPTIxRRGGSW
    # FAGFGRxRNxGTKLFCISGHVNNPCTVEEEMSIPLxELIxxHCGGVRGGWDNLLGIxPGGSSVPVLPKSICDDVLMDFDALRDAxSGLGTxAVIVMDK
    # STDIxRAISRxxxFYxHESCGQCTPCREGTxWLxxMMxRFExGxxHxxEIDQIxELTKQIEGHTICALGDAAAWPVQGLIRHFRPExExRxAxF....
    # .......]
    #

Example output:
    jcf7180000788868        Scipio  protein_match   29004   30616   .       +       .       ID=Scipio.match.4
    jcf7180000788868        Scipio  match_part      29004   29234   0.644   +       0       ID=4604;Parent=Scipio.match.4;Target=jgi|Conco1|10118|gm1.8358_g 42 122;Mismatches=54 56 57 60 63 66 77 78 82 83 89 104 115 118;Gap=M74 I1 M2 I3 M1
    jcf7180000788868        Scipio  match_part      29288   29420   0.644   +       0       ID=4604;Parent=Scipio.match.4;Target=jgi|Conco1|10118|gm1.8358_g 123 166;Mismatches=144 146 152 160 163 165
    jcf7180000788868        Scipio  match_part      29475   29784   0.644   +       2       ID=4604;Parent=Scipio.match.4;Target=jgi|Conco1|10118|gm1.8358_g 167 270;Mismatches=177 180 181 184 191 192 224 246 256 268
    jcf7180000788868        Scipio  match_part      29839   30075   0.644   +       1       ID=4604;Parent=Scipio.match.4;Target=jgi|Conco1|10118|gm1.8358_g 271 349;Mismatches=281 284 310 314 315 331
    jcf7180000788868        Scipio  match_part      30136   30182   0.644   +       1       ID=4604;Parent=Scipio.match.4;Target=jgi|Conco1|10118|gm1.8358_g 350 364;Mismatches=359
    jcf7180000788868        Scipio  match_part      30312   30616   0.644   +       2       ID=4604;Parent=Scipio.match.4;Target=jgi|Conco1|10118|gm1.8358_g 365 466;Mismatches=365 377 383 384 385 388 403 406 407 410 414 416 417 419 420 426 459 461 463 465

Future development possibilities:
- Add an output mode for exporting alignments as gene models or match/match_parts

"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser( description='Converts Scipio GFF output to GFF3 with grouped match_parts')

    # output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to parse' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    
    args = parser.parse_args()
    match_part_rows = list()
    match_atts = dict()
    next_match_num = 1

    ofh = open(args.output_file, 'w')

    for line in open(args.input_file, 'r'):
        if line.startswith('##query'):
            if len(match_part_rows) > 0:
                export_match( ofh, match_part_rows, match_atts )
                match_part_rows = list()
                match_atts = dict()
                continue
        elif line.startswith('#'):
            continue
        
        cols = line.split("\t")

        if len(cols) != 9:
            continue

        # if this is the first match part, populate the overall match attributes we can.
        if len(match_atts) == 0:
            match_atts['molecule']  = cols[0]
            match_atts['algorithm'] = cols[1]
            match_atts['strand'] = cols[6]
            match_atts['min_coord'] = cols[3]
            match_atts['max_coord'] = cols[4]
            match_atts['id'] = "{0}.match.{1}".format( match_atts['algorithm'], next_match_num )
            next_match_num += 1

        if cols[3] < match_atts['min_coord']:
            match_atts['min_coord'] = cols[3]
        
        if cols[4] > match_atts['max_coord']:
            match_atts['max_coord'] = cols[4]

        cols[2] = 'match_part'

        # set a Parent feature
        cols[8] = cols[8].replace(';Query', ';Parent=' + match_atts['id'] + ';Target')

        line = "\t".join(cols)
        
        match_part_rows.append( line )
    

    # make sure to do the last one
    if len(match_part_rows) > 0:
        export_match( ofh, match_part_rows, match_atts )

    


def export_match( out, match_rows, atts ):
    
    # write the match_feature
    data_column = "ID={0}".format(atts['id'])
    out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
            atts['molecule'], atts['algorithm'], 'protein_match', atts['min_coord'], atts['max_coord'], \
            '.', atts['strand'], '.', data_column
            ) )

    # write the match_parts
    for line in match_rows:
        #data_column = "ID={0};Parent={1};Target={2} {3} {4}".format(mRNA_id, gene_id, hit_id, hit_min, hit_max)
        out.write(line)
        #out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( row ) )



    return True

if __name__ == '__main__':
    main()





















