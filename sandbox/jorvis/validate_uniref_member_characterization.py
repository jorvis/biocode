#!/usr/bin/env python3

"""
Example XML entry:

<entry id="UniRef100_Q6GZX2" updated="2017-10-25">
<name>Cluster: Uncharacterized protein 3R</name>
<property type="member count" value="3"/>
<property type="common taxon" value="Ranavirus"/>
<property type="common taxon ID" value="10492"/>
<representativeMember>
  <dbReference type="UniProtKB ID" id="003R_FRG3G">
    <property type="UniProtKB accession" value="Q6GZX2"/>
    <property type="UniParc ID" value="UPI00003B0FD6"/>
  <sequence length="438" checksum="075C8FA17B3C5C56">
MARPLLGKTSSVRRRLESLSACSIFFFLRKFCQKMASLVFLNSPVYQMSNILLTERRQVD
RAMGGSDDDGVMVVALSPSDFKTVLGSALLAVERDMVHVVPKYLQTPGILHDMLVLLTPI
FGEALSVDMSGATDVMVQQIATAGFVDVDPLHSSVSWKDNVSCPVALLAVSNAVRTMMGQ
PCQVTLIIDVGTQNILRDLVNLPVEMSGDLQVMAYTKDPLGKVPAVGVSVFDSGSVQKGD
AHSVGAPDGLVSFHTHPVSSAVELNYHAGWPSNVDMSSLLTMKNLMHVVVAEEGLWTMAR
TLSMQRLTKVLTDAEKDVMRAAAFNLFLPLNELRVMGTKDSNNKSLKTYFEVFETFTIGA
LMKHSGVTPTAFVDRRWLDNTIYHMGFIPWGRDMRFVVEYDLDGTNPFLNTVPTLMSVKR
KAKIQEMFDNMVSRMVTS
  </sequence>
</representativeMember>
<member>
  <dbReference type="UniProtKB ID" id="A0A223PJ07_9VIRU">
    <property type="UniProtKB accession" value="A0A223PJ07"/>
    <property type="UniParc ID" value="UPI00003B0FD6"/>
    <property type="protein name" value="Uncharacterized protein"/>
    <property type="NCBI taxonomy" value="2026840"/>
    <property type="source organism" value="Rana catesbeiana virus 2"/>
    <property type="length" value="438"/>
  </dbReference>
</member>
<member>
  <dbReference type="UniProtKB ID" id="W8SPC6_FRG3V">
    <property type="UniProtKB accession" value="W8SPC6"/>
    <property type="UniParc ID" value="UPI00003B0FD6"/>
    <property type="protein name" value="IIV6 orf229L like protein"/>
    <property type="NCBI taxonomy" value="10493"/>
    <property type="source organism" value="Frog virus 3 (FV-3)"/>
    <property type="length" value="438"/>
  </dbReference>
</member>
</entry>

Actions:

For each representative member, check if it's my trembl sqlite3 at all.  If so,
check if it's characterized.

Then check if each member both exists in the sqlite3 db and that its characterization
matches the representative member.
"""


import argparse
import gzip
import os
import re
import sqlite3

def main():
    parser = argparse.ArgumentParser( description='Just testing for now')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the uniref100.xml file.' )
    parser.add_argument('-t', '--trembl_db', type=str, required=True, help='Path to an input SQLite3 db for trembl' )
    parser.add_argument('-s', '--sprot_db', type=str, required=True, help='Path to an input SQLite3 db for sprot' )
    args = parser.parse_args()

    trembl_conn = sqlite3.connect(args.trembl_db)
    trembl_curs = trembl_conn.cursor()

    sprot_conn = sqlite3.connect(args.sprot_db)
    sprot_curs = sprot_conn.cursor()

    records_processed = 0
    is_characterized = 0
    is_compressed = False
    line_num = 0

    characterized_go_codes = ['EXP', 'IDA', 'IMP', 'IGI', 'IPI', 'IEP']

    if args.input.endswith('.gz'):
        ifh = gzip.open(args.input, 'rb')
        is_compressed = True
    else:
        ifh = open(args.input)
    
    rep_id = None
    member_ids = list()
    in_rep_section = False
    in_member_section = False

    print("INFO: Parsing XML file ...")
    for line in ifh:
        line_num += 1

        if is_compressed:
            line = line.decode()
        
        if line.startswith('<entry'):
            records_processed += 1

            if records_processed % 1000 == 0:
                print("{0} records processed ...".format(records_processed))

        elif line.startswith('<representativeMember'):
            in_rep_section = True

        elif line.startswith('</representativeMember'):
            in_rep_section = False

        elif line.startswith('<member>'):
            in_member_section = True
        
        elif line.startswith('</member>'):
            in_member_section = False

        elif line.startswith('  <dbReference type="UniProtKB ID" id="'):
            m = re.match('  <dbReference type="UniProtKB ID" id="(.+)"', line)
            if m:
                if in_member_section:
                    member_ids.append(m.group(1))
                elif in_rep_section:
                    rep_id = m.group(1)
                else:
                    raise Exception("Found an ID reference line outside of member or rep sections on line {0}".format(line_num))

        elif line.startswith('</entry'):
            # did we find a rep_id?
            if rep_id is None:
                raise Exception("Logic error.  On line {0} failed to find a rep ID before a closing entry".format(line_num))

            # do stuff
            rep_char = get_characterized_status(trembl_curs, sprot_curs, rep_id)
            # print("Rep ID: {0} ({1})\n  Members:".format(rep_id, rep_char))
            
            for member_id in member_ids:
                member_char = get_characterized_status(trembl_curs, sprot_curs, member_id)
                #print("\t{0} ({1})".format(member_id, member_char))

                #if member_char != rep_char:
                if member_char == 0 and rep_char == 1:
                    print("Characterization disagreement.  Repchar: Rep: {0} ({1}), Member: {2} ({3})".format(rep_id, rep_char, member_id, member_char))
                elif member_char == 1 and rep_char == 0:
                    print("Characterization disagreement.  Memberchar: Rep: {0} ({1}), Member: {2} ({3})".format(rep_id, rep_char, member_id, member_char))
            # now reset
            rep_id = None
            member_ids = list()
        
    trembl_curs.close()
    print("INFO: Complete.")
    

def get_characterized_status(trembl_cursor, sprot_cursor, id):
    # returns 1, 0 or None (if the accession wasn't in the database at all)
    dsl = "SELECT is_characterized FROM entry_acc WHERE id = ?"
    for row in trembl_cursor.execute(dsl, (id,)):
        return row[0]

    # if we got this far, it wasn't in the trembl db
    dsl = "SELECT is_characterized FROM entry_acc WHERE id = ?"
    for row in sprot_cursor.execute(dsl, (id,)):
        return row[0]

    return None
    

if __name__ == '__main__':
    main()







