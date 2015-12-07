#!/usr/bin/env python3

"""

Work tracked in DACC-584

"""

import argparse
import os
from pwd import getpwuid
import re
import shutil

def main():
    parser = argparse.ArgumentParser( description='Generates a shell script of Attributor runs given a range of pipeline IDs')

    ## output file to be written
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output shell file to be created' )
    parser.add_argument('-od', '--output_directory', type=str, required=True, help='Directory where output files will be written' )
    args = parser.parse_args()

    # inclusive
    pipeline_min = 10927248802
    pipeline_max = 11214274766

    # CONFIG
    project_area = '/usr/local/projects/dacc'
    owners = ['hhuot', 'cmccracken']
    config_template = '/usr/local/projects/jorvis/dacc_annotation/assign_functional_annotation.SRS017697.config'
    # assumes all child names are like: SRS147134.rapsearch2.m8.gz
    rapsearch2_base = '/local/hmp/dacc/t3/hhs/genome/microbiome/wgs/analysis/hmgi/rapsearch'
    attributor_path = '/home/hhuot/git/Attributor/assign_functional_annotation.py'
    # names are like either SRS045739.metagenemark.gff3 or SRS045763.metagenemark3.gff3
    gff_base_dir = '/local/hmp/dacc/t3/hhs/genome/microbiome/wgs/analysis/hmorf/gff3'

    # key = pipeline ID, value = SRS ID
    pipelines = dict()

    ## make sure the config template exists
    if not os.path.exists(config_template):
        raise Exception("ERROR: config template file not found: {0}".format(config_template))

    for pipeline_id in os.listdir("{0}/workflow/runtime/pipeline".format(project_area)):
        m = re.match("^\d+$", pipeline_id)
        if m:
            pipeline_id = int(pipeline_id)
            if pipeline_id > 11214274766 or pipeline_id < 10927248802:
                continue
        else:
            continue

        pipeline_path = "{0}/workflow/runtime/pipeline/{1}/pipeline.xml".format(project_area, pipeline_id)

        owner = find_owner(pipeline_path)
        if owner not in owners:
            continue
        
        pipeline_comment_path = "{0}/workflow/runtime/pipeline/{1}/pipeline.xml.comment".format(project_area, pipeline_id)
        if os.path.exists(pipeline_comment_path):
            comment = open(pipeline_comment_path).read()

            m = re.search("(SRS\d+)", comment)
            if m:
                srs_id = m.group(1)
                pipelines[pipeline_id] = srs_id
            else:
                raise Exception("Pipeline found without an SRS comment: {0}".format(pipeline_comment_path))

    print("INFO: found {0} pipelines".format(len(pipelines)))

    batch_fh = open(args.output_file, 'wt')

    for pipeline_id in pipelines:
        config_path = "{0}/{1}.config".format(args.output_directory, pipelines[pipeline_id])
        tmp_config_path = "/tmp/{0}.config".format(pipeline_id)

        ## get the source FASTA path
        split_fasta_config = "{0}/workflow/runtime/split_multifasta/{1}_default/split_multifasta.default.final.config".format(project_area, pipeline_id)
        fasta_path = None
        for line in open(split_fasta_config):
            line = line.rstrip()
            m = re.match("\$\;INPUT_FILE\$\;\s*=\s*(.+)", line)
            if m:
                fasta_path = m.group(1)

        ## get the source GFF path
        if os.path.exists("{0}/{1}.metagenemark.gff3".format(gff_base_dir, pipelines[pipeline_id])):
            gff3_path = "{0}/{1}.metagenemark.gff3".format(gff_base_dir, pipelines[pipeline_id])
        elif os.path.exists("{0}/{1}.metagenemark3.gff3".format(gff_base_dir, pipelines[pipeline_id])):
            gff3_path = "{0}/{1}.metagenemark.gff3".format(gff_base_dir, pipelines[pipeline_id])
        else:
            raise Exception("ERROR: failed to find GFF3 file for {0}".format(pipelines[pipeline_id]))

        if fasta_path is None:
            raise Exception("ERROR: failed to find FASTA path for pipeline {0}".format(pipeline_id))
        
        ## copy the config file
        shutil.copy(config_template, config_path)
        
        ## modify it with these paths
        ofh = open(tmp_config_path, 'wt')
        last_label = None

        for line in open(config_path):
            line = line.rstrip()
            m = re.match("   - label: (.+)", line)
            if m:
                last_label = m.group(1)

            m = re.match("     path: (.+)", line)
            if m:
                if last_label == 'coding_hmm_lib__equivalog':
                    ofh.write("     path: {0}/output_repository/hmmpfam3/{1}_default/hmmpfam3.htab.list\n".format(project_area, pipeline_id))
                    continue
                elif last_label == 'rapsearch2__uniref100__trusted_full_full':
                    ofh.write("     path: {0}/{1}.rapsearch2.m8\n".format(rapsearch2_base, pipelines[pipeline_id]))
                    continue
                elif last_label == 'rapsearch2__uniref100__trusted_partial_partial':
                    ofh.write("     path: {0}/{1}.rapsearch2.m8\n".format(rapsearch2_base, pipelines[pipeline_id]))
                    continue
                elif last_label == 'coding_hmm_lib__equivalog_domain':
                    ofh.write("     path: {0}/output_repository/hmmpfam3/{1}_default/hmmpfam3.htab.list\n".format(project_area, pipeline_id))
                    continue
                elif last_label == 'tmhmm':
                    ofh.write("     path: {0}/output_repository/tmhmm/{1}_default/tmhmm.raw.list\n".format(project_area, pipeline_id))
                    continue
                elif last_label == 'lipoprotein_motif':
                    ofh.write("     path: {0}/output_repository/lipoprotein_motif/{1}_default/lipoprotein_motif.bsml.list\n".format(project_area, pipeline_id))
                    continue

            m = re.match("   polypeptide_fasta: (\S+)", line)
            if m:
                ofh.write("   polypeptide_fasta: {0}\n".format(fasta_path))
                continue

            m = re.match("   gff3: (\S+)", line)
            if m:
                ofh.write("   gff3: {0}\n".format(gff3_path))
                continue

            ofh.write("{0}\n".format(line))
        
        ofh.close()
        shutil.move(tmp_config_path, config_path)
        batch_fh.write("{0} -f gff3 -c {1} -o {2}/{3}\n".format(attributor_path, config_path, args.output_directory,
                                                        pipelines[pipeline_id]))
        


def find_owner(filename):
    return getpwuid(os.stat(filename).st_uid).pw_name

if __name__ == '__main__':
    main()







