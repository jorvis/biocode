#!/usr/local/packages/Python-3.3.2/bin/python3

"""

Started at 3:30 PM:
$ /usr/local/packages/bowtie2-2.2.4/bowtie2 -x /usr/local/projects/dacc/jorvis/read_to_metaref_alignments/mumi.20150105.genomic -1 /usr/local/scratch/jorvis/dacc/read_to_metaref_seed_alignments/phase2/SRS143214/SRS143214.denovo_duplicates_marked.trimmed.1.fastq -2 /usr/local/scratch/jorvis/dacc/read_to_metaref_seed_alignments/phase2/SRS143214/SRS143214.denovo_duplicates_marked.trimmed.2.fastq -U /usr/local/scratch/jorvis/dacc/read_to_metaref_seed_alignments/phase2/SRS143214/SRS143214.denovo_duplicates_marked.trimmed.singleton.fastq -S /usr/local/scratch/jorvis/dacc/read_to_metaref_seed_alignments/phase2/SRS143214/SRS143214.vs_metaref_seeds.bowtie2.sam -a >& run.out

509710 reads; of these:
  509710 (100.00%) were paired; of these:
    325656 (63.89%) aligned concordantly 0 times
    72035 (14.13%) aligned concordantly exactly 1 time
    112019 (21.98%) aligned concordantly >1 times
    ----
    325656 pairs aligned concordantly 0 times; of these:
      511 (0.16%) aligned discordantly 1 time
    ----
    325145 pairs aligned 0 times concordantly or discordantly; of these:
      650290 mates make up the pairs; of these:
        630220 (96.91%) aligned 0 times
        12267 (1.89%) aligned exactly 1 time
        7803 (1.20%) aligned >1 times

509710 reads?  Files were removed from scratch during runtime

Test execution (4.7GB tarball):
./generate_read_to_metaref_seed_alignment.py -r /local/projects-t2/dacc/dmz_Illumina/PHASEII/stool/SRS148424.tar.bz2 -s SRS148424

Test SGE execution:
qsub -P owhite-dacc-irc -q threaded.q -v PATH -pe thread 8 -l mem_free=25G -wd /usr/local/scratch/jorvis/dacc/read_to_metaref_seed_alignments -o SRS140663.process.log -e SRS140663.process.stderr -b y /usr/local/projects/dacc/bin/generate_read_to_metaref_seed_alignment.py -r /local/projects-t2/dacc/dmz_Illumina/PHASEII/anterior_nares/SRS140663.tar.bz2 -s SRS140663 -c 8



"""

import argparse
import datetime
import os
import subprocess


def main():
    READ_BASE_DIR = "/local/projects-t2/dacc/dmz_Illumina"
    #WORK_BASE_DIR = "/usr/local/scratch/jorvis/dacc/read_to_metaref_seed_alignments/partition20"
    #WORK_BASE_DIR = "/local/scratch2/dacc/read_to_metaref_seed_alignments/partition20"
    WORK_BASE_DIR = "/local/hmp/dacc/restore/read_to_metaref_seed_alignments/partition20"
    BOWTIE_PATH = "/usr/local/packages/bowtie2-2.2.4/bowtie2"
    COMPLETION_BASE = "/local/hmp/dacc/restore/read_to_metaref_seed_alignments/complete/microbes/"

    # any reads >= this percentage of Ns will be removed
    FASTQ_FILTERING_N_PCT_CUTOFF = 80

    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-r', '--read_file', type=str, required=True, help='Read file (in tar.bz2 format)' )
    parser.add_argument('-s', '--sample', type=str, required=True, help='sample_base_name, like SRS016516' )
    parser.add_argument('-c', '--cpu_cores', type=int, required=False, default=8, help='CPU count to use for bowtie step' )
    args = parser.parse_args()

    ## decompress it
    # tar -xjf posterior_fornix/SRS016516.tar.bz2 -C /tmp/
    # ls /tmp/SRS016516/
    cmd = "tar -xjf {0} -C {1}/".format(args.read_file, WORK_BASE_DIR)
    run_command(cmd)
    
    ## creates three files:
    # SRS016516.denovo_duplicates_marked.trimmed.1.fastq
    # SRS016516.denovo_duplicates_marked.trimmed.2.fastq
    # SRS016516.denovo_duplicates_marked.trimmed.singleton.fastq

    scratch_dir = "{0}/{1}".format(WORK_BASE_DIR, args.sample)
    f_reads = "{0}/{1}.denovo_duplicates_marked.trimmed.1.fastq".format(scratch_dir, args.sample)
    r_reads = "{0}/{1}.denovo_duplicates_marked.trimmed.2.fastq".format(scratch_dir, args.sample)
    s_reads = "{0}/{1}.denovo_duplicates_marked.trimmed.singleton.fastq".format(scratch_dir, args.sample)

    # touch each of these files so they aren't removed while running!!
    run_command("touch {0}".format(f_reads))
    run_command("touch {0}".format(r_reads))
    run_command("touch {0}".format(s_reads))  

    ## run my script to filter out reads with Ns
    f_reads_trimmed = "{0}/{1}.1.ntrimmed.fastq".format(scratch_dir, args.sample)
    r_reads_trimmed = "{0}/{1}.2.ntrimmed.fastq".format(scratch_dir, args.sample)
    s_reads_trimmed = "{0}/{1}.singletons.ntrimmed.fastq".format(scratch_dir, args.sample)
    cmd = "/home/jorvis/git/biocode/fastq/filter_fastq_by_N_content.py -l {0} -r {1} -s {2} -lo {3} -ro {4} -so {5} -p {7} -or {8}/{6}.ntrimming.report".format(f_reads, r_reads, s_reads, f_reads_trimmed, r_reads_trimmed, s_reads_trimmed, args.sample, FASTQ_FILTERING_N_PCT_CUTOFF, scratch_dir)
    run_command(cmd)

    # run bowtie2p
    # /usr/local/packages/bowtie2-2.2.4/bowtie2 -1 <m1> -2 <m2> -U <r>} -S <sam>

    mapped_list = "{0}/bam_files_to_merge.mapped.list".format(scratch_dir)
    mapped_list_to_merge = open(mapped_list, 'wt')
    unmapped_list = "{0}/bam_files_to_merge.unmapped.list".format(scratch_dir)
    unmapped_list_to_merge = open(unmapped_list, 'wt')

    # each of these could be different threads
    for index_i in range(1, 21):
        sam_file_base = "{0}/{1}.vs_metaref_seeds.fragment{2}.bowtie2".format(scratch_dir, args.sample, index_i)
        cmd = "{0} -p {5} -a -x /usr/local/projects/dacc/jorvis/read_to_metaref_alignments/partition20/mumi.20150105.genomic.withNs.fna.part{6} -1 {1} -2 {2} -U {3} -S {4}.sam".format(BOWTIE_PATH, f_reads_trimmed, r_reads_trimmed, s_reads_trimmed, sam_file_base, args.cpu_cores, index_i)
        run_command(cmd)

        # Convert SAM to BAM
        cmd = "samtools view -bS {0}.sam > {0}.bam".format(sam_file_base)
        run_command(cmd)

        # Delete SAM
        cmd = "rm {0}.sam".format(sam_file_base)
        run_command(cmd)

        # Sort BAM
        # May have to set the -m option here to limit the memory used here.
        cmd = "samtools sort -@ 4 -m 3G {0}.bam {0}.sorted".format(sam_file_base)
        run_command(cmd)

        # Delete unsorted BAM
        cmd = "rm {0}.bam".format(sam_file_base)
        run_command(cmd)

        # Write a file of just the mapped reads
        cmd = "samtools view -h -F 4 -b {0}.sorted.bam > {0}.sorted.mapped.bam".format(sam_file_base)
        run_command(cmd)
        mapped_list_to_merge.write("{0}.sorted.mapped.bam\n".format(sam_file_base))

        # Write a file of just the unmapped reads
        cmd = "samtools view -h -f 4 -b {0}.sorted.bam > {0}.sorted.unmapped.bam".format(sam_file_base)
        run_command(cmd)
        unmapped_list_to_merge.write("{0}.sorted.unmapped.bam\n".format(sam_file_base))

        # Delete the full file
        cmd = "rm {0}.sorted.bam".format(sam_file_base)
        run_command(cmd)

    # merge the mapped BAM files (attempt using samtools):
    #"{0}/{1}.vs_metaref_seeds.fragment{2}.bowtie2".format(scratch_dir, args.sample, index_i)
    #cmd = "samtools merge {0}/{1}.vs_metaref_seeds.mapped.bowtie2.sorted.bam {0}/{1}.*.fragment.*.mapped.bam".format(scratch_dir, args.sample)
    #run_command(cmd)

    # merge the mapped BAM files (attempt using Picard):
    #cmd = "java -jar -Xmx12g /usr/local/packages/picard-tools-1.115/MergeSamFiles.jar TMP_DIR={0}/{1} MERGE_SEQUENCE_DICTIONARIES=true USE_THREADING=true OUTPUT={0}/{1}.picard.mapped.merged.bam ".format(scratch_dir, args.sample)
    #for idx in range(1, 21):
    #    partial_bamfile = "{0}/{1}.vs_metaref_seeds.fragment{2}.bowtie2.sorted.mapped.bam".format(scratch_dir, args.sample, idx)
    #    cmd += "INPUT={0} ".format(partial_bamfile)

    mapped_list_to_merge.close()
    unmapped_list_to_merge.close()

    # merge the mapped BAM files (my own script), then delete them:
    cmd = "/home/jorvis/git/biocode/general/merge_bam_files.py -i {0} -o {1}/{2}.mapped.merged".format(mapped_list, scratch_dir, args.sample)
    run_command(cmd)
    run_command("rm {0}/*.fragment*.mapped.bam".format(scratch_dir))

    # merge the unmapped BAM files (my own script), then delete them:
    cmd = "/home/jorvis/git/biocode/general/merge_bam_files.py -i {0} -o {1}/{2}.unmapped.merged".format(unmapped_list, scratch_dir, args.sample)
    run_command(cmd)
    run_command("rm {0}/*.fragment*.unmapped.bam".format(scratch_dir))

    # Sort the merged BAM file, then delete the unsorted one
    cmd = "samtools sort {0}/{1}.mapped.merged.bam {0}/{1}.mapped.merged.sorted".format(scratch_dir, args.sample)
    run_command(cmd)
    run_command("rm {0}/{1}.mapped.merged.bam".format(scratch_dir, args.sample))

    cmd = "samtools sort -@ 4 -m 3G {0}/{1}.unmapped.merged.bam {0}/{1}.unmapped.merged.sorted".format(scratch_dir, args.sample)
    run_command(cmd)
    run_command("rm {0}/{1}.unmapped.merged.bam".format(scratch_dir, args.sample))

    # Done - migrate results
    result_dir = "{0}/{1}".format(COMPLETION_BASE, args.sample)
    run_command("mkdir {0}".format(result_dir))
    
    run_command("mv {2}/{0}.mapped.merged.sorted.bam {2}/{0}.unmapped.merged.sorted.bam {1}/".format(args.sample, result_dir, scratch_dir))
    run_command("rm -rf {0}".format(scratch_dir))
    

    
def run_command(cmd):
    print("INFO: [{1}] Running command: {0}\n".format(cmd, datetime.datetime.now()), flush=True)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
       raise Exception("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))
    


if __name__ == '__main__':
    main()
