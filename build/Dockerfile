FROM ubuntu:latest

MAINTAINER Joshua Orvis <jorvis@gmail.com>
LABEL Description="This image is used to run any bioinformatics utilities/scripts in the Biocode collection" Version="1.2"

RUN apt-get -y update && apt-get install -y git libmldbm-perl python3
ENV PYTHONPATH /opt/biocode/lib
ENV PATH /opt/biocode/fasta:/opt/biocode/fastq:/opt/biocode/gff:/opt/biocode/genbank:/opt/biocode/general:/opt/biocode/gtf:/opt/biocode/hmm:/opt/biocode/sam_bam:/opt/biocode/sysadmin:/opt/biocode/taxonomy:$PATH
RUN git clone https://github.com/jorvis/biocode.git /opt/biocode
