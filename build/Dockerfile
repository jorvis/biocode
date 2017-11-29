FROM ubuntu:yakkety

MAINTAINER Joshua Orvis <jorvis@gmail.com>
LABEL Description="This image is used to run any bioinformatics utilities/scripts in the Biocode collection"

# From the best practices guide:
#  https://docs.docker.com/engine/userguide/eng-image/dockerfile_best-practices/#run
RUN apt-get -y update && apt-get install -y \
    libmldbm-perl \
    pandoc \
    python3 \
    python3-pip \
 && rm -rf /var/lib/apt/lists/*

RUN cpan JORVIS/TIGR/TIGR-0.03.tar.gz

RUN pip3 install pypandoc jinja2
RUN pip3 install biocode==0.4.8
