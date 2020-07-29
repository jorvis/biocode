FROM ubuntu:18.04

MAINTAINER Joshua Orvis <jorvis@gmail.com>
LABEL Description="This image is used to run any bioinformatics utilities/scripts in the Biocode collection"

# https://serverfault.com/a/949998
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && apt-get install tzdata
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
RUN dpkg-reconfigure --frontend noninteractive tzdata

# From the best practices guide:
#  https://docs.docker.com/engine/userguide/eng-image/dockerfile_best-practices/#run
# libblas-dev liblapack-dev are for R modules
RUN apt-get install -y wget libmldbm-perl python3 python3-pip r-base zlib1g-dev libblas-dev liblapack-dev libxml2-dev \
 && rm -rf /var/lib/apt/lists/*

RUN cpan JORVIS/TIGR/TIGR-0.03.tar.gz

RUN Rscript -e "install.packages('tidyr', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('ggplot2', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('ade4', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('ape', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('igraph', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('vegan', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('readr', repos='http://cran.rstudio.com/')"
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("phyloseq")'

RUN wget -O /usr/bin/report_top20_metaphlan_table.R https://raw.githubusercontent.com/jorvis/biocode/master/general/report_top20_metaphlan_table.R
RUN chmod +x /usr/bin/report_top20_metaphlan_table.R

RUN wget -O /usr/bin/visualize_go_slim_mapping.R https://raw.githubusercontent.com/jorvis/biocode/master/general/visualize_go_slim_mapping.R
RUN chmod +x /usr/bin/visualize_go_slim_mapping.R

RUN pip3 install jinja2
RUN pip3 install biocode==0.9.0
RUN pip3 install numpy
