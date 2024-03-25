FROM --platform=linux/amd64 rocker/tidyverse:4.2
MAINTAINER rokita@chop.edu
WORKDIR /rocker-build/

#RUN RSPM="https://packagemanager.rstudio.com/cran/2023-09-20" \
#  && echo "options(repos = c(CRAN='$RSPM'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

COPY script/install_bioc.r .
COPY script/install_github.r .

### Install apt-getable packages to start
#########################################
RUN apt-get update  
RUN apt-get install -y --no-install-recommends apt-utils dialog
RUN apt-get install -y --no-install-recommends \
  libxt6 \
  bzip2 

# Install dev libraries and curl
RUN apt update && apt install -y zlib1g-dev \
	libncurses5-dev \
	libbz2-dev \
	liblzma-dev \
	libcurl4-openssl-dev \
	libssl-dev \
	curl \
	cmake
	

# install R packages
RUN ./install_bioc.r \
	Biobase \
	BiocManager \
	circlize \
	corrplot \
	cowplot \
	cutpointr \
	ggalluvial \
	ggpubr \
	ggthemes \
	ggstatsplot \
	ggfortify \
	ggrepel \
	GenomicFeatures \
	GSVA \
	msigdbr \
	reshape2 \
	R.utils \
	survival \
	survminer
	
## R packages for tp53_nf1_score
RUN ./install_bioc.r \
    GenomicRanges \
    optparse \
    broom \
    data.table 
    
# Install pip3 and low-level python installation reqs
RUN apt-get update
RUN apt-get -y --no-install-recommends install \
    python3-pip  python3-dev
RUN ln -s /usr/bin/python3 /usr/bin/python    
RUN pip3 install \
    "Cython==0.29.15" \
    "setuptools==46.3.0" \
    "six==1.14.0" \
    "wheel==0.34.2" 

# Install python3 tools and dependencies
RUN pip3 install \
    "numpy==1.24.3" \
    "pandas==2.0.1" \
    "matplotlib==3.7.1" \
    "scikit-learn==1.2.2" \
    "seaborn==0.12.2" \
    "rpy2==3.5.0" \
    "utils==1.0.1" 
    
RUN installGithub.r \
  jokergoo/ComplexHeatmap \
	clauswilke/colorblindr 

# Required for mapping segments to genes
# Add bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -zxvf bedtools-2.28.0.tar.gz && rm -f bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make && \
    mv bin/* /usr/local/bin && \
    cd .. && rm -rf bedtools2

# Add bedops per the BEDOPS documentation
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.37/bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    tar -jxvf bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    rm -f bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    mv bin/* /usr/local/bin
    
# add annoFusedata
RUN R -e "remotes::install_github('d3b-center/annoFusedata', ref = '321bc4f6db6e9a21358f0d09297142f6029ac7aa', dependencies = TRUE)"

ADD Dockerfile .
    
