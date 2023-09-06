FROM --platform=linux/amd64 rocker/tidyverse:4.2
MAINTAINER gengz@chop.edu
WORKDIR /rocker-build/

COPY script/install_bioc.r .

COPY script/install_github.r .

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog
RUN apt-get install -y --no-install-recommends libxt6

# Install dev libraries and curl
RUN apt update && apt install -y zlib1g-dev \
	libncurses5-dev \
	libbz2-dev \
	liblzma-dev \
	libcurl4-openssl-dev \
	libssl-dev curl
	

# install R packages
RUN ./install_bioc.r \
	Biobase \
	BiocManager \
	corrplot \
	cowplot \
	cutpointr \
	ggpubr \
	ggthemes \
	ggstatsplot \
	ggfortify \
	ggrepel \
	GenomicFeatures \
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
    
RUN installGithub.r jokergoo/ComplexHeatmap 

# Required for mapping segments to genes
# Add bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -zxvf bedtools-2.28.0.tar.gz && rm -f bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make && \
    mv bin/* /usr/local/bin && \
    cd .. && rm -rf bedtools2
    
# add annoFusedata
RUN R -e "remotes::install_github('d3b-center/annoFusedata', ref = '321bc4f6db6e9a21358f0d09297142f6029ac7aa', dependencies = TRUE)"


# Specify the version of circlize, same to what we used in OpenPedCan
RUN R -e "remotes::install_github('jokergoo/circlize', ref = 'b7d86409d7f893e881980b705ba1dbc758df847d', dependencies = TRUE)"
    
ADD Dockerfile .
    