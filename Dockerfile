# Based on https://github.com/Bioconductor/bioconductor_docker/blob/master/Dockerfile
FROM rocker/rstudio:3.6.3

# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*

# issues with '/var/lib/dpkg/available' not found
# this will recreate
RUN dpkg --clear-avail

# This is to avoid the error
# 'debconf: unable to initialize frontend: Dialog'
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && \
	apt-get -y --no-install-recommends install --fix-missing \
	libxml2-dev \
	libmagick++-dev \
	libudunits2-dev \
	libgdal-dev gdal-bin libgdal20 \
	mrbayes iqtree \
	python3-pip \
	curl \
	default-jre \
	&& rm -rf /var/lib/apt/lists/* \
	&& install2.r -s -e -r https://mran.microsoft.com/snapshot/2020-02-28 -r http://bioconductor.org/packages/3.10/bioc BH BiocManager DBI DescTools KernSmooth LearnBayes MASS Matrix R.methodsS3 R.oo R.utils R6 RColorBrewer RCurl RNeXML RSQLite Rcpp RcppArmadillo XML ade4 adegenet adephylo animation ape askpass assertthat backports bit bit64 bitops blob boot callr class classInt cli cluster clusterGeneration coda colorspace combinat crayon curl dbplyr deldir desc digest dplyr e1071 ellipse ellipsis evaluate expm fansi farver fastmap fastmatch formatR futile.logger futile.options gdata ggplot2 glue gmodels gtable gtools hms htmltools httpuv httr igraph isoband jsonlite labeling lambda.r later lattice lazyeval lifecycle magick magrittr maps matrixStats mcmcse memoise mgcv mime mnormt munsell mvtnorm nlme numDeriv openssl pegas permute phangorn phylobase phylotools phytools pillar pixmap pkgbuild pkgconfig pkgload plogr plotrix plyr praise prettyunits processx progress promises ps purrr quadprog rappdirs raster rentrez reshape2 rlang rncl rprojroot rstudioapi scales scatterplot3d segmented seqRFLP seqinr sf shiny snow sourcetools sp spData spdep stringi stringr sys testthat tibble tidyr tidyselect units utf8 uuid vctrs vegan viridisLite withr xml2 xtable \
	&& R -e 'BiocManager::install(c("AnnotationDbi", "BSgenome", "Biobase", "BiocFileCache", "BiocGenerics", "BiocParallel", "BiocVersion", "Biostrings", "DelayedArray", "GenomeInfoDb", "GenomeInfoDbData", "GenomicAlignments", "GenomicFeatures", "GenomicRanges", "IRanges", "Rhtslib", "Rsamtools", "S4Vectors", "SummarizedExperiment", "VariantAnnotation", "XVector", "biomaRt", "genbankr", "rtracklayer", "zlibbioc"))' \
	&& python3 -m pip install setuptools wheel pandas numpy beautifulsoup4 requests lxml \
	&& mkdir /usr/local/seqgen && cd /usr/local/seqgen \
	&& curl -LO https://github.com/rambaut/Seq-Gen/archive/1.3.4.tar.gz \
	&& tar -xvzf 1.3.4.tar.gz && rm /usr/local/seqgen/1.3.4.tar.gz \
	&& cd Seq-Gen-1.3.4/source \
	&& make \
	&& export PATH=${PATH}:/usr/local/seqgen/Seq-Gen-1.3.4/source \
	&& mkdir -p /VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/ \
	&& curl -LO https://github.com/PangeaHIV/VirusTreeSimulator/raw/master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar


COPY . /transmissionpairs

# Install alamos-extract for load_hiv, add symbolic link for RStudio
WORKDIR /transmissionpairs/DataCollation/LANLRetrieval/alamos-extract
RUN python3 setup.py develop && ln -s /transmissionpairs /home/rstudio/transmissionpairs

WORKDIR /transmissionpairs
