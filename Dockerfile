# SGG modified from https://github.com/Bioconductor/bioconductor_docker/blob/master/Dockerfile

FROM rocker/rstudio:3.6.2

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
	gdb \
	libxml2-dev \
	python-pip \
	libz-dev \
	liblzma-dev \
	libbz2-dev \
	libpng-dev \
	libmariadb-dev \
	libjpeg62-turbo-dev \
	libudunits2-dev \
	libmagick++-dev \
	imagemagick \
	&& rm -rf /var/lib/apt/lists/*

RUN echo "R_LIBS=/usr/local/lib/R/host-site-library:\${R_LIBS}" > /usr/local/lib/R/etc/Renviron.site \
	&& echo "options(defaultPackages=c(getOption('defaultPackages'),'BiocManager'))" >> /usr/local/lib/R/etc/Rprofile.site

ADD install.R /tmp/

RUN R -f /tmp/install.R

# DEVEL: Add sys env variables to DEVEL image
RUN curl -O https://raw.githubusercontent.com/Bioconductor/BBS/master/3.11/Renviron.bioc \
	&& cat Renviron.bioc | grep -o '^[^#]*' | sed 's/export //g' >>/etc/environment \
	&& cat Renviron.bioc >> /root/.bashrc \
	&& rm -rf Renviron.bioc

# Init command for s6-overlay
CMD ["/init"]