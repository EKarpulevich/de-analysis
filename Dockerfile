FROM r-base:4.4.1

ARG HOBOTNICA_PUBLIC_VERSION=7.0.1

ENV TERM=xterm-256color

WORKDIR /usr/gitc

RUN set -eux; \
        apt-get update; \
        apt-get install -y \
            libfreetype6-dev \
            build-essential\
            curl \
            wget \
            gfortran \
            autoconf \
            automake \
            bash \
            libbz2-dev \
            xz-utils \
            zlib1g-dev \
            gcc \
            make \
            liblzma-dev \
            musl-dev \
            ncurses-dev \
            libcurl4-openssl-dev \
            libxml2-dev \
            libmagick++-dev \
            libssl-dev \
            libudunits2-dev \
            cargo \
            libharfbuzz-dev \
            libfribidi-dev \
            libproj-dev \
            libgeos-dev \
            libgdal-dev \
            gdal-bin \
            proj-data \
            proj-bin;



RUN R --help
COPY install.R /root/
RUN cd /root; Rscript install.R
COPY run.R /root/
RUN mkdir /root/source
COPY source/ /root/source
RUN mkdir /root/hobotnica
COPY hobotnica/ /root/hobotnica
RUN mkdir /root/separate_tools
COPY separate_tools/ /root/separate_tools

