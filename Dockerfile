FROM ubuntu:16.04 as builder

USER  root

ARG VER_BIODBHTS="2.10"
ARG VER_HTSLIB="1.9"
ARG VER_SAMTOOLS="1.9"
ENV VER_VCFTOOLS="0.1.16"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
curl \
ca-certificates \
libperlio-gzip-perl \
make \
bzip2 \
gcc \
psmisc \
time \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
libcurl4-gnutls-dev \
libncurses5-dev \
libgnutls-dev \
libgd-dev \
libdb-dev \
g++

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN mkdir -p $OPT/bin

WORKDIR /install_tmp

ADD build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

COPY . .
RUN bash build/opt-build-local.sh $OPT

FROM ubuntu:16.04

LABEL maintainer="cgphelp@sanger.ac.uk"\
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Sanger Institute" \
      version="5.7.0" \
      description="vafCorrect"

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL C

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
curl \
ca-certificates \
libperlio-gzip-perl \
bzip2 \
psmisc \
time \
zlib1g \
liblzma5 \
libncurses5 \
exonerate \
openssl \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

#RUN apt-get install -yq libdevel-nytprof-perl

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
