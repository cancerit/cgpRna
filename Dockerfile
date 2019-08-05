FROM  quay.io/wtsicgp/dockstore-cgpmap:3.1.4 as builder
USER  root

# Tool version ENVs, some of them are also used in the second build stage, make sure version are consistent between the two stages.
## VAGrENG dependcies
ENV VER_VCFTOOLS "0.1.16"
ENV VER_Set_IntervalTree "0.12"
ENV VER_BEDTOOLS "2.25.0-1"
## CancerIT dependencies
ENV VER_CGPVCF "v2.2.1"
ENV VER_GRASS "v2.1.1"
ENV VER_VAGRENT "v3.3.3"
## cgpRna dependencies
ENV VER_File_ShareDir_Install "0.13"
ENV VER_Config_IniFiles "3.000002"
ENV VER_STAR "2.5.0c"
ENV VER_STARFUSION "v0.1.1"
ENV VER_TOPHAT "2.1.0"
ENV VER_DEFUSE "v0.8.2"
ENV SOURCE_FATOTWOBIT "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit"
ENV SOURCE_BLAT "https://hgwdev.gi.ucsc.edu/~kent/src/blatSrc35.zip"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
g++ \
make \
gcc \
pkg-config \
zlib1g-dev \
software-properties-common \
zip \
unzip \
libpng-dev \
libboost-all-dev
# libboost-all-dev is required to compile defuse. Its installation installs python2 as well, which is required for building Bedtools

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$OPT/biobambam2/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# build tools from other repos
ADD build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

# build the tools in this repo, separate to reduce build time on errors
COPY . .
RUN bash build/opt-build-local.sh $OPT

FROM ubuntu:16.04

LABEL maintainer="cgphelp@sanger.ac.uk" \
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Trust Sanger Institute" \
      version="2.3.3" \
      description="cgpRna docker"

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
p11-kit \
software-properties-common \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

# dependecy tool versions, some of them are also used in the first stage, make sure they are consistent between stages.
## VAGrENG dependcies
ENV VER_VCFTOOLS "0.1.16"
ENV VER_BEDTOOLS "2.25.0-1"
## cgpRna dependencies
ENV VER_BOWTIE1 "1.1.2-3"
ENV VER_BOWTIE2 "2.2.6-2"
ENV VER_BLAST "2.2.31-4"
ENV VER_GMAP "2015-12-31.v7-1"
ENV VER_RSEQC "3.0.0"
ENV VER_HTSEQ "0.7.2"

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$OPT/biobambam2/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

COPY build/opt-build-sys-dependencies.sh ./
COPY build/config-defuse.sh ./
RUN bash opt-build-sys-dependencies.sh && rm -f opt-build-sys-dependencies.sh && bash config-defuse.sh $OPT && rm -f config-defuse.sh

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
