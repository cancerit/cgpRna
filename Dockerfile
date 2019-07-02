# Tool version args, set them here so that both building stages can use them
## VAGrENG dependcies
ARG VER_VCFTOOLS="0.1.16"
ARG VER_Set_IntervalTree="0.12"
ARG VER_BEDTOOLS="2.25.0-1"
## CancerIT dependencies
ARG VER_CGPVCF="v2.2.1"
ARG VER_GRASS="v2.1.1"
ARG VER_VAGRENT="v3.3.3"
## cgpRna dependencies
ARG VER_File_ShareDir_Install="0.13"
ARG VER_Config_IniFiles="3.000002"
ARG VER_STAR="2.5.0c"
ARG VER_STARFUSION="v0.1.1"
ARG VER_BOWTIE1="1.1.2-3"
ARG VER_BOWTIE2="2.2.6-2"
ARG VER_TOPHAT="2.1.0+dfsg-1build1"
ARG VER_BLAST="2.2.31-4"
ARG VER_DEFUSE="v0.8.2"
ARG VERSION_DEFUSE="0.7.0"
ARG SOURCE_FATOTWOBIT="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit"
ARG SOURCE_BLAT="http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip"
ARG VER_GMAP="2015-12-31.v7-1"
ARG VER_RSEQC="3.0.0"
ARG VER_HTSEQ="0.7.2"

FROM  quay.io/wtsicgp/dockstore-cgpmap:3.1.4 as builder
USER  root

# Declare versions that will be installed in this stage
ARG VER_VCFTOOLS
ARG VER_Set_IntervalTree
ARG VER_BEDTOOLS
ARG VER_CGPVCF
ARG VER_GRASS
ARG VER_VAGRENT
ARG VER_File_ShareDir_Install
ARG VER_Config_IniFiles
ARG VER_STAR
ARG VER_STARFUSION
ARG VER_DEFUSE
ARG VERSION_DEFUSE
ARG SOURCE_FATOTWOBIT
ARG SOURCE_BLAT

ENV VER_VCFTOOLS $VER_VCFTOOLS
ENV VER_Set_IntervalTree $VER_Set_IntervalTree
ENV VER_BEDTOOLS $VER_BEDTOOLS
ENV VER_CGPVCF $VER_CGPVCF
ENV VER_GRASS $VER_GRASS
ENV VER_VAGRENT $VER_VAGRENT
ENV VER_File_ShareDir_Install $VER_File_ShareDir_Install
ENV VER_Config_IniFiles $VER_Config_IniFiles
ENV VER_STAR $VER_STAR
ENV VER_STARFUSION $VER_STARFUSION
ENV VER_DEFUSE $VER_DEFUSE
ENV VERSION_DEFUSE $VERSION_DEFUSE
ENV SOURCE_FATOTWOBIT $SOURCE_FATOTWOBIT
ENV SOURCE_BLAT $SOURCE_BLAT

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
      version="1.0.0" \
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

# declare so can use the defaults set at the top of this file in the ENV instructions
ARG VER_VCFTOOLS
ARG VER_BEDTOOLS
ARG VER_BOWTIE1
ARG VER_BOWTIE2
ARG VER_TOPHAT
ARG VER_BLAST
ARG VER_GMAP
ARG VER_RSEQC
ARG VER_HTSEQ

ENV VER_VCFTOOLS $VER_VCFTOOLS
ENV VER_BEDTOOLS $VER_BEDTOOLS
ENV VER_BOWTIE1 $VER_BOWTIE1
ENV VER_BOWTIE2 $VER_BOWTIE2
ENV VER_TOPHAT $VER_TOPHAT
ENV VER_BLAST $VER_BLAST
ENV VER_GMAP $VER_GMAP
ENV VER_RSEQC $VER_RSEQC
ENV VER_HTSEQ $VER_HTSEQ

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
