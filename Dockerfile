# The two args will be passed to the later build stage as envs, which are used by opt-build-sys-dependencies.sh
ARG VER_RSEQC="3.0.0"
ARG VER_HTSEQ="0.7.2"

FROM  quay.io/wtsicgp/dockstore-cgpmap:3.1.4 as builder
USER  root

# tool versions used by opt-build.sh
## VAGrENG dependcies
ENV VER_VCFTOOLS="0.1.16"
ENV VER_Set_IntervalTree="0.12"
ENV SOURCE_BEDTOOLS="https://github.com/arq5x/bedtools2/releases/download/v2.22.1/bedtools-2.22.1.tar.gz"
## CancerIT dependencies
ENV VER_CGPVCF="v2.2.1"
ENV VER_GRASS="v2.1.1"
ENV VER_VAGRENT="v3.3.3"
## cgpRna dependencies
ENV VER_File_ShareDir_Install="0.13"
ENV VER_Config_IniFiles="3.000002"
ENV VER_STAR="2.5.0c"
ENV VER_STARFUSION="v0.1.1"
ENV SOURCE_BOWTIE1="https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.1/bowtie-1.1.1-linux-x86_64.zip/download"
ENV VERSION_BOWTIE1="1.1.1"
ENV SOURCE_BOWTIE2="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip/download"
ENV VERSION_BOWTIE2="2.2.3"
ENV VER_TOPHAT="2.1.0"
ENV SOURCE_BLASTN="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz"
ENV VER_DEFUSE="v0.8.2"
ENV VERSION_DEFUSE="0.7.0"
ENV SOURCE_FATOTWOBIT="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit"
ENV SOURCE_BLAT="http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip"
ENV SOURCE_GMAP="http://research-pub.gene.com/gmap/src/gmap-gsnap-2015-09-10.tar.gz"

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
ARG VER_RSEQC
ARG VER_HTSEQ
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
