#! /bin/bash

set -xe

# Add python ppa
add-apt-repository -y ppa:deadsnakes/ppa
# Add R key
echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | apt-key add -

# install python3, R, bedtools, bowtie1, bowtie2, blast, gmap and other packages for RSeQC to install
apt-get update
apt-get install -yq --no-install-recommends \
python3.7 python3.7-dev \
r-base r-base-dev \
bedtools=${VER_BEDTOOLS} \
bowtie=${VER_BOWTIE1} \
bowtie2=${VER_BOWTIE2} \
ncbi-blast+=${VER_BLAST} \
gmap=${VER_GMAP} \
libcurl4-gnutls-dev zlib1g-dev
apt-get upgrade -yq gcc

# install R packages:
Rscript -e 'install.packages("ada")'  # required by Defuse

# Replace python3
rm /usr/bin/python3 && ln -s /usr/bin/python3.7 /usr/bin/python3

# install pip3
curl -s https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3 get-pip.py
rm -f get-pip.py

# for bx-python (one of RSeQC package dependencies)
apt-get install -yq --no-install-recommends liblzo2-dev

# install RSeQC and HTSeq
pip3 install RSeQC=="$VER_RSEQC"
pip3 install HTSeq=="$VER_HTSEQ"

# if use HTSeq to plot
pip3 install matplotlib

# cleanning
apt-get autoremove -yq
