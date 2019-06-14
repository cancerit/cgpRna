#! /bin/bash

set -xe

# Add python ppa
add-apt-repository -y ppa:deadsnakes/ppa
# Add R key
echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | apt-key add -

# install python3, R and other packages for RSeQC to install
apt-get update
apt-get install -yq --no-install-recommends python3.7 python3.7-dev r-base r-base-dev libcurl4-gnutls-dev zlib1g-dev
apt-get upgrade -yq gcc

# Replace python3
rm /usr/bin/python3 && ln -s /usr/bin/python3.7 /usr/bin/python3

# install pip3
curl -s https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3 get-pip.py
rm -f get-pip.py

# for bx-python (one of RSeQC package dependencies)
apt-get install -yq --no-install-recommends liblzo2-dev

# if use HTSeq to plot
apt-get install -yq --no-install-recommends python3-matplotlib

# install RSeQC and HTSeq
pip3 install RSeQC=="$VER_RSEQC"
pip3 install HTSeq=="$VER_HTSEQ"

# cleanning
apt-get autoremove -yq
