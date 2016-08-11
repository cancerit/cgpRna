#!/bin/bash

########## LICENCE ##########
# Copyright (c) 2014-2016 Genome Research Ltd.
#
# Author: Cancer Genome Project <cgpit@sanger.ac.uk>
#
# This file is part of cgpRna.
#
# cgpRna is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’."
########## LICENCE ##########

SOURCE_STAR="https://github.com/alexdobin/STAR/archive/2.5.0c.tar.gz"
SOURCE_STARFUSION="https://github.com/STAR-Fusion/STAR-Fusion/archive/v0.1.1.tar.gz"
SOURCE_RSEQC="http://sourceforge.net/projects/rseqc/files/RSeQC-2.6.3.tar.gz/download"
SOURCE_BOWTIE1="https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.1/bowtie-1.1.1-linux-x86_64.zip/download"
VERSION_BOWTIE1="1.1.1"
SOURCE_BOWTIE2="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip/download"
VERSION_BOWTIE2="2.2.3"
SOURCE_TOPHAT="http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.0.Linux_x86_64.tar.gz"
SOURCE_BLASTN="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz"
SOURCE_DEFUSE="https://bitbucket.org/dranew/defuse/get/v0.7.0.tar.gz"
VERSION_DEFUSE="0.7.0"
SOURCE_GMAP="http://research-pub.gene.com/gmap/src/gmap-gsnap-2015-09-10.tar.gz"
SOURCE_BLAT="http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip"
SOURCE_FATOTWOBIT="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit"
SOURCE_BEDTOOLS="https://github.com/arq5x/bedtools2/releases/download/v2.21.0/bedtools-2.21.0.tar.gz"
SOURCE_HTSEQ="https://pypi.python.org/packages/3c/6e/f8dc3500933e036993645c3f854c4351c9028b180c6dcececde944022992/HTSeq-0.6.1p1.tar.gz"

done_message () {
    if [ $? -eq 0 ]; then
        echo " done."
        if [ "x$1" != "x" ]; then
            echo $1
        fi
    else
        echo " failed.  See setup.log file for error messages." $2
        echo "    Please check INSTALL file for items that should be installed by a package manager"
        exit 1
    fi
}

get_distro () {
  EXT=""
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
  elif [[ $2 == *.zip* ]] ; then
    EXT="zip"
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi
  if hash curl 2>/dev/null; then
    curl -sS -o $1.$EXT -L $2
  else
    wget -nv -O $1.$EXT $2
  fi
}

get_file () {
# output, source
  if hash curl 2>/dev/null; then
    curl --insecure -sS -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/cgpRna"
  exit 0
fi

INST_PATH=$1


CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

# get current directory
INIT_DIR=`pwd`

# re-initialise log file
echo > $INIT_DIR/setup.log

# log information about this system
(
    echo '============== System information ===='
    set -x
    lsb_release -a
    uname -a
    sw_vers
    system_profiler
    grep MemTotal /proc/meminfo
    set +x
    echo; echo
) >>$INIT_DIR/setup.log 2>&1

set -eu

# cleanup inst_path
mkdir -p $INST_PATH/bin
mkdir -p $INST_PATH/config
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
unset PERL5LIB
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"

# Set PYTHONPATH as well so that RSeQC can be installed
unset PYTHONPATH
PYTHONROOT=$INST_PATH/lib/python2.7/site-packages
mkdir -p $PYTHONROOT
export PYTHONPATH="$PYTHONROOT"

CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' PCAP`
if [[ "x$CHK" == "x" ]] ; then
  echo "PREREQUISITE: Please install PCAP-core before proceeding: https://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases"
  exit 1;
fi

CHK2=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Sanger::CGP::Vagrent`
if [[ "x$CHK2" == "x" ]] ; then
	echo "PREREQUISITE: Please install VAGrENT before proceeding: https://github.com/cancerit/VAGrENT/releases"
	exit 1;
fi

CHK3=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Sanger::CGP::Grass`
if [[ "x$CHK3" == "x" ]] ; then
	echo "PREREQUISITE: Please install Grass before proceeding: https://github.com/cancerit/grass/releases"
	exit 1;
fi

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

## grab cpanm:
rm -f $SETUP_DIR/cpanm
get_file $SETUP_DIR/cpanm http://xrl.us/cpanm
chmod +x $SETUP_DIR/cpanm

perlmods=( "Set::IntervalTree" )

set -e
for i in "${perlmods[@]}" ; do
  echo -n "Installing build prerequisite $i..."
  (
    set -x
    perl $SETUP_DIR/cpanm -v --mirror http://cpan.metacpan.org -l $INST_PATH $i
    set +x
    echo; echo
  ) >>$INIT_DIR/setup.log 2>&1
  done_message "" "Failed during installation of $i."
done

# Install STAR
echo -n "Installing STAR ..."
if [ -e $SETUP_DIR/star.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "star" $SOURCE_STAR
  mkdir -p star
  tar --strip-components 1 -C star -zxf star.tar.gz
  cp star/bin/Linux_x86_64_static/STAR $INST_PATH/bin/.
  touch $SETUP_DIR/star.success
  rm -rf $SETUP_DIR/star
  rm -f $SETUP_DIR/star.tar.gz
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build STAR."

# Install STAR-Fusion
echo -n "Installing STAR-Fusion ..."
if [ -e $SETUP_DIR/starfusion.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "starfusion" $SOURCE_STARFUSION
  mkdir -p starfusion
  tar --strip-components 1 -C starfusion -zxf starfusion.tar.gz
  cp starfusion/STAR-Fusion $INST_PATH/bin/.
  cp starfusion/lib/* $PERLROOT/.
  touch $SETUP_DIR/starfusion.success
  rm -rf $SETUP_DIR/starfusion
  rm -f $SETUP_DIR/starfusion.tar.gz
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build STAR-Fusion."

# Install bowtie1
echo -n "Installing bowtie1 ..."
if [ -e $SETUP_DIR/bowtie1.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "bowtie1" $SOURCE_BOWTIE1
  unzip -qu bowtie1.zip
  cd $SETUP_DIR/bowtie-$VERSION_BOWTIE1
  cp bowtie* $INST_PATH/bin/.
  touch $SETUP_DIR/bowtie1.success
  rm -f $SETUP_DIR/bowtie1.zip
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build bowtie1."

# Install bowtie2
echo -n "Installing bowtie2 ..."
if [ -e $SETUP_DIR/bowtie2.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "bowtie2" $SOURCE_BOWTIE2
  unzip -qu bowtie2.zip
  cd $SETUP_DIR/bowtie2-$VERSION_BOWTIE2
  cp bowtie2* $INST_PATH/bin/.
  touch $SETUP_DIR/bowtie2.success
  rm -f $SETUP_DIR/bowtie2.zip
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build bowtie2."

# Install tophat
echo -n "Installing tophat ..."
if [ -e $SETUP_DIR/tophat.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "tophat" $SOURCE_TOPHAT
  mkdir -p tophat
  tar --strip-components 1 -C tophat -zxf tophat.tar.gz
  cd tophat
  rm ./AUTHORS ./README 
  cp -r ./* $INST_PATH/bin/.
  touch $SETUP_DIR/tophat.success
  rm -rf $SETUP_DIR/tophat
  rm -f $SETUP_DIR/tophat.tar.gz
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build tophat."

# Install blastn
echo -n "Installing blastn ..."
if [ -e $SETUP_DIR/blastn.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "blastn" $SOURCE_BLASTN
  mkdir -p blastn
  tar --strip-components 1 -C blastn -zxf blastn.tar.gz
  cp blastn/bin/blastn $INST_PATH/bin/.
  touch $SETUP_DIR/blastn.success
  rm -rf $SETUP_DIR/blastn
  rm -f $SETUP_DIR/blastn.tar.gz
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build blastn."

# Install defuse
echo -n "Installing defuse ..."
if [ -e $SETUP_DIR/defuse.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "defuse" $SOURCE_DEFUSE
  mkdir -p defuse
  tar --strip-components 1 -C defuse -zxf defuse.tar.gz
  cd ./defuse/scripts
  perl_path=`which perl`
  sed -i "s|^#!/usr/bin/perl|#!${perl_path}|" *.pl
  cd ../tools
  include_search=`grep "#include <map>" ./Common.h | wc -l`
  if [ $include_search -eq 0 ]; then
    sed -i 's/#include <vector>/#include <map>\n#include <vector>/' ./Common.h
  fi
  make -j$CPU
  cd ../../
  mkdir -p $INST_PATH/bin/defuse_install
  cp -r ./defuse/* $INST_PATH/bin/defuse_install
  cp ./defuse/scripts/*pm $INST_PATH/lib/perl5
  if [ -e $INST_PATH/bin/defuse.pl ]; then
    rm $INST_PATH/bin/defuse.pl
  fi
  ln -s $INST_PATH/bin/defuse_install/scripts/defuse.pl $INST_PATH/bin/defuse.pl
  touch $SETUP_DIR/defuse.success
  rm -rf $SETUP_DIR/defuse
  rm -f $SETUP_DIR/defuse.tar.gz
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build defuse."

# Install faToTwoBit
echo -n "Installing faToTwoBit ..."
if [ -e $SETUP_DIR/faToTwoBit.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  wget $SOURCE_FATOTWOBIT
  chmod +x faToTwoBit
  cp faToTwoBit $INST_PATH/bin/.
  touch $SETUP_DIR/faToTwoBit.success
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to install faToTwoBit."

# Install blat
echo -n "Installing blat ..."
if [ -e $SETUP_DIR/blat.success ] || [ -e $INST_PATH/bin/blat ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "blat" $SOURCE_BLAT
  unzip -qu blat.zip
  cd $SETUP_DIR/blatSrc
	BINDIR=$SETUP_DIR/blat/bin
  export BINDIR
  export MACHTYPE
  mkdir -p $BINDIR
  make -j$CPU
  cp $BINDIR/blat $INST_PATH/bin/.
  touch $SETUP_DIR/blat.success
  rm -f $SETUP_DIR/blat.zip
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build blat."

# Install bedtools
echo -n "Installing bedtools ..."
if [ -e $SETUP_DIR/bedtools.success ] || [ -e $INST_PATH/bin/bedtools ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "bedtools2" $SOURCE_BEDTOOLS
  mkdir -p bedtools2
  tar --strip-components 1 -C bedtools2 -zxf bedtools2.tar.gz
  make -C bedtools2 -j$CPU
  cp bedtools2/bin/* $INST_PATH/bin/.
  touch $SETUP_DIR/bedtools.success
  rm -rf $SETUP_DIR/bedtools2
  rm -f $SETUP_DIR/bedtools2.tar.gz
)>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build bedtools."

# Install gmap
echo -n "Installing gmap ..."
if [ -e $SETUP_DIR/gmap.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "gmap" $SOURCE_GMAP
  mkdir -p gmap
  tar --strip-components 1 -C gmap -zxf gmap.tar.gz
  cd gmap
  ./configure --prefix=$INST_PATH --with-gmapdb=$INST_PATH
  make -j$CPU
  make install
  touch $SETUP_DIR/gmap.success
  rm -rf $SETUP_DIR/gmap
  rm -f $SETUP_DIR/gmap.tar.gz
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build gmap."

# Install RSeQC using PYTHONPATH location set above
echo -n "Installing RSeQC ..."
if [ -e $SETUP_DIR/rseqc.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "rseqc" $SOURCE_RSEQC
  mkdir -p rseqc
  tar --strip-components 1 -C rseqc -zxf rseqc.tar.gz
  cd $SETUP_DIR/rseqc &&
  python ./setup.py install --prefix=$INST_PATH
  touch $SETUP_DIR/rseqc.success
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build RSeQC."

# install HTSeq
echo -n "Installing HTSeq ..."
if [ -e $SETUP_DIR/htseq.success ]; then
  echo -n " previously installed ...";
else
(
  cd $SETUP_DIR
  get_distro "htseq" $SOURCE_HTSEQ
  mkdir -p htseq
  tar --strip-components 1 -C htseq -zxf htseq.tar.gz 
  cd $SETUP_DIR/htseq &&
  python ./setup.py install --prefix=$INST_PATH
  touch $SETUP_DIR/htseq.success
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build HTSeq."


#add bin path for install tests
export PATH="$INST_PATH/bin:$PATH"

cd $INIT_DIR/perl

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
(
  set -x
  perl $SETUP_DIR/cpanm -v --mirror http://cpan.metacpan.org --notest -l $INST_PATH/ --installdeps . < /dev/null
  set +x
) >>$INIT_DIR/setup.log 2>&1
done_message "" "Failed during installation of core dependencies."

# Install cgpRna code
echo -n "Installing cgpRna..."
(
  cd $INIT_DIR/perl
  perl Makefile.PL INSTALL_BASE=$INST_PATH
  make
  make test
  make install
  cp $INIT_DIR/perl/config/*.ini $INST_PATH/config/
) >>$INIT_DIR/setup.log 2>&1
done_message "" "cgpRna install failed."

# clean-up
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo "Please add the following to beginning of PYTHONPATH:"
echo "  $PYTHONROOT"
echo
echo "If you intend to use the fusion pipeline, open the defuse config file: $INST_PATH/bin/defuse_install/scripts/config.txt and update the following values..."
echo "source_directory = $INST_PATH/bin/defuse_install"
echo "Then further down in the section titled # Paths to external tools..."
echo "samtools_bin = $INST_PATH/bin/samtools"
echo "bowtie_bin = $INST_PATH/bin/bowtie"
echo "bowtie_build_bin = $INST_PATH/bin/bowtie-build"
echo "blat_bin = $INST_PATH/bin/blat"
echo "fatotwobit_bin = $INST_PATH/bin/faToTwoBit"
echo "r_bin = <full path to where the R executable is installed in your environment>"
echo "rscript_bin = <full path to where the Rscript executable is installed in your environment>"
echo "gmap_bin = $INST_PATH/bin/gmap"
echo "gmap_build_bin = $INST_PATH/bin/gmap_build"
echo
echo "Finally, open the file: $INST_PATH/config/defuse.ini and update the following parameter to..."
echo "defuseversion=$VERSION_DEFUSE"

exit 0
