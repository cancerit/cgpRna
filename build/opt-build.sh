#! /bin/bash

set -xe

if [[ -z "${TMPDIR}" ]]; then
  TMPDIR=/tmp
fi

set -u

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path such as /opt/ICGC"
  exit 1
fi

# get path to this script
SCRIPT_PATH=`dirname $0`;
SCRIPT_PATH=`(cd $SCRIPT_PATH && pwd)`

# get the location to install to
INST_PATH=$1
mkdir -p $1
INST_PATH=`(cd $1 && pwd)`
echo $INST_PATH

# get current directory
INIT_DIR=`pwd`

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR/distro # don't delete the actual distro directory until the very end
mkdir -p $INST_PATH/bin
mkdir -p $INST_PATH/R-lib
mkdir -p $INST_PATH/python-lib
cd $SETUP_DIR

# make sure tools installed can see the install loc of libraries
set +u
export LD_LIBRARY_PATH=`echo $INST_PATH/lib:$LD_LIBRARY_PATH | perl -pe 's/:\$//;'`
export PATH=`echo $INST_PATH/bin:$PATH | perl -pe 's/:\$//;'`
export MANPATH=`echo $INST_PATH/man:$INST_PATH/share/man:$MANPATH | perl -pe 's/:\$//;'`
export PERL5LIB=`echo $INST_PATH/lib/perl5:$PERL5LIB | perl -pe 's/:\$//;'`
set -u


# install R packages
Rscript -e "install.packages(\"ada\", \"$INST_PATH/R-lib\")"  # required by Defuse

# install python packages
OPT_BK=$OPT # Somehow OPT affects compilation of numpy
unset OPT
pip3 install --upgrade --ignore-installed --root=$INST_PATH numpy
pip3 install --upgrade --ignore-installed --root=$INST_PATH \
  RSeQC=="$VER_RSEQC" \
  HTSeq=="$VER_HTSEQ" \
  matplotlib==3.0
OPT=$OPT_BK

## vcftools
if [ ! -e $SETUP_DIR/vcftools.success ]; then
  curl -sSL --retry 10 https://github.com/vcftools/vcftools/releases/download/v${VER_VCFTOOLS}/vcftools-${VER_VCFTOOLS}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 2 -C distro -xzf distro.tar.gz
  cd distro
  ./configure --prefix=$INST_PATH --with-pmdir=lib/perl5
  make -j$CPU
  make install
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/vcftools.success
fi

## add File::ShareDir::Install for VAGrENT
if [ ! -e $SETUP_DIR/File_ShareDir_Install.success ]; then
  cpanm -l $INST_PATH --mirror http://cpan.metacpan.org File::ShareDir::Install@$VER_File_ShareDir_Install
  touch $SETUP_DIR/File_ShareDir_Install.success.success
fi

### VAGrENT
if [ ! -e $SETUP_DIR/vagrent.success ]; then
  curl -sSL --retry 10 https://github.com/cancerit/VAGrENT/archive/${VER_VAGRENT}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -xzf distro.tar.gz
  cd distro
  cpanm --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH --installdeps .
  cpanm -v --no-interactive --mirror http://cpan.metacpan.org -l $INST_PATH .
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/vagrent.success
fi

### cgpVcf
if [ ! -e $SETUP_DIR/cgpVcf.success ]; then
  curl -sSL --retry 10 https://github.com/cancerit/cgpVcf/archive/${VER_CGPVCF}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -xzf distro.tar.gz
  cd distro
  cpanm --no-interactive --notest --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
  cpanm -v --no-interactive --mirror http://cpan.metacpan.org -l $INST_PATH .
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/cgpVcf.success
fi

### Grass
if [ ! -e $SETUP_DIR/grass.success ]; then
  curl -sSL --retry 10 https://github.com/cancerit/grass/archive/${VER_GRASS}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -xzf distro.tar.gz
  cd distro
  cpanm --no-interactive --notest --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
  cpanm -v --no-interactive --mirror http://cpan.metacpan.org -l $INST_PATH .
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/grass.success
fi

# Install STAR
if [ ! -e $SETUP_DIR/star.success ]; then
  curl -sSL --retry 10 https://github.com/alexdobin/STAR/archive/${VER_STAR}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -zxf distro.tar.gz
  cp distro/bin/Linux_x86_64_static/STAR $INST_PATH/bin/.
  touch $SETUP_DIR/star.success
fi

# Install STAR-Fusion
if [ ! -e $SETUP_DIR/starfusion.success ]; then
  curl -sSL --retry 10 https://github.com/STAR-Fusion/STAR-Fusion/archive/${VER_STARFUSION}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -zxf distro.tar.gz
  cp distro/STAR-Fusion $INST_PATH/bin/.
  cp distro/lib/* $INST_PATH/lib/perl5/.
  touch $SETUP_DIR/starfusion.success
fi

# Install tophat
if [ ! -e $SETUP_DIR/tophat.success ]; then
  curl -sSL --retry 10 https://ccb.jhu.edu/software/tophat/downloads/tophat-${VER_TOPHAT}.Linux_x86_64.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -zxf distro.tar.gz
  rm -f distro/AUTHORS distro/LICENSE distro/README
  cp -r distro/* $INST_PATH/bin/.
  # patch tophat-fusion-post
  patch $INST_PATH/bin/tophat-fusion-post $SCRIPT_PATH/patches/tophat-fusion-post.patch
  chmod +x $INST_PATH/bin/tophat-fusion-post
  touch $SETUP_DIR/tophat.success
fi

# Install defuse
if [ ! -e $SETUP_DIR/defuse.success ]; then
  curl -sSL --retry 10 "https://bitbucket.org/dranew/defuse/get/${VER_DEFUSE}.tar.gz" > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -zxf distro.tar.gz
  cd distro/tools
  include_search=`grep "#include <map>" ./Common.h | wc -l`
  if [ $include_search -eq 0 ]; then
    sed -i 's/#include <vector>/#include <map>\n#include <vector>/' ./Common.h
  fi
  make -j$CPU
  cd ../../
  cp -r distro $INST_PATH/bin/defuse_install
  cp distro/scripts/*pm $INST_PATH/lib/perl5
  if [ -e $INST_PATH/bin/defuse.pl ]; then
    rm $INST_PATH/bin/defuse.pl
  fi
  ln -s $INST_PATH/bin/defuse_install/scripts/defuse_run.pl $INST_PATH/bin/defuse.pl
  touch $SETUP_DIR/defuse.success
fi

# Install faToTwoBit
if [ ! -e $SETUP_DIR/faToTwoBit.success ]; then
  curl -sSL --retry 10 $SOURCE_FATOTWOBIT > faToTwoBit
  chmod +x faToTwoBit
  cp faToTwoBit $INST_PATH/bin/
  touch $SETUP_DIR/faToTwoBit.success
fi

# Install blat
if [ ! -e $SETUP_DIR/blat.success ]; then
  curl -sSL --retry 10 $SOURCE_BLAT > blat.zip
  rm -rf blatSrc
  unzip -qu blat.zip
  cd blatSrc
	export BINDIR=$SETUP_DIR/blat/bin
  # Blat does not recognise the startand MACHTYPE in system, so have to save it and use a shorter version of it instead
  MACHTYPE_BK=$MACHTYPE
  export MACHTYPE=x86_64
  mkdir -p $BINDIR
  make -j$CPU
  cp $BINDIR/blat $INST_PATH/bin/
  unset BINDIR
  # Restore the factory setting of MACHTYPE
  export MACHTYPE=$MACHTYPE_BK
  cd $SETUP_DIR
  touch $SETUP_DIR/blat.success
fi

## add a few perl modules for cgpRna
if [ ! -e $SETUP_DIR/Set_IntervalTree.success ]; then
  cpanm -l $INST_PATH --mirror http://cpan.metacpan.org Set::IntervalTree@$VER_Set_IntervalTree
  cpanm -l $INST_PATH --mirror http://cpan.metacpan.org Config::IniFiles@$VER_Config_IniFiles
  touch $SETUP_DIR/Set_IntervalTree.success
fi
