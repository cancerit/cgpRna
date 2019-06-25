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
cd $SETUP_DIR

# make sure tools installed can see the install loc of libraries
set +u
export LD_LIBRARY_PATH=`echo $INST_PATH/lib:$LD_LIBRARY_PATH | perl -pe 's/:\$//;'`
export PATH=`echo $INST_PATH/bin:$PATH | perl -pe 's/:\$//;'`
export MANPATH=`echo $INST_PATH/man:$INST_PATH/share/man:$MANPATH | perl -pe 's/:\$//;'`
export PERL5LIB=`echo $INST_PATH/lib/perl5:$PERL5LIB | perl -pe 's/:\$//;'`
set -u

# Install cgpRna
cd $INIT_DIR/perl
perl Makefile.PL INSTALL_BASE=$INST_PATH
make
make test
make install
mkdir -p $INST_PATH/config
cp $INIT_DIR/perl/config/*.ini $INST_PATH/config/

# config defuse
sed -i -e "/defuseversion/c defuseversion=$VER_DEFUSE" $INST_PATH/config/defuse.ini
echo -e \
"source_directory = $INST_PATH/bin/defuse_install\n"\
"samtools_bin = $(which samtools)\n"\
"bowtie_bin = $(which bowtie)\n"\
"bowtie_build_bin = $(which bowtie-build)\n"\
"blat_bin = $(which blat)\n"\
"fatotwobit_bin = $(which faToTwoBit)\n"\
"r_bin = $(which R)\n"\
"rscript_bin = $(which Rscript)\n"\
"gmap_bin = $(which gmap)\n"\
"gmap_build_bin = $(which gmap_build)"\
> $INST_PATH/config/defuse_running_env_config.ini
echo "updateconfig=$(realpath $INST_PATH)/config/defuse_running_env_config.ini" >> $INST_PATH/config/defuse.ini
