#! /bin/bash

set -xe

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path such as /opt/ICGC"
  exit 1
fi

INST_PATH=$1
INST_PATH=$(realpath $INST_PATH)

if [ ! -d $INST_PATH/bin/defuse_install ]; then
  echo "Could not find bin/defuse_install in $INST_PATH: Wrong installation path?"
  exit 1
fi

if [ ! -d $INST_PATH/config ]; then
  echo "Could not find folder 'config' in $INST_PATH: Wrong installation path?"
  exit 1
fi

# how many spaces in the sed matching patter is important, do not change without close look.
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
