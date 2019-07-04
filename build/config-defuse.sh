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

echo -e \
"source_directory = $INST_PATH/bin/defuse_install\n"\
"samtools_bin = samtools\n"\
"bowtie_bin = bowtie\n"\
"bowtie_build_bin = bowtie-build\n"\
"blat_bin = blat\n"\
"fatotwobit_bin = faToTwoBit\n"\
"r_bin = R\n"\
"rscript_bin = Rscript\n"\
"gmap_bin = gmap\n"\
"gmap_build_bin = gmap_build"\
> $INST_PATH/config/defuse_running_env_config.ini
echo "updateconfig=$(realpath $INST_PATH)/config/defuse_running_env_config.ini" >> $INST_PATH/config/defuse.ini
