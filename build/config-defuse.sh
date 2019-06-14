#! /bin/bash

set -xe

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path such as /opt/ICGC"
  exit 1
fi

INST_PATH=$1
INST_PATH=$(realpath $INST_PATH)
CONFIG_FILE="$INST_PATH/bin/defuse_install/scripts/config.txt"

if [ ! -f $CONFIG_FILE ]; then
  echo "Could not find config file: $CONFIG_FILE. Wrong installation path?"
  exit 1
fi

# how many spaces in the sed matching patter is important, do not change without close look.
sed -i \
-e "/source_directory  /c source_directory = $INST_PATH/bin/defuse_install" \
-e "/samtools_bin/c samtools_bin = $(which samtools)" \
-e "/bowtie_bin/c bowtie_bin = $(which bowtie)" \
-e "/bowtie_build_bin/c bowtie_build_bin = $(which bowtie-build)" \
-e "/blat_bin/c blat_bin = $(which blat)" \
-e "/fatotwobit_bin/c fatotwobit_bin = $(which faToTwoBit)" \
-e "/r_bin/c r_bin = $(which R)" \
-e "/rscript_bin/c rscript_bin = $(which Rscript)" \
-e "/gmap_bin/c gmap_bin = $(which gmap)" \
-e "/gmap_build_bin/c gmap_build_bin = $(which gmap_build)" \
$CONFIG_FILE
