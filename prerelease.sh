#!/bin/bash

set -eu # exit on first error or undefined value in subtitution
set -o pipefail

# get current directory
INIT_DIR=`pwd`

# get location of this file
MY_PATH="`dirname \"$0\"`"              # relative
MY_PATH="`( cd \"$MY_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$MY_PATH" ] ; then
  # error; for some reason, the path is not accessible
  # to the script (e.g. permissions re-evaled after suid)
  echo Failed to determine location of script >2
  exit 1  # fail
fi
# change into the location of the script
cd $MY_PATH/perl
rm -rf blib
rm -rf docs

echo '### Running perl tests ###'

export HARNESS_PERL_SWITCHES=-MDevel::Cover=-db,reports,-ignore,'t/.*\.t'
mkdir -p docs/reports_text
prove -w -I lib
echo
echo '### Generating test/pod coverage reports ###'
# removed 'condition' from coverage as '||' 'or' doesn't work properly
cover -coverage branch,subroutine,pod -report_c0 50 -report_c1 85 -report_c2 100 -report html_basic reports -silent
cover -coverage branch,subroutine,pod -report text reports -silent > docs/reports_text/coverage.txt
rm -rf reports/structure reports/digests reports/cover.13 reports/runs
cp reports/coverage.html reports/index.html
mv reports docs/reports_html
unset HARNESS_PERL_SWITCHES

echo '### Generating POD ###'
mkdir -p docs/pod_html
perl -MPod::Simple::HTMLBatch -e 'Pod::Simple::HTMLBatch::go' lib:bin docs/pod_html > /dev/null

echo '### Archiving docs folder ###'
tar cz -f docs.tar.gz docs

# generate manifest, and cleanup
echo '### Generating MANIFEST ###'
# delete incase any files are moved, the make target just adds stuff
rm -f MANIFEST
# cleanup things which could break the manifest
rm -rf install_tmp
perl Makefile.PL > /dev/null
make manifest &> /dev/null
rm -f Makefile MANIFEST.bak pm_to_blib MYMETA.*

# change back to original dir
cd $INIT_DIR
