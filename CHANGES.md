# CHANGES

# next

* Added run-cgprna subcommands: `tophat-fusion`, `start-fusion` and `defuse` to run infuse pipeline.
* Re-organized expected file structure in Star reference bundle. Previous Star reference bundle on the FTP server will not work for this version.
* Added CWL files to run Infuse pipeline in cgpRna
* Uploaded new set of reference bundle files for GRCh38 to the FTP server.

## 2.4.1

* added a patched topaht-fusion-post script, in which the contig name bug #37 is fixed

## 2.4.0

* Revised Dockerfile so that the builder stage is properly used and the image size is reduced. Resolved #28.
* Added a CLI, so that user can complete a step in RNA-seq data workflow with just one command, which also eases the development of CWL files.
   * CLI is written in Python
   * currently 4 subcommands are implemented:
      * `map`: uses star_mapping.pl to map and marks duplicates after mapping.
      * `stats`: generates mapping stats using bam_stats and RSeQC.
      * `bigwig`: generates bigwig file using bamToBw.pl
      * `counts`: counts reads using htseq-count.
* Built a new set of reference files for CLI to use. They're available on ftp://ftp.sanger.ac.uk/pub/cancer/support-files/cgpRna_container/.
* Added CWL files:
   * added a workflow to map sample by lanes, generate mapping stats for lanes, merge lane bams, generate bigwig file and count reads.
   * added CWL tools/workflows for the workflow above to use.
   * added example JSON for most of the CWL files.

## 2.3.4

* RG tags are converted to shell safe strings before passing to Star. partially resolve #30.

## 2.3.3

* fixed verison numbers

## 2.3.2

* fixed blatSrc url in setup.sh

## 2.3.1

* fixed blatSrc url in Docker container build script

## 2.3.0

* dockerised cgpRna. Within the docker container, version of some dependent tools have been changed:

  1. Python3(**3.7.X**) is used in the container, thus:
     - RSeQC is updated from version ***2.6.4*** to ***3.0.0***.
     - However, version of HTSeq is not changed.

  1. Defuse is updated from ***v0.7.0*** to ***v0.8.2*** due to "Possible precedence issue with control flow operator" warning with the version of Perl installed in the container. The fix is [here](https://bitbucket.org/dranew/defuse/commits/b979855999b8106f5dc9f9e54f86935c7bf4f62f). However it is not merged in any ***v0.7.x*** versions (also because there's only ONE ***v0.7.x*** version --***v0.7.0***), hence ***v0.8.2***, which is the latest version at the time, was chosen. 

  1. To utilise apt packages for the ease of their installations, the following tool version changes were made:
     - BedTools: ***2.21.0*** to ***2.25.0-1*** in apt,
     - bowtie1: ***1.1.1*** to ***1.1.2-3*** in apt,
     - bowtie2: ***2.2.3*** to ***2.2.6-2*** in apt,
     - blast: ***2.2.30*** to ***2.2.31-4*** in apt,
     - gmap: ***2015-09-10*** to ***2015-12-31.v7-1*** in apt. 

* added an extra option "-updateconfig" to *defuse_fusion.pl*. It takes a file as input, which content is used to update *defuse-config.txt*. It'll search *defuse.ini* to find its default value. This addition shouldn't break any existing usage.

* *defuse_fusion.pl* will always create a *defuse-config.txt* file in its temp folder and use it to run *defuse.pl*. In the temp file, **dataset_directory** is corrected using related command line input values. If "-updateconfig" presents, values in the file will be used to overwrite corresponding values in the temp file. **Note:** if ***dataset_directory*** key presents in the "-updateconfig" file, this will be the final value in the temp config file.

* tophat-fusion-post now skips read dist step, which our infuse pipeline does not care and the `tophat-fusion-post.py` does not insert generated read distributions in final html report anyway neither.

## 2.2.2

* Change tabix query call to query_full

## 2.2.1

* Handle files not being written by tophat-fusion-post when no events generated.

## 2.2.0

* RSeqQC  - Updated to 2.6.4, fixes known issues with low levels of mapped data.
  * Also resolves 2.6.3 that is no longer available from central repos.
* HTSeq - updated to 0.7.2 to ensure compatilbility with current pysam.

## 2.1.9

* Fix missed compilation failure

## 2.1.8

* Fix issues with empty output from tophatpost when running tophat_filter

## 2.1.6-2.1.7

* Fixing issues with data exhibiting no fusions.

## 2.1.5

* updating cpanm source url

## 2.1.4

* Fixing installation of tophat. Since we no longer use a compiled repository. Bug introduced in 2.1.2.

## 2.1.3

* Correct name from sample.bas to sample.bam.bas

## 2.1.2

* Fixed tophat-fusion-post issue in local tophat fork added tophat v2.1.0a to setup.sh

## 2.1.1

* Added HTseq installer for read counting
