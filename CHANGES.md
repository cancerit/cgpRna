# CHANGES

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
