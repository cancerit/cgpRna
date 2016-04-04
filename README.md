LICENCE
=======

Copyright (c) 2014-2015 Genome Research Ltd.

Author: Cancer Genome Project <cgpit@sanger.ac.uk>

This file is part of cgpRna.

cgpRna is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."

cgpRna
======

cgpRna provides pipelines, for RNA-Seq data, that implement commonly used mapping 
and analysis programs, such as TopHat and rna-star.
At the present time (May 2016), only pipelines for mapping (with STAR), lane QC 
and fusion gene detection are included in this codebase but this will be added 
to over time with; differential expression, gene/transcript quantification, splice 
variant analysis and allele specific expression.

### Dependencies/Installation

Please install Perl packages [PCAP-core](https://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases) and [VAGrENT](https://github.com/cancerit/VAGrENT/releases) first.

Prerequisites for the [RSeQC](http://rseqc.sourceforge.net/#installation) software are:
* gcc
* [python2.7](https://www.python.org/downloads/) (The minimum version the pipeline has been tested with is python-2.7.6
* [R](https://www.r-project.org/)
* [numpy](http://www.numpy.org/)
 
Once that is done and your $PATH environment variable has been updated so that newly installed software can be found, run the following to install cgpRna:

./setup.sh path_to_install_to

N.B. the path_to_install_to should be the same as the install location used for PCAP-core and VAGrENT above.

### Tools installed by setup.sh

* Some CPAN hosted libraries, see perl/Makefile.PL
* [STAR](https://github.com/alexdobin/STAR/releases)
* [Tophat](https://ccb.jhu.edu/software/tophat/index.shtml)
* [deFuse](https://bitbucket.org/dranew/defuse)
* [RSeQC](http://rseqc.sourceforge.net)
* [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) N.B. both bowtie and bowtie2 are installed and can be used with Tophat
* [blat](http://hgwdev.cse.ucsc.edu/~kent/src/) Unless already in the install location bin directory
* [gmap](http://research-pub.gene.com/gmap/) The aligner used by deFuse
* [faToTwoBit](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/) deFuse dependency
* [bedtools](https://github.com/arq5x/bedtools2/) Unless already in the install location bin directory
* [blastn](ftp://ftp.ncbi.nlm.nih.gov/blast/) Used by tophat-fusion post

N.B. samtools is also a dependency but this is installed by PCAP-Core which should have already been installed (see above).

If you are planning to use the fusion pipeline, specifically defuse_fusion.pl, the deFuse config.txt file will need to be updated with the installed locations of a number of tools.
These paths are printed to screen if the setup.sh script completes successfully so make a note of the locations and update the file as instructed.
