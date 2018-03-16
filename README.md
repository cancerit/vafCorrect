# cgpVAF

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][travis-master]][travis-base] | [![Develop Badge][travis-develop]][travis-base] |


Calculates the Variant Allele Fraction for variants sites in VCF and/or BED file 
This script performs comparative analysis of variant sites in multiple tumour/normal samples in an individual.
Also facilitates the merging of varinats sites across the samples in a sample group defined by 
set of related samples in an individual and provides unbiased pileup[MNV] and exonerate[Indel] output for each variant site.


# LICENCE

Copyright (c) 2017 Genome Research Ltd.

Author: Cancer Genome Project <cgpit@sanger.ac.uk>

This file is part of vafCorrect.

vafCorrect is free software: you can redistribute it and/or modify it under
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
2009, 2010, 2011, 2012’.

### Dependencies/Install
Some of the code included in this package has dependencies on several C packages:

 * [Samtools](https://github.com/samtools/samtools) - max 0.1.20 until perl bindings are updated
 * [vcftools](http://vcftools.sourceforge.net/)
 * [Bio-HTS] (https://github.com/Ensembl/Bio-DB-HTS)
 * [Exonerate](http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)


Please use `setup.sh` to install the dependencies.  Please be aware that this expects basic C
compilation libraries and tools to be available, most are listed in `INSTALL`.


## Creating a release

#### Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

#### Cutting the release
1. Update `perl/lib/Sanger/CGP/Vaf.pm` to the correct version.
2. Update `CHANGES.md` to show major items.
3. Run `./perl/prerelease.sh`
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to draft a release.
