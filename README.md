# cgpVAF

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][travis-master-badge]][travis-repo] | [![Develop Badge][travis-develop-badge]][travis-repo] |

* [vafCorrect](#vafcorrect)
* [Quick installation](#quick-installation)
	* [Skipping all external dependencies](#skipping-all-external-dependencies)
	* [Skipping exonerate install](#skipping-exonerate-install)
* [INSTALL](#install)
	* [Dependencies](#dependencies)
* [Creating a release](#creating-a-release)
	* [Preparation](#preparation)
	* [Cutting the release](#cutting-the-release)
* [LICENCE](#licence)

## vafCorrect

Calculates the Variant Allele Fraction for variants sites in VCF and/or BED file
This script performs comparative analysis of variant sites in multiple tumour/normal samples in an individual.
Also facilitates the merging of varinats sites across the samples in a sample group defined by
set of related samples in an individual and provides unbiased pileup[MNV] and exonerate[Indel] output for each variant site.

## Quick installation

```bash
./setup.sh path_to_install_to
```

### Skipping all external dependencies

If you want to only install vafCorrect and use existing versions of
tools from your path run as:

```bash
./setup.sh path_to_install_to 1
```

`vcftools` and `htslib` still need to install in this instances as they have perl module bindings
which are required dependencies.

### Skipping exonerate install

Central install via package manager of 2.2.0 is adequate. To skip just exonerate install run:

```bash
./setup.sh path_to_install_to 2
```

## INSTALL

### Dependencies

Some of the code included in this package has dependencies on several C packages:

* [Samtools]
* [vcftools]
* [Bio-HTS]
* [Exonerate]

`setup.sh` should install these dependencies.

```bash
./setup.sh /some/install/location
```

Please be aware that this expects basic C compilation libraries and
tools to be available, most are listed in `INSTALL`.

## Creating a release

### Preparation

* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

### Cutting the release

1. Update `perl/lib/Sanger/CGP/Vaf.pm` to the correct version.
2. Update `CHANGES.md` to show major items.
3. Run `./perl/prerelease.sh`
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to draft a release.

## LICENCE

```
Copyright (c) 2017-2019 Genome Research Ltd.

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

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
2009, 2010, 2011, 2012’."
```

<!-- References -->
[Samtools]: https://github.com/samtools/samtools
[vcftools]: http://vcftools.sourceforge.net/
[Bio-HTS]: https://github.com/Ensembl/Bio-DB-HTS
[Exonerate]: http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
[vafCorrect-releases]: https://github.com/cancerit/vafCorrect/releases
[travis-master-badge]: https://travis-ci.org/cancerit/vafCorrect.svg?branch=master
[travis-develop-badge]: https://travis-ci.org/cancerit/vafCorrect.svg?branch=dev
[travis-repo]: https://travis-ci.org/cancerit/vafCorrect
