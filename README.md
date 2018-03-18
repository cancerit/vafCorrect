# cgpVAF

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][travis-master]][travis-base] | [![Develop Badge][travis-develop]][travis-base] |


<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [vafCorrect](#cgpVaf)
- [INSTALL](#install)
	- [Dependencies](#package-dependencies)
- [Creating a release](#creating-a-release)
	- [Preparation](#preparation)
	- [Cutting the release](#cutting-the-release)
- [LICENSE](#License)

<!-- /TOC -->

## vafCorrect

Calculates the Variant Allele Fraction for variants sites in VCF and/or BED file 
This script performs comparative analysis of variant sites in multiple tumour/normal samples in an individual.
Also facilitates the merging of varinats sites across the samples in a sample group defined by 
set of related samples in an individual and provides unbiased pileup[MNV] and exonerate[Indel] output for each variant site.

## INSTALL 

### Dependencies
Some of the code included in this package has dependencies on several C packages:

 * [Samtools] ( max 0.1.20 until perl bindings are updated)
 * [vcftools]
 * [Bio-HTS] 
 * [Exonerate]

 `setup.sh` should install these dependencies.  

```
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

<!-- References -->
[Samtools]: https://github.com/samtools/samtools 
[vcftools]: http://vcftools.sourceforge.net/
[Bio-HTS]: https://github.com/Ensembl/Bio-DB-HTS
[Exonerate]: http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
[vafCorrect-releases]: https://github.com/cancerit/vafCorrect/releases
[travis-master-badge]: https://travis-ci.org/cancerit/vafCorrect.svg?branch=master
[travis-develop-badge]: https://travis-ci.org/cancerit/vafCorrect.svg?branch=develop
[travis-repo]: https://travis-ci.org/cancerit/vafCorrect
