# cgpVAF

Calculates the Variant Allele Fraction for variants sites in VCF and/or BED file
This script performs comparative analysis of variant sites in multiple tumour/normal samples in an individual.
Also facilitates the merging of varinats sites across the samples in a sample group defined by set of related samples in an individual and provides unbiased pileup[MNV] and exonerate[Indel] output for each variant site.

## Dependencies
* Please install htslib and make sure environment variable HTSLIB\_DIR is set pointing to directory containing libhts.a
## Installation
* Download current installer version from git repository
*  ./setup.sh <INSTALL_PATH> 


