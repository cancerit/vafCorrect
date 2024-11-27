# CHANGES

## 5.7.3

* Exit with a clearer error message when file does not pass the minimum read threshold
* Update Dockerfile and setup.sh

## 5.7.2

* OS and Bio::DB::HTS updates

## 5.7.1

* Fixed bug that prevented vcf+bed locations being possible

## 5.7.0

* Switch to low level HTS API to speed up fetch calls.
* Collate mate read calls into single fetch for overlapping ranges.
* ~46% faster - see [#68](https://github.com/cancerit/vafCorrect/pull/68).

## 5.6.1

* Ensure line read rather than $_ are used in split
* Update license years

## 5.6.0

* Releases multiple fixes listed in [milestone 1](https://github.com/cancerit/vafCorrect/milestone/1)

## 5.5.0

* Adds docker build - please use this in preference to other methods.  Should be functional with singularity.

## 5.4.1

* Remove help text about adding to allowed perl lib paths as not being used
  * corrected to indicate install skip for deps
* Changed so that following are installed even when 'skip' is used
  * htslib/Bio::DB::HTS
  * vcftools
* Removed hiding of errors during install
  * Some failures not surfaced so user needed to know to check `$?` or the setup.log file
* `Bio::DB::HTS` requires current bioperl.

## 5.4.0
* Reduced memory footprint for augment vcf file
* Augmentation can now run per chromosome
* SAM flag filters are now controlled by vafConstantas.pm file
* .bas file now used by default to calculate library size
* added support for cram files

## 5.3.9
* Change max reads read when determining read length to 50k
* Modify _get_bam_header_data to get the max read length of ALL bam files encountered and therefore to produce the correct lib_size
* Converted internal SAM flag constants to point to `Bio::DB::HTS` constants (in the future  internal constants could be replaced)
* Added command line options for filter include and exclude of reads when checking read length
* Added very simple test for `_get_bam_header_data`
* All contribute to fixing #45
* All hard tabs now spaces. Fixes whitespace part of #49

## 5.3.8

* Updated CHANAGES.md with corerct version

## 5.3.7

* corrected condition in vafcorrect for erroneous PASS flag setting

## 5.3.6

* Correct query_full to ensure `undef` return is caught.
* vaf output precison raised from two to four decimal places.

## 5.3.5

* Update `samtols` and `hstlib` to v1.7
* Update `Bio::DB::HTS` to v2.10
* Change tabix->query call to tabix->query_full
* fixed bug in accessing chr prefixed chromosome names for HDR regions

## 5.3.4

* corrected typo in setup.sh removed safelist branches from .travis
* Added checks for INST_METHOD
* Changed CHANGES.md syntax to markdown linter
* Exonerate install is now skipped if INST_METHOD is 2

## 5.3.2

* Fixed bug which caused incorrent assignment of mutant reads overlapping deletion regions to Unknon category - if there is any mismatch in 3prime region from deletionstart site

## 5.2.2

* Heade line for SAMPLE is now generated from command line
* Fixed bug for missing header data when running with bedonly option

## 5.2.1

* added condition to test vcf extension when bed file is not present

## 5.2.0

* added functionality to process bed file per chromosome
* progress file check during concatenation step ignores already processed chromosomes
* vcf file extension is optional if --bedOnly option is selected

## 5.1.3

* Fixed broken INFO tags in vcf and tsv
* tags -t option is now only applicable to tsv file

## 5.1.2

* forced default shell to bash for syste command

## 5.1.1

* Added custom sorting for vcf files

## 5.1.0

* User can now run per chromosome analysis [useful to run in parallel ] and concatenate the resulting files

## 5.0.0

* Removed exonerate parameter to filter alignments based on score
* Removed condition to discard alignment length smaller than read length
* Added user defined option to specify exonerat percent value
* Increased padding to extended reference sequence -- this will allow to map longer reads in amplicon data.

## 4.5.6

* Update to Bio::DB::HTS install to use already installed htslib

## 4.5.5

* updated version tag

## 4.5.4

* Removed default high depth regions filter option, users can optionally provide this file on command line (refer wiki)

## 4.5.3 * Updated licensing policy

* Added option to supply external ignore depth file
* removed legacy code to run external commands
* Added prerelease.sh to

## 4.5.2

* Updated .gitignore to allow test data to be added to git

## 4.5.1

* Updated test data set files

## 4.5.0

* Added option to restrict reads based on mapping quality threshold

## 4.4.2

* Exonerate --percent threshold reduced to 90 from 95 to allow 4 mismatches

## 4.4.1

* More robust check for mismatches at varinat positions in the alignment

## 4.4.0

* added additional condition to check mismatch at the variant region in the alignment
* reads mapping on alternate sequence and have mismatch at variant region were categorised into UNK [ unknown ] reads
* updated option to correctly read bed file

## 4.3.6

* Creates augmented vcf file even if no records were augmented
* Corrected processlog in vcf output header

## 4.3.5

* Minor change to test script

## 4.3.4

* Appropriate warning message when no vcf file found for a sample

## 4.3.3

* updated changes.md and verision information

## 4.3.2

* corrected tabix syntax as for Bio::DB::HTS::Tabix, added test for high depth region overlap

## 4.3.1

* removed legacy script from Makefile

## 4.3.0

* removed database dependecy code for external release

## 4.2.3

* replaced Tabix with Bio::DB::HTS::Tabix

## 4.2.2

* changed Bio::DB::HTS version to 1.12

## 4.2.1

* changed Bio::DB::HTS version to 1.11

## 4.2.0

* Added support to Bio::DB::HTS

## 4.1.17

* Updated SQL query to get input data to generate commands

## 4.1.16

* corrected condition to get unmatched normal

## 4.1.15

* Fixed vcf header, as it was copied twice in the outfile

## 4.1.14

* Added explicit check for normal sample using SIP_ATTRIBUTES table

## 4.1.13

* updated defualt project ini folder to cgpVaf

## 4.1.12

* Updated verison

## 4.1.11

* Added condition to ignore writing progress file if augment only option is selected

## 4.1.10

* Improved and corrected setup.sh and location of Makefile.PL

## 4.1.8

*added test to check final result removal
*chnaged bio:db:sam version to 1.42

## 4.1.7

*changed installation path to bin from perl/bin in setup.sh

## 4.1.6

*fixed test to created progress file path

## 4.1.5

*updated condition to create outfile only if it is absent

## 4.1.4

*updated file postfix to <variant>_vaf

## 4.1.3

*Added tumour sample name as postfix to tmp folder

## 4.1.2

*Added setup.sh

## 4.1.1

*Fixed bug in unmapped mate retrival

## 4.1.0

*Added option to accept vcf file as command line input

## 4.0.1

* Fixed readme file

## 4.0

* Contains pre-calculated VAF values in FORMAT field
* Per chromosome progress tracking i.e., failed jobs can be resumed from last unsuccessful chromosome
* Added setup script script

## 3.2.4

* Added condition to check snps that fall at the end of chromosome

## 3.2.3

* Added option to create config file for user defined list of samples
