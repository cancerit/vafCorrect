
## v5.1.1 #######
* Added custom sorting for vcf files

## v5.1.0 #######
* User can now run per chromosome analysis [useful to run in parallel ] and concatenate the resulting files

## v5.0.0 #######
* Removed exonerate parameter to filter alignments based on score
* Removed condition to discard alignment length smaller than read length
* Added user defined option to specify exonerat percent value  
* Increased padding to extended reference sequence -- this will allow to map longer reads in amplicon data.

## v4.5.6 #######
* Update to Bio::DB::HTS install to use already installed htslib

## v4.5.5 #######
* updated version tag

## v4.5.4 #######
* Removed default high depth regions filter option, users can optionally provide this file on command line (refer wiki)

## v4.5.3 ########
* Updated licensing policy
* Added option to supply external ignore depth file
* removed legacy code to run external commands
* Added prerelease.sh to

## v4.5.2 ########

* Updated .gitignore to allow test data to be added to git

## v4.5.1 ########

* Updated test data set files

## v4.5.0 ########

* Added option to restrict reads based on mapping quality threshold

## v4.4.2 ########

* Exonerate --percent threshold reduced to 90 from 95 to allow 4 mismatches

## v4.4.1 ########

* More robust check for mismatches at varinat positions in the alignment

## v4.4.0 ########

* added additional condition to check mismatch at the variant region in the alignment
* reads mapping on alternate sequence and have mismatch at variant region were categorised into UNK [ unknown ] reads
* updated option to correctly read bed file

## v4.3.6 ########

* Creates augmented vcf file even if no records were augmented
* Corrected processlog in vcf output header

## v4.3.5 ########

* Minor change to test script

## v4.3.4 ########

* Appropriate warning message when no vcf file found for a sample

## v4.3.3 ########

* updated changes.md and verision information

## V4.3.2 ########

* corrected tabix syntax as for Bio::DB::HTS::Tabix, added test for high depth region overlap

## V4.3.1 ########

* removed legacy script from Makefile

## V4.3.0 ########

* removed database dependecy code for external release

## V4.2.3 ########

* replaced Tabix with Bio::DB::HTS::Tabix

## V4.2.2 ########

* changed Bio::DB::HTS version to 1.12

## V4.2.1 ########

* changed Bio::DB::HTS version to 1.11

## V4.2.0 ########

* Added support to Bio::DB::HTS

## V4.1.17 ########

* Updated SQL query to get input data to generate commands

## V4.1.16 ########

* corrected condition to get unmatched normal

## V4.1.15 ########

* Fixed vcf header, as it was copied twice in the outfile

## V4.1.14 ########

* Added explicit check for normal sample using SIP_ATTRIBUTES table

## V4.1.13 ########

* updated defualt project ini folder to cgpVaf

## V4.1.12 ########

* Updated verison

## V4.1.11 ########

* Added condition to ignore writing progress file if augment only option is selected

## V4.1.10 ########

* Improved and corrected setup.sh and location of Makefile.PL

## V4.1.8 ########

*added test to check final result removal
*chnaged bio:db:sam version to 1.42

## V4.1.7 ########

*changed installation path to bin from perl/bin in setup.sh

## V4.1.6 ########

*fixed test to created progress file path

## V4.1.5 ########

*updated condition to create outfile only if it is absent

## V4.1.4 ########

*updated file postfix to <varinat>_vaf

## V4.1.3 ########

*Added tumour sample name as postfix to tmp folder

## V4.1.2 ########

*Added setup.sh

## V4.1.1 ########

*Fixed bug in unmapped mate retrival

## V4.1.0 ########

*Added option to accept vcf file as command line input

## V4.0.1 ########

* Fixed readme file

## V4.0 ########

* Contains pre-calculated VAF values in FORMAT field
* Per chromosome progress tracking i.e., failed jobs can be resumed from last unsuccessful chromosome
* Added setup script script

### 3.2.4 ########

* Added condition to check snps that fall at the end of chromosome

### 3.2.3 ########
*	Added option to create config file for user defined list of samples
