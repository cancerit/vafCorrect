
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
