#!/usr/bin/perl-5.16.3
use strict;
use Test::More;
use Test::Fatal;
use File::Spec;
use FindBin qw($Bin);
use File::Temp qw/ :seekable /;
use Data::Dumper;

#use Sanger::CGP::Config::Config qw(vcfCommons);
use Sanger::CGP::Vaf::Data::ReadVcf;
use Sanger::CGP::Vaf::Process::Variant;
use Sanger::CGP::Vaf::VafConstants;
use Bio::DB::HTS::Tabix;
use English qw( -no_match_vars );
use warnings FATAL => 'all';
use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use Const::Fast qw(const);
use File::Path qw(mkpath);

const my $MODULE1 => ' Sanger::CGP::Vaf::Data::ReadVcf';
const my $MODULE2 => 'Sanger::CGP::Vaf::Process::Variant';
const my $MODULE3 => 'Sanger::CGP::Vaf::VafConstants';
const my $test_data => "$Bin/testData";
const my $test_output => "$Bin/testData/test_output_vaf";
const my $test_project => '1086';
const my $test_ext1 => '.cave.annot.vcf.gz';
const my $test_ext2 => '.caveman_c.annot.vcf.gz';
const my $test_ext3 => '.pindel.annot.vcf.gz';
const my $test_variant_type1 => 'snp';
const my $test_variant_type2 => 'indel';
const my $test_genome => $test_data."/"."genome.fa";
const my $test_reads => $test_data."/"."temp.reads";

const my @test_samples => qw(PD20515a PD20515c);
const my @all_samples => qw(PD20515b PD20515a PD20515c);
const my $normal_sample => 'PD20515b'; 
const my $test_bam1 => "$Bin/testData/PD20515a.bam";
const my $test_bam2 => "$Bin/testData/PD20515c.bam";
const my $test_bam3 => "$Bin/testData/PD20515b.bam";
const my $test_bed => "$Bin/testData/test.bed";

subtest 'Initialisation checks' => sub {
  use_ok($MODULE1);
  use_ok($MODULE2);
  use_ok($MODULE3);
};

const my $tags => $Sanger::CGP::Vaf::VafConstants::SNP_TAGS;

mkpath("$test_output/tmpvcf_$test_samples[0]");

my $options={
	'd'	=>	$test_data,
	'o'	=>	$test_output,
	'a'	=>	$test_variant_type1,
	'g' =>  $test_genome,
	'tn'=>  \@test_samples,
	'nn'=>  $normal_sample,
	'e'	=>	$test_ext2,
	'ao' => 0,
	't'=>"VT,VC",
	'r'=>0,
	'c'=>005,
	'tmp' => "$test_output/tmpvcf_$test_samples[0]",
	#'b' => "$test_data/test.bed",
	#'bo' => 1
	};

my (%data_for_all_samples);
	my $info={'VT'=>'Sub','VC' =>'intronic'};
	
	$data_for_all_samples{'PD20515a'}{'1:16901544:C:T'}={'INFO'=>$info,'FILTER'=>'UM;MN;MQ;TI;HSD', 'RD'=>0};																									
	$data_for_all_samples{'PD20515a'}{'1:16901564:G:A'}={'INFO'=>$info,'FILTER'=>'UM;MN;TI;HSD', 'RD'=>0};
	$data_for_all_samples{'PD20515a'}{'1:16902712:T:C'}={'INFO'=>$info,'FILTER'=>'UM;MN;MQ', 'RD'=>0};
	$data_for_all_samples{'PD20515c'}{'1:16903781:C:T'}={'INFO'=>$info,'FILTER'=>'UM;MN;HSD', 'RD'=>0};
	$data_for_all_samples{'PD20515c'}{'1:16907525:G:C'}={'INFO'=>$info,'FILTER'=>'UM;MN;MQ', 'RD'=>0};
	$data_for_all_samples{'PD20515c'}{'1:2212488:A:G'}={'INFO'=>$info,'FILTER'=>'PASS', 'RD'=>0};
	#bed locations
	my $bed_info={'VT' => '-','VC' => '-'};
	#$data_for_all_samples{'PD20515a'}{'1:148021700:G:A'}={'INFO'=>$bed_info,'FILTER'=>'NA','RD'=>1};
	#$data_for_all_samples{'PD20515a'}{'1:148594340:G:A'}={'INFO'=>$bed_info,'FILTER'=>'NA','RD'=>1};
	#$data_for_all_samples{'PD20515c'}{'1:148021700:G:A'}={'INFO'=>$bed_info,'FILTER'=>'NA','RD'=>1};
	#$data_for_all_samples{'PD20515c'}{'1:148594340:G:A'}={'INFO'=>$bed_info,'FILTER'=>'NA','RD'=>1};
														
const my $normal_bam => $options->{'d'}.'/'.$normal_sample.'.bam';

const my @chr => qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT);

subtest 'AbstractVcf' => sub {
	my $vcf_obj = Sanger::CGP::Vaf::Data::ReadVcf->new($options);
	is_deeply($vcf_obj->getNormalName,$normal_sample,'AbstractVcf:getNormalName');
	is_deeply($vcf_obj->getTumourName,\@test_samples,'AbstractVcf:getTumourName');
	#is_deeply($vcf_obj->getAllSampleNames,\@all_samples,'AbstractVcf:getAllSampleNames');
	is_deeply($vcf_obj->getNormalBam,$normal_bam,'AbstractVcf:getNormalBam');
	undef $vcf_obj;	
};



my $expected_store_results;
my $data_for_all_samples_res;
my $unique_locations;
my $expected_unique_locations =  {
   '1:16902712:T:C' => 'PD20515a-UM;MN;MQ',
   '1:16901564:G:A' => 'PD20515a-UM;MN;TI;HSD',
   '1:16901544:C:T' => 'PD20515a-UM;MN;MQ;TI;HSD',
   '1:16903781:C:T' => 'PD20515c-UM;MN;HSD',
   '1:16907525:G:C' => 'PD20515c-UM;MN;MQ',
   '1:2212488:A:G' => 'PD20515c-PASS'
 };

subtest 'ReadVcf' => sub {
	my $vcf_obj = Sanger::CGP::Vaf::Data::ReadVcf->new($options);
	$vcf_obj->getAllSampleNames;
	#diag(@{$vcf_obj->{'allSamples'}}[2]);
	my($info_tag_val,$updated_info_tags,$vcf_file_obj)=$vcf_obj->getVcfHeaderData;
	my($bam_objects,$bas_files)=$vcf_obj->_get_bam_object;
	my($bam_header_data,$lib_size)=$vcf_obj->_get_bam_header_data($bam_objects,$bas_files);
	#create variant object
	my $variant=Sanger::CGP::Vaf::Process::Variant->new( 
		'location' 		=> undef,
		'varLine' 		=> undef,
		'varType' 		=> $vcf_obj->{'_a'},
		'libSize' 		=> defined $lib_size?$lib_size:undef,
		'samples' 		=> $vcf_obj->{'allSamples'},
		'tumourName'	=> $vcf_obj->getTumourName,
		'normalName'	=> $vcf_obj->getNormalName,
		'vcfStatus' 	=> $vcf_obj->{'vcf'},
		'noVcf'    		=> defined $vcf_obj->{'noVcf'}?$vcf_obj->{'noVcf'}:undef,
		'outDir'			=> $vcf_obj->getOutputDir,
		'passedOnly'  => $vcf_obj->{'_r'},
		'tmp'					=> $options->{'tmp'},
		'tabix_hdr' 		=> Bio::DB::HTS::Tabix->new(filename => "$Bin/hdr/seq.cov".$vcf_obj->{'_c'}.'.ONHG19_sorted.bed.gz')
		);

	#diag(Dumper $variant);
	is_deeply($vcf_obj->getChromosomes,\@chr,'ReadVcf:getChromosomes');
	my($bed_locations)=$vcf_obj->getBedHash;
	my ($chromosomes)=$vcf_obj->getChromosomes;	
	my($progress_fhw,$progress_data)=$vcf_obj->getProgress;
	if(!defined $options->{'bo'}){
		($data_for_all_samples_res,$unique_locations)=$vcf_obj->getMergedLocations('1',$updated_info_tags,$vcf_file_obj);
		is_deeply($unique_locations,$expected_unique_locations,'ReadVcf:getMergedLocations_unique_locations');
	}
	#diag(Dumper %data_for_all_samples);
	#test will not pass for bed only options
	#is_deeply($data_for_all_samples_res,\%data_for_all_samples,'ReadVcf:getMergedLocations');	
	if(defined $options->{'b'} ){
			($bed_locations)=$vcf_obj->filterBedLocations($unique_locations,$bed_locations);	
	}	
 	my $store_results;
 	($store_results)=$vcf_obj->processMergedLocations($data_for_all_samples_res,$unique_locations,$variant,$bam_header_data,$bam_objects,$store_results,'1',$tags,$info_tag_val,$progress_fhw,$progress_data);
		
	if(defined $bed_locations) {
		my($data_for_all_samples,$unique_locations)=$vcf_obj->populateBedLocations($bed_locations,$updated_info_tags);
		#diag(Dumper $data_for_all_samples,$unique_locations);
		($store_results)=$vcf_obj->processMergedLocations($data_for_all_samples,$unique_locations,$variant,$bam_header_data,$bam_objects,$store_results,'bed_file_data',$tags,$info_tag_val,$progress_fhw,$progress_data);	
		diag(Dumper $store_results);
	}
	# if augment option is not provided - results will go in tmp files (store results will be empty)...
	is_deeply($store_results,$expected_store_results,'ReadVcf:processMergedLocations');
	my($outfile_name_no_ext)=$vcf_obj->writeFinalFileHeaders($info_tag_val,$tags);
	$vcf_obj->catFiles($options->{'tmp'},'vcf',$outfile_name_no_ext);
	$vcf_obj->catFiles($options->{'tmp'},'tsv',$outfile_name_no_ext);
	my($outfile_gz,$outfile_tabix)=$vcf_obj->compressVcf("$outfile_name_no_ext.vcf");	
	
	
	# test 
	#my $expected_ref_seq=qw(CACCAGCAACTACCTCAGCCAGTCAGCTCCGTTCTACCTCTGTCATCTCAGATGAGAAGAGCAGGCCAGTATCTCTGGCCTTACCTGAAATATCTTAAGGCCGTAATTTACATTTTAGGCATGAATGATTTTCTAAAACCCACGATCAGAGTTTCTCTGGGAATCGGCGTCTGGCTTAGGAACACATTCATTTGTTTGACAAATACCTTCCCAAAACTATTTTAAAACACAGCTGCTGGGCGGGACGCAGTGGCTCACACCTGTAATCCCAGGACTTTGGGTGGCCGAGGCGGGTGAATCACTTGAGGTCAGCAGTTCAAGACCAGCTTGGCCAACATAGTGAAATCCTGTCTCTACTAAAAATACAAAAATTAGCCGGGTGTGGCAGTGCATGCCTATAATCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATCGTTTGAACCTGGGAGGCGGAGGTTGCAGTGAGCCGAGATTGCACCACTGCACTCCAGCCTGGGCGACAGAACAAGACTCTGTCTCAAAAAAGTAAATAAATAAATAAATAAATAAAGCTTCATATCAGCATTTCCTTTTTGGGAACTATACTATTCATCTGAATTAGCATATATATATATATGGGGCCGGACACAGCGGCTCACACCTGTAATCTCAAAACTTTGGAAGGCCAAAACAGGTGGTTCACCGGAGGTCAGGTGTTTTGAGACATGTCTGGCCAACGTGGTGAAACCCCATCTCTACTAAAAATACCAAAATTAGCCAGGCGTGGTGGTACGCCGCACCTGTAATCCCAGCTACTCAGGATGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCAAAGATTGCAGTGAGCCGAGATCACGCCATTGCACTCCAGCAGGGGTGACAGACTGAGACTCCATCTCAAAAAAGAAGTCTACCACATTTTACTCTGAGACAAGGAAATGTCCACAGGGAAGTGGCCACACACAGAAGTTAACCTAAAAGACAATGAATTCAGAGGACGGACATGAACAAATGTGCAATTTAAAACACAGGCCAGGTGCAGTGGCAACCCCTATAATCCCAGAGCTTTGGGAGGCCAAGGCGGGCTCATCACATGAGGTCAGGACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAGTACAAAAATTAGCCAGGCGTGGTGGCACATGCCTGAAATCCCAGCTACTCGGG);
	#my $expected_reconstructed_alt_seq=qw(CACCAGCAACTACCTCAGCCAGTCAGCTCCGTTCTACCTCTGTCATCTCAGATGAGAAGAGCAGGCCAGTATCTCTGGCCTTACCTGAAATATCTTAAGGCCGTAATTTACATTTTAGGCATGAATGATTTTCTAAAACCCACGATCAGAGTTTCTCTGGGAATCGGCGTCTGGCTTAGGAACACATTCATTTGTTTGACAAATACCTTCCCAAAACTATTTTAAAACACAGCTGCTGGGCGGGACGCAGTGGCTCACACCTGTAATCCCAGGACTTTGGGTGGCCGAGGCGGGTGAATCACTTGAGGTCAGCAGTTCAAGACCAGCTTGGCCAACATAGTGAAATCCTGTCTCTACTAAAAATACAAAAATTAGCCGGGTGTGGCAGTGCATGCCTATAATCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATCGTTTGAACCTGGGAGGCGGAGGTTGCAGTGAGCCGAGATTGCACCACTGCACTCCAGCCTGGGCGACAGAACAAGACTCTGTCTCAAAAAAGTAAATAAATAAATAAATAAATAAAGCTTCATATCAGCATTTCCTTTTTGGGAACTATACTATTCATCTGAATTAGCATATATATATATGGGGCCGGACACAGCGGCTCACACCTGTAATCTCAAAACTTTGGAAGGCCAAAACAGGTGGTTCACCGGAGGTCAGGTGTTTTGAGACATGTCTGGCCAACGTGGTGAAACCCCATCTCTACTAAAAATACCAAAATTAGCCAGGCGTGGTGGTACGCCGCACCTGTAATCCCAGCTACTCAGGATGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCAAAGATTGCAGTGAGCCGAGATCACGCCATTGCACTCCAGCAGGGGTGACAGACTGAGACTCCATCTCAAAAAAGAAGTCTACCACATTTTACTCTGAGACAAGGAAATGTCCACAGGGAAGTGGCCACACACAGAAGTTAACCTAAAAGACAATGAATTCAGAGGACGGACATGAACAAATGTGCAATTTAAAACACAGGCCAGGTGCAGTGGCAACCCCTATAATCCCAGAGCTTTGGGAGGCCAAGGCGGGCTCATCACATGAGGTCAGGACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAGTACAAAAATTAGCCAGGCGTGGTGGCACATGCCTGAAATCCCAGCTACTCGGG);
	#my $resulted_ref_seq = Sanger::CGP::VcfCompare::VcfMergeAndPileup::_get_dna_segment($expected_obj->{'PD20515c'},$g_pu_new->{'chr'},$g_pu_new->{'pos_5p'},$g_pu_new->{'pos_3p'});
	
};

# test individula subsroutines .....




 subtest 'ReadVcf_processMergedLocations' => sub {
 
 	my $expected_g_pu={
   'chr' => '1',
   'region' => '1:16902712-16902712',
   'end' => '16902712',
   'alt_seq' => 'C',
   'ins_flag' => undef,
   'del_flag' => undef,
   'ref_seq' => 'T',
   'start' => '16902712'
 };

  my $vcf_obj = Sanger::CGP::Vaf::Data::ReadVcf->new($options);
	$vcf_obj->getAllSampleNames;
	my($info_tag_val,$updated_info_tags,$vcf_file_obj)=$vcf_obj->getVcfHeaderData;
	my($bam_objects,$bas_files)=$vcf_obj->_get_bam_object;
	#lib size information only applicable for indels....
	my($bam_header_data,$lib_size)=$vcf_obj->_get_bam_header_data($bam_objects,$bas_files);
	#create variant object
	my $variant=Sanger::CGP::Vaf::Process::Variant->new( 
		'location' 		=> undef,
		'varLine' 		=> undef,
		'varType' 		=> $vcf_obj->{'_a'},
		'libSize' 		=> defined $lib_size?$lib_size:undef,
		'samples' 		=> $vcf_obj->{'allSamples'},
		'tumourName'	=> $vcf_obj->getTumourName,
		'normalName'	=> $vcf_obj->getNormalName,
		'vcfStatus' 	=> $vcf_obj->{'vcf'},
		'noVcf'    		=> defined $vcf_obj->{'noVcf'}?$vcf_obj->{'noVcf'}:undef,
		'outDir'			=> $vcf_obj->getOutputDir,
		'passedOnly'  => $vcf_obj->{'_r'},
		'tabix_hdr' 		=> Bio::DB::HTS::Tabix->new(filename => "$Bin/hdr/seq.cov".$vcf_obj->{'_c'}.'.ONHG19_sorted.bed.gz')
		);
	 		
 	  $variant->setLocation('1:16902712:T:C');
		$variant->setVarLine('PD20515a-UM;MN;MQ');
		my($g_pu)=$variant->formatVarinat();
		is_deeply($g_pu,$expected_g_pu,'ReadVcf_processMergedLocations:g_pu');
 		$g_pu=$variant->populateHash($g_pu,'PD20515a',$bam_header_data);
 		$g_pu=$variant->getPileup($bam_objects->{'PD20515a'},$g_pu);
 };


subtest 'CleanTestResults' => sub {
	  my $vcf_obj = Sanger::CGP::Vaf::Data::ReadVcf->new($options);
	  my($cleaned1) = $vcf_obj->cleanTempdir($options->{'tmp'}); 
		is_deeply($cleaned1,$options->{'tmp'},'CleanTestResults:tmpdir');
		my($cleaned2) = $vcf_obj->cleanTempdir($test_output);
		is_deeply($cleaned2,$test_output,'CleanTestResults:testOutdir');
};


done_testing();
