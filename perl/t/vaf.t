#!/usr/bin/perl-5.16.3
use strict;
use Test::More;
use Test::Fatal;
use File::Spec;
use FindBin qw($Bin);
use File::Temp qw/ :seekable /;
use Data::Dumper;

use Sanger::CGP::Config::Config qw(vcfCommons);
use Sanger::CGP::Vaf::Data::ReadVcf;
use Sanger::CGP::Vaf::Process::Variant;
use Bio::DB::Sam;

use Config::IniFiles;
use Vcf;;
use English qw( -no_match_vars );
use warnings FATAL => 'all';
use Carp;
use Getopt::Long;
use Set::IntSpan;
use Math::Round;
use Pod::Usage;
use File::Path;
use Const::Fast qw(const);
use List::Util qw(first);
use Try::Tiny qw(try catch finally);
use Capture::Tiny qw(:all);

const my $MODULE1 => ' Sanger::CGP::Vaf::Data::ReadVcf';
const my $MODULE2 => 'Sanger::CGP::Vaf::Process::Variant';
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
};

my $options={
	'd'	=>	$test_data,
	'o'	=>	$test_output,
	'a'	=>	$test_variant_type1,
	'g' =>  $test_genome,
	'tn'=>  \@test_samples,
	'nn'=>  $normal_sample,
	'e'	=>	$test_ext2,
	't'=>"VT,VC",
	'r'=>0,
	'c'=>005
	};

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

subtest 'ReadVcf' => sub {
	my $vcf_obj = Sanger::CGP::Vaf::Data::ReadVcf->new($options);
	$vcf_obj->getAllSampleNames;
	#diag(@{$vcf_obj->{'allSamples'}}[2]);
	my($info_tag_val,$updated_info_tags,$vcf_file_obj)=$vcf_obj->getVcfHeaderData;
	my($variant,$bam_header_data,$bam_objects)=$vcf_obj->getVarinatObject($info_tag_val);
	diag(Dumper $variant);
	is_deeply($vcf_obj->getChromosomes,\@chr,'ReadVcf:getChromosomes');
	my($bed_locations)=$vcf_obj->getBedHash;
	my ($chromosomes)=$vcf_obj->getChromosomes;
			
	#my($progress_fhw,$progress_data)=$vcf_obj->getProgress;
};


done_testing();
