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

const my $MODULE => ' Sanger::CGP::Vaf::Data::ReadVcf';
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

my @test_samples=qw(PD20515a PD20515c);
const my $normal_sample=> 'PD20515b'; 
const my $test_bam1 => "$Bin/testData/PD20515a.bam";
const my $test_bam2 => "$Bin/testData/PD20515c.bam";
const my $test_bam3 => "$Bin/testData/PD20515b.bam";
const my $test_bed => "$Bin/testData/test.bed";

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
};

my $options={
	'd'	=>	$test_data,
	'o'	=>	$test_data,
	'a'	=>	$test_variant_type1,
	'g' =>  $test_genome,
	'tn'=>  "@test_samples",
	'nn'=>  $normal_sample,
	'e'	=>	$test_ext2,
	't'	=>	"VT,VC",
	'b'	=>	$test_bed,
	'ao' => 0,
	};


subtest 'vcf_object' => sub {
	my $vcf_obj = Sanger::CGP::Vaf::Data::ReadVcf->new($options);
	
};


done_testing();
