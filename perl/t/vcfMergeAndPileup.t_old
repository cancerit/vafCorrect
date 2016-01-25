#!/usr/bin/perl-5.16.3
use strict;
use Test::More;
use Test::Fatal;

use File::Spec;
use FindBin qw($Bin);
use File::Temp qw/ :seekable /;
use Data::Dumper;

$main::SQL_LIB_LOC = '.'; # this suppresses warnings about uninitialised values
use Sanger::CGP::Config::Config qw(vcfCommons);
use Sanger::CGP::VcfCompare::VcfMergeAndPileup;
use Sanger::CGP::WholeGenome::Base;
use Sanger::CGP::Database::Conn;

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

const my $SEP => "/";
const my $test_mode => 0;
const my $SPANNING_REGION => 1;
const my $MAX_INSERT_SIZE_FACTOR => 1;
const my $PROPER_PAIRED => 0x2;
const my $UNMAPPED => 0x4;
const my $REVERSE_STRAND => 0x10;
const my $MATE_REVERSE_STRAND => 0x20;
const my $NOT_PRIMARY_ALIGN => 0x100;
const my $MATE_UNMAPPED => 0x0008;
const my $READ_PAIRED => 0x0001;
const my $FIRST_IN_PAIR => 0x40;
const my $MAX_PILEUP_DEPTH => '1000000';
const my $SUPP_ALIGNMENT => 0x800;
const my $DUP_READ => 0x400;
const my $VENDER_FAIL => 0x200;

const my $MODULE => 'Sanger::CGP::VcfCompare::VcfMergeAndPileup';
const my $test_data => "$Bin/testData";
const my $test_project => '1086';
const my $test_ext1 => 'cave.annot.vcf.gz';
const my $test_ext2 => 'caveman_c.annot.vcf.gz';
const my $test_ext3 => 'pindel.annot.vcf.gz';
const my $test_variant_type1 => 'snp';
const my $test_variant_type2 => 'indel';
const my $test_config_file => join('/',$test_data,'1086_VcfMergeConfig.ini');
const my $test_genome => $test_data."/"."genome.fa";
const my $test_reads => $test_data."/"."temp.reads";

my @test_samples=qw(PD20515a PD20515c);
const my $normal_sample=> 'PD20515b';
const my $test_bam1 => "$Bin/testData/PD20515a.bam";
const my $test_bam2 => "$Bin/testData/PD20515c.bam";
const my $test_bam3 => "$Bin/testData/PD20515b.bam";
const my $test_bed => "$Bin/testData/test.bed";

my $base = Sanger::CGP::WholeGenome::Base->new;
$base->db_type('live');
my $conn = $base->connection;

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
 };

my $options={
	'd'	=>	$test_data,
	'i'	=>	$test_config_file,
	'o'	=>	$test_data,
	'a'	=>	$test_variant_type1,
	'e'	=>	$test_ext2,
	't'	=>	"VT,VC",
	'b'	=>	$test_bed,
	'ao' => 0,
	's' => 1,
	
	};

	my $expected_output_options={
	'd'=>$test_data,
	'i'=>$test_config_file,
	'o' =>"$test_data/output",
	'a'=>$test_variant_type1,
	'e'=>$test_ext2,
	't'=>"VT,VC",
	'g'=>$test_genome,
	'b'	=>	$test_bed,
	'c' =>  '005',
	'r' => '1',
	'ao' => 0,
	's' => 1,
	'bo' => 0
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
	$data_for_all_samples{'PD20515a'}{'1:148021700:G:A'}={'INFO'=>$bed_info,'FILTER'=>'NA','RD'=>1};
	$data_for_all_samples{'PD20515a'}{'1:148594340:G:A'}={'INFO'=>$bed_info,'FILTER'=>'NA','RD'=>1};
	$data_for_all_samples{'PD20515c'}{'1:148021700:G:A'}={'INFO'=>$bed_info,'FILTER'=>'NA','RD'=>1};
	$data_for_all_samples{'PD20515c'}{'1:148594340:G:A'}={'INFO'=>$bed_info,'FILTER'=>'NA','RD'=>1};
																							

  my $expected_obj;
	foreach my $sample(@test_samples) {
	my $expected_sam = Bio::DB::Sam->new(-bam=>"$Bin/testData/$sample.bam", 
																				-fasta =>$test_genome);
	$expected_sam->max_pileup_cnt($MAX_PILEUP_DEPTH);
	$expected_obj->{$sample}=$expected_sam;
	
	}
	

subtest 'get_options' => sub {
	my $resulting_output=Sanger::CGP::VcfCompare::VcfMergeAndPileup::get_options($options);
	is_deeply($resulting_output,$expected_output_options,'get_options:test');
};

subtest 'unique_locations' => sub {

	my $expected_output1={ '1:16901544:C:T' => 'PD20515a-UM;MN;MQ;TI;HSD', '1:16901564:G:A' => 'PD20515a-UM;MN;TI;HSD', '1:16902712:T:C' => 'PD20515a-UM;MN;MQ',
												'1:16903781:C:T' => 'PD20515c-UM;MN;HSD', '1:16907525:G:C' => 'PD20515c-UM;MN;MQ','1:2212488:A:G' => 'PD20515c-PASS',
												 '1:148021700:G:A' => "$test_bed-BEDFILE", '1:148594340:G:A' => "$test_bed-BEDFILE"
												};
																				
	my $expected_output5={'PD20515b'=>'|PD20515a|PD20515c'};
	
	my ($resulting_output1,$resulting_output2,$resulting_output3,$resulting_output4,$resulting_output5)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::get_unique_locations(\@test_samples,$expected_output_options);
	is_deeply($resulting_output1,$expected_output1,'unique_locations:unique_locations_hash');
	is_deeply($resulting_output2,\%data_for_all_samples,'unique_locations:all_locations_hash');
	is_deeply($resulting_output4,1,'unique_locations:vcf_flag');
	is_deeply($resulting_output5,$expected_output5,'unique_locations:normal sample');

};

subtest 'run_and_consolidate_mpileup_SNP' => sub {
	my $normal_sample_hash={'PD20515b'=>1};
	my $expected_normal='PD20515b';
	my $expected_flag=1;
	#Sanger::CGP::VcfCompare::VcfMergeAndPileup::load_sql($conn);
	#my ($resulting_output,$normal_flag)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::_create_symlink_for_normal($expected_output_options->{'d'},$normal_sample_hash,$conn);
	#is_deeply($resulting_output,$expected_normal,'run_and_consolidate_mpileup: _create_symlik_for_normal');
	#is_deeply($normal_flag,$expected_flag,'run_and_consolidate_mpileup: _create_symlik_for_normal');
	#diag("------>$resulting_output");
 # different bam objects returned ???
	
	#my ($resulting_bam_obj)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::_get_bam_object(\@test_samples,$expected_output_options->{'d'},$test_genome,);
	#is_deeply($resulting_bam_obj,$expected_obj,'run_and_consolidate_mpileup: _get_bam_object');

	unshift(@test_samples,'PD20515b');
	my $test_location='1:16901544:C:T';
	
	my $expected_info={'NA' => 2, 'NS' => 2, 'NC' => 1, 'NP' => 0, 'VT'=>'Sub','VC' =>'intronic',  };
	my $exp_NFS=2;
	my $exp_OFS_val='1';
	my $exp_original_flag='UM;MN;MQ;TI;HSD';
	
	# test _get_original_results
	my $vcf_file_status={'PD20515a' => 1, 'PD20515c' => 1 ,'PD20515b' => 0};	
	my $no_vcf_counter=0;
	my $expexted_depth=0;
	my ($original_vcf_info, $NFS,$original_flag,$max_depth)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::_get_original_results(\@test_samples, \%data_for_all_samples, $test_location,$vcf_file_status,$no_vcf_counter,$expected_output_options,$test_location);
	is_deeply($original_vcf_info,$expected_info,'run_and_consolidate_mpileup:_get_original_results:original_vcf_info');
	is_deeply($NFS,$exp_NFS,'run_and_consolidate_mpileup:_get_original_results:NFS');
	is_deeply($max_depth,$expexted_depth,'run_and_consolidate_mpileup:_get_original_results:OFS_val');
	is_deeply($original_flag->{'PD20515a'},$exp_original_flag,'run_and_consolidate_mpileup:_get_original_results:original_flag');
	
	# test _get_region
	my $exp_chr='1';
	my $exp_start='16901544';
	my $exp_end='16901544';
	my $exp_region='1:16901544-16901544';
	my $exp_alt='T';
	my $exp_ref='C';
	my $exp_insertion_flag=undef; # for snp data
	my $exp_del_flag=undef;	# for snp data
	
	my $exp_output={'chr'=>$exp_chr, 'start' => $exp_start, 'end' => $exp_end,
				'region' => $exp_region, 'ins_flag' => $exp_insertion_flag,
				'del_flag' => $exp_del_flag, 'alt_seq' => $exp_alt, 'ref_seq' => $exp_ref
			};

	my ($g_pu)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::_get_region($test_location,$test_variant_type1);
	is_deeply($g_pu,$exp_output,'run_and_consolidate_mpileup:_get_region');
	
	
	my $exp_g_pu = {'chr' => 1,
                'start' => '2212488',
                'end' => '2212488',
                'region'=> '1:2212488-2212488',
                'depth' => 0,    # populated later in format pileup line sub
                'ref_p' => 1,
                'ref_n' => 1,
                'alt_p' => 3,
                'alt_n' => 2                         
              };
	
	my $initialized_g_pu = {'chr' => 1,
                'start' => '2212488',
                'end' => '2212488',
                'region'=> '1:2212488-2212488',
                'depth' => 0,
                'ref_p' => 0,
                'ref_n' => 0,
                'alt_p' => 0,
                'alt_n' => 0                         
              };

	
	
	my ($resulted_g_pu)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::_get_pileup($expected_obj->{'PD20515c'},$initialized_g_pu);
	is_deeply($resulted_g_pu,$exp_g_pu,'run_and_consolidate_mpileup:_get_pileup');
	
	};
	
subtest 'run_and_consolidate_mpileup_indels' => sub {
	my $bam_header_data;
	my %expected_read_length;
	my $vcf_file_status={'PD20515a' => 1, 'PD20515c' => 1 ,'PD20515b' => 0};	
	
	my $populate_hash;
	my $bam_header->{'PD1234_test'}={'read_length'=>75, 
																		'lib_size'=>370
																	};

	my $expected_hash={
								'ref_p' => 0,
                'ref_n' => 0,
                'alt_p' => 0,
                'alt_n' => 0,
                'read_length' => 75,
                'lib_size'=>370,
                'sample' => 'PD1234_test',
                'amb' => '0'
                };
	my($resulted_hash)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::_populate_hash($populate_hash,'PD1234_test',$bam_header);
	is_deeply($resulted_hash,$expected_hash,'run_and_consolidate_mpileup:_populate_hash');
	


	my $g_pu_new={
								'chr' => 1,
                'start' => '1492875',
                'end' => '1492878',
                'pos_5p' => '1492273',
                'pos_3p' => '1493477',
                'ins_flag' =>0,
                'del_flag' =>1,
                'alt_seq' => 'C',
                'ref_p' => 0,
                'ref_n' => 0,
                'alt_p' => 0,
                'alt_n' => 0,
                'hdr' => 0,
                'lib_size' => 370,
                'region' => '1:1492875-1492878',
                'read_length' => 75,
                'E' => 0
              };
              
  $g_pu_new->{'ref_pos_5p'}=($g_pu_new->{'start'} - $g_pu_new->{'pos_5p'}) + 1 ;
  $g_pu_new->{'ref_pos_3p'}=$g_pu_new->{'ref_pos_5p'} - ( $g_pu_new->{'end'} - $g_pu_new->{'start'} ) -1 ;
  $g_pu_new->{'alt_pos_3p'}=$g_pu_new->{'ref_pos_5p'} - length( $g_pu_new->{'alt_seq'}) -1 ;
  
	my $expected_ref_seq=qw(CACCAGCAACTACCTCAGCCAGTCAGCTCCGTTCTACCTCTGTCATCTCAGATGAGAAGAGCAGGCCAGTATCTCTGGCCTTACCTGAAATATCTTAAGGCCGTAATTTACATTTTAGGCATGAATGATTTTCTAAAACCCACGATCAGAGTTTCTCTGGGAATCGGCGTCTGGCTTAGGAACACATTCATTTGTTTGACAAATACCTTCCCAAAACTATTTTAAAACACAGCTGCTGGGCGGGACGCAGTGGCTCACACCTGTAATCCCAGGACTTTGGGTGGCCGAGGCGGGTGAATCACTTGAGGTCAGCAGTTCAAGACCAGCTTGGCCAACATAGTGAAATCCTGTCTCTACTAAAAATACAAAAATTAGCCGGGTGTGGCAGTGCATGCCTATAATCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATCGTTTGAACCTGGGAGGCGGAGGTTGCAGTGAGCCGAGATTGCACCACTGCACTCCAGCCTGGGCGACAGAACAAGACTCTGTCTCAAAAAAGTAAATAAATAAATAAATAAATAAAGCTTCATATCAGCATTTCCTTTTTGGGAACTATACTATTCATCTGAATTAGCATATATATATATATGGGGCCGGACACAGCGGCTCACACCTGTAATCTCAAAACTTTGGAAGGCCAAAACAGGTGGTTCACCGGAGGTCAGGTGTTTTGAGACATGTCTGGCCAACGTGGTGAAACCCCATCTCTACTAAAAATACCAAAATTAGCCAGGCGTGGTGGTACGCCGCACCTGTAATCCCAGCTACTCAGGATGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCAAAGATTGCAGTGAGCCGAGATCACGCCATTGCACTCCAGCAGGGGTGACAGACTGAGACTCCATCTCAAAAAAGAAGTCTACCACATTTTACTCTGAGACAAGGAAATGTCCACAGGGAAGTGGCCACACACAGAAGTTAACCTAAAAGACAATGAATTCAGAGGACGGACATGAACAAATGTGCAATTTAAAACACAGGCCAGGTGCAGTGGCAACCCCTATAATCCCAGAGCTTTGGGAGGCCAAGGCGGGCTCATCACATGAGGTCAGGACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAGTACAAAAATTAGCCAGGCGTGGTGGCACATGCCTGAAATCCCAGCTACTCGGG);
	my $expected_reconstructed_alt_seq=qw(CACCAGCAACTACCTCAGCCAGTCAGCTCCGTTCTACCTCTGTCATCTCAGATGAGAAGAGCAGGCCAGTATCTCTGGCCTTACCTGAAATATCTTAAGGCCGTAATTTACATTTTAGGCATGAATGATTTTCTAAAACCCACGATCAGAGTTTCTCTGGGAATCGGCGTCTGGCTTAGGAACACATTCATTTGTTTGACAAATACCTTCCCAAAACTATTTTAAAACACAGCTGCTGGGCGGGACGCAGTGGCTCACACCTGTAATCCCAGGACTTTGGGTGGCCGAGGCGGGTGAATCACTTGAGGTCAGCAGTTCAAGACCAGCTTGGCCAACATAGTGAAATCCTGTCTCTACTAAAAATACAAAAATTAGCCGGGTGTGGCAGTGCATGCCTATAATCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATCGTTTGAACCTGGGAGGCGGAGGTTGCAGTGAGCCGAGATTGCACCACTGCACTCCAGCCTGGGCGACAGAACAAGACTCTGTCTCAAAAAAGTAAATAAATAAATAAATAAATAAAGCTTCATATCAGCATTTCCTTTTTGGGAACTATACTATTCATCTGAATTAGCATATATATATATGGGGCCGGACACAGCGGCTCACACCTGTAATCTCAAAACTTTGGAAGGCCAAAACAGGTGGTTCACCGGAGGTCAGGTGTTTTGAGACATGTCTGGCCAACGTGGTGAAACCCCATCTCTACTAAAAATACCAAAATTAGCCAGGCGTGGTGGTACGCCGCACCTGTAATCCCAGCTACTCAGGATGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCAAAGATTGCAGTGAGCCGAGATCACGCCATTGCACTCCAGCAGGGGTGACAGACTGAGACTCCATCTCAAAAAAGAAGTCTACCACATTTTACTCTGAGACAAGGAAATGTCCACAGGGAAGTGGCCACACACAGAAGTTAACCTAAAAGACAATGAATTCAGAGGACGGACATGAACAAATGTGCAATTTAAAACACAGGCCAGGTGCAGTGGCAACCCCTATAATCCCAGAGCTTTGGGAGGCCAAGGCGGGCTCATCACATGAGGTCAGGACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAGTACAAAAATTAGCCAGGCGTGGTGGCACATGCCTGAAATCCCAGCTACTCGGG);
	my $resulted_ref_seq = Sanger::CGP::VcfCompare::VcfMergeAndPileup::_get_dna_segment($expected_obj->{'PD20515c'},$g_pu_new->{'chr'},$g_pu_new->{'pos_5p'},$g_pu_new->{'pos_3p'});
	is_deeply($resulted_ref_seq,$expected_ref_seq,'run_and_consolidate_mpileup_indels:_get_dna_segment_ref');
	
	my $resulted_reconstructed_alt_seq =Sanger::CGP::VcfCompare::VcfMergeAndPileup::_get_alt_seq($expected_obj->{'PD20515c'},$g_pu_new);
	is_deeply($resulted_reconstructed_alt_seq,$expected_reconstructed_alt_seq,'run_and_consolidate_mpileup_indels:_get_dna_segment_alt');

  open(my $Read_FH_new,">$options->{'o'}/temp.reads");

	Sanger::CGP::VcfCompare::VcfMergeAndPileup::_fetch_features($expected_obj->{'PD20515c'}, $g_pu_new,$Read_FH_new);


			my $ref_seq_file="$options->{'o'}/temp.ref"; 
			open (my $ref_FH,'>'.$ref_seq_file);
			print $ref_FH ">ref\n$expected_ref_seq\n>alt\n$expected_reconstructed_alt_seq\n";
			close($ref_FH);

 #my($exonearte_results)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::_get_exonerate_output($ref_seq_file,"$options->{'o'}/temp.reads");
  
 my $expected_g_pu_new={
								'chr' => 1,
                'start' => '1492875',
                'end' => '1492878',
                'pos_5p' => '1492273',
                'pos_3p' => '1493477',
                'ins_flag' =>0,
                'del_flag' =>1,
                'alt_seq' => 'C',
                'ref_p' => 0,
                'ref_n' => 0,
                'alt_p' => 1,
                'alt_n' => 1,
                'hdr' => 0,
                'region' => '1:1492875-1492878',
                'lib_size' => 370,
                'read_length' => 75,
                'long_indel' => 0,
                'E' => 0
              };
              
  $expected_g_pu_new->{'ref_pos_5p'}=($expected_g_pu_new->{'start'} - $expected_g_pu_new->{'pos_5p'}) + 1 ;
  $expected_g_pu_new->{'ref_pos_3p'}=$expected_g_pu_new->{'ref_pos_5p'} - ( $expected_g_pu_new->{'end'} - $expected_g_pu_new->{'start'} ) -1 ;
  $expected_g_pu_new->{'alt_pos_3p'}=$expected_g_pu_new->{'ref_pos_5p'} - length( $expected_g_pu_new->{'alt_seq'}) -1 ;

              
 my ($resulted_g_pu)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::_do_exonerate($ref_seq_file,"$options->{'o'}/temp.reads",$g_pu_new);
is_deeply($resulted_g_pu,$expected_g_pu_new,'run_and_consolidate_mpileup:_do_exonerate');

};


subtest 'run_and_consolidate_mpileup:common' => sub {

my $original_flag={'PD1234_test' => 'PASS'};
my $g_pu_format={
								'ref_p' => 2,
                'ref_n' => 3,
                'alt_p' => 8,
                'alt_n' => 2,
                'sample'=>'PD1234_test',
                'amb' => '0'
              };
my $expected_pileup={'MTR'=>10,
										'WTR'=>5,
										'DEP'=>15,
										'MDR'=>3,
										'WDR'=>3,
										'OFS'=>'PASS',
										};

	my($resulted_g_pu)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::_format_pileup_line($original_flag,$g_pu_format,$expected_output_options);
	is_deeply($resulted_g_pu, $expected_pileup,'run_and_consolidate_mpileup:common_format_pileup_line');

};


sub print_hash {
	my $hash=shift;
	foreach my $key (sort keys %$hash) {
		if(defined($hash->{$key}))
		{	
			diag("$key---->$hash->{$key}---\n");	
		}
	}
}


 
done_testing();
