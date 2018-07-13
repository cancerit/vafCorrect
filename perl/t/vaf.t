use strict;
use Test::More;
use Test::Fatal;
use File::Spec;
use FindBin qw($Bin);
use File::Temp qw/ :seekable /;
use Data::Dumper;
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

const my @MODULES => qw( Sanger::CGP::Vaf::Data::ReadVcf
                        Sanger::CGP::Vaf::Process::Variant
                        Sanger::CGP::Vaf::VafConstants);
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

const my @test_samples => qw(samplea samplec);
const my @all_samples => qw(sampleb samplea samplec);
const my $normal_sample => 'sampleb'; 
const my $test_bam1 => "$Bin/testData/samplea.bam";
const my $test_bam2 => "$Bin/testData/samplec.bam";
const my $test_bam3 => "$Bin/testData/sampleb.bam";
const my $test_bed => "$Bin/testData/test.bed";
const my $dummy_chr => 'bed_file_data';
const my $test_chr => '1';

subtest 'Initialisation checks' => sub {
    foreach my $mod(@MODULES){
        use_ok($mod);
    }
};

const my $tags => $Sanger::CGP::Vaf::VafConstants::SNP_TAGS;

mkpath("$test_output/tmpvcf_$test_samples[0]");

system("gunzip -c $test_genome.gz >$test_genome");

my $options={
    'd'    =>    $test_data,
    'o'    =>    $test_output,
    'a'    =>    $test_variant_type1,
    'g' =>  $test_genome,
    'tn'=>  \@test_samples,
    'nn'=>  $normal_sample,
    'e'    =>    $test_ext2,
    'ao' => 0,
    't'=>"VT,VC",
    'r'=>0,
    'c'=>005,
    'be' => ".bam",
    'tmp' => "$test_output/tmpvcf_$test_samples[0]",
    'finc' => $Sanger::CGP::Vaf::VafConstants::DEFAULT_READLEN_INCLUDE,
    'fexc' => $Sanger::CGP::Vaf::VafConstants::DEFAULT_READLEN_EXCLUDE,
    #'b' => "$test_data/test.bed",
    #'bo' => 1
    };

my (%data_for_all_samples);
    my $info={'VT'=>'Sub','VC' =>'intronic'};
    
    $data_for_all_samples{'samplea'}{'1:16901544:C:T'}={'INFO'=>$info,'FILTER'=>'UM;MN;MQ;TI;HSD', 'RD'=>0};
    $data_for_all_samples{'samplea'}{'1:16901564:G:A'}={'INFO'=>$info,'FILTER'=>'UM;MN;TI;HSD', 'RD'=>0};
    $data_for_all_samples{'samplea'}{'1:16902712:T:C'}={'INFO'=>$info,'FILTER'=>'UM;MN;MQ', 'RD'=>0};
    $data_for_all_samples{'samplec'}{'1:16903781:C:T'}={'INFO'=>$info,'FILTER'=>'UM;MN;HSD', 'RD'=>0};
    $data_for_all_samples{'samplec'}{'1:16907525:G:C'}={'INFO'=>$info,'FILTER'=>'UM;MN;MQ', 'RD'=>0};
    $data_for_all_samples{'samplec'}{'1:2212488:A:G'}={'INFO'=>$info,'FILTER'=>'PASS', 'RD'=>0};
    #bed locations
    my $bed_info={'VT' => '-','VC' => '-'};
                                                            
const my $normal_bam => $options->{'d'}.'/'.$normal_sample.'.bam';

const my @chr => qw(1);

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
   '1:16902712:T:C' => 'samplea-UM;MN;MQ',
   '1:16901564:G:A' => 'samplea-UM;MN;TI;HSD',
   '1:16901544:C:T' => 'samplea-UM;MN;MQ;TI;HSD',
   '1:16903781:C:T' => 'samplec-UM;MN;HSD',
   '1:16907525:G:C' => 'samplec-UM;MN;MQ',jw32@sanger.ac.uk
   '1:2212488:A:G' => 'samplec-PASS'
 };

subtest '_get_bam_header_data' => sub {
    $options->{'a'} = $test_variant_type2;
    my $vcf_obj = Sanger::CGP::Vaf::Data::ReadVcf->new($options);
    $vcf_obj->getAllSampleNames;
    my($info_tag_val,$vcf_file_obj)=$vcf_obj->getVcfHeaderData;
    my($bam_objects,$bas_files)=$vcf_obj->_get_bam_object;
    my($bam_header_data,$lib_size)=$vcf_obj->_get_bam_header_data($bam_objects,$bas_files);
    is($lib_size,401);
    $options->{'a'} = $test_variant_type1;
};

subtest 'ReadVcf' => sub {
    my $vcf_obj = Sanger::CGP::Vaf::Data::ReadVcf->new($options);
    $vcf_obj->getAllSampleNames;
    #diag(@{$vcf_obj->{'allSamples'}}[2]);
    my($info_tag_val,$vcf_file_obj)=$vcf_obj->getVcfHeaderData;
    my($bam_objects,$bas_files)=$vcf_obj->_get_bam_object;
    my($bam_header_data,$lib_size)=$vcf_obj->_get_bam_header_data($bam_objects,$bas_files);
    #create variant object
    my $variant=Sanger::CGP::Vaf::Process::Variant->new( 
        'location'         => undef,
        'varLine'         => undef,
        'varType'         => $vcf_obj->{'_a'},
        'libSize'         => defined $lib_size?$lib_size:undef,
        'samples'         => $vcf_obj->{'allSamples'},
        'tumourName'    => $vcf_obj->getTumourName,
        'normalName'    => $vcf_obj->getNormalName,
        'vcfStatus'     => $vcf_obj->{'vcf'},
        'noVcf'            => defined $vcf_obj->{'noVcf'}?$vcf_obj->{'noVcf'}:undef,
        'outDir'            => $vcf_obj->getOutputDir,
        'passedOnly'  => $vcf_obj->{'_r'},
        'tmp'                    => $options->{'tmp'},
    'tabix_hdr'   => defined  $vcf_obj->{'_hdr'}?Bio::DB::HTS::Tabix->new(filename => $vcf_obj->{'_hdr'}):undef,
        
        );

    #diag(Dumper $variant);
    is_deeply($vcf_obj->getChromosomes,\@chr,'ReadVcf:getChromosomes');

    my ($chromosomes)=$vcf_obj->getChromosomes([$test_chr]);    
    my ($progress_hash)=$vcf_obj->getProgress($chromosomes);
    my ($data_for_all_samples_res,$unique_locations)=$vcf_obj->getMergedLocations($test_chr,$vcf_file_obj);
    is_deeply($unique_locations,$expected_unique_locations,'ReadVcf:getMergedLocations_unique_locations');
    #diag(Dumper %data_for_all_samples);
    #test will not pass for bed only options
    #is_deeply($data_for_all_samples_res,\%data_for_all_samples,'ReadVcf:getMergedLocations');    
    if(defined $options->{'b'} ){
        my($bed_locations)=$vcf_obj->getBedHash($test_chr);
        if( $options->{'bo'} == 1 && (defined $data_for_all_samples_res) ) {
                ($data_for_all_samples_res, $unique_locations)=$vcf_obj->filterBedLocationsFromVCF($data_for_all_samples_res, $unique_locations, $bed_locations);    
            }else{
                 ($data_for_all_samples_res,$unique_locations)=$vcf_obj->populateBedLocations($data_for_all_samples_res,$unique_locations,$bed_locations);
            }
    }
        
     my $store_results;
     my($progress_fhw,$progress_data)=@{$progress_hash->{$test_chr}};
     ($store_results)=$vcf_obj->processMergedLocations($data_for_all_samples_res,$unique_locations,$variant,$bam_header_data,$bam_objects,$store_results,$test_chr,$tags,$info_tag_val,$progress_fhw,$progress_data);
    # if augment option is not provided - results will go in tmp files (store results will be empty)...
    is_deeply($store_results,$expected_store_results,'ReadVcf:processMergedLocations');
    my($outfile_name_no_ext)=$vcf_obj->writeFinalFileHeaders($info_tag_val,$tags);
    $vcf_obj->catFiles($options->{'tmp'},'vcf',$outfile_name_no_ext);
    $vcf_obj->catFiles($options->{'tmp'},'tsv',$outfile_name_no_ext);
    my($outfile_gz,$outfile_tabix)=$vcf_obj->gzipAndIndexVcf("$outfile_name_no_ext.vcf");    
    
    
    # test 
    #my $expected_ref_seq=qw(CACCAGCAACTACCTCAGCCAGTCAGCTCCGTTCTACCTCTGTCATCTCAGATGAGAAGAGCAGGCCAGTATCTCTGGCCTTACCTGAAATATCTTAAGGCCGTAATTTACATTTTAGGCATGAATGATTTTCTAAAACCCACGATCAGAGTTTCTCTGGGAATCGGCGTCTGGCTTAGGAACACATTCATTTGTTTGACAAATACCTTCCCAAAACTATTTTAAAACACAGCTGCTGGGCGGGACGCAGTGGCTCACACCTGTAATCCCAGGACTTTGGGTGGCCGAGGCGGGTGAATCACTTGAGGTCAGCAGTTCAAGACCAGCTTGGCCAACATAGTGAAATCCTGTCTCTACTAAAAATACAAAAATTAGCCGGGTGTGGCAGTGCATGCCTATAATCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATCGTTTGAACCTGGGAGGCGGAGGTTGCAGTGAGCCGAGATTGCACCACTGCACTCCAGCCTGGGCGACAGAACAAGACTCTGTCTCAAAAAAGTAAATAAATAAATAAATAAATAAAGCTTCATATCAGCATTTCCTTTTTGGGAACTATACTATTCATCTGAATTAGCATATATATATATATGGGGCCGGACACAGCGGCTCACACCTGTAATCTCAAAACTTTGGAAGGCCAAAACAGGTGGTTCACCGGAGGTCAGGTGTTTTGAGACATGTCTGGCCAACGTGGTGAAACCCCATCTCTACTAAAAATACCAAAATTAGCCAGGCGTGGTGGTACGCCGCACCTGTAATCCCAGCTACTCAGGATGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCAAAGATTGCAGTGAGCCGAGATCACGCCATTGCACTCCAGCAGGGGTGACAGACTGAGACTCCATCTCAAAAAAGAAGTCTACCACATTTTACTCTGAGACAAGGAAATGTCCACAGGGAAGTGGCCACACACAGAAGTTAACCTAAAAGACAATGAATTCAGAGGACGGACATGAACAAATGTGCAATTTAAAACACAGGCCAGGTGCAGTGGCAACCCCTATAATCCCAGAGCTTTGGGAGGCCAAGGCGGGCTCATCACATGAGGTCAGGACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAGTACAAAAATTAGCCAGGCGTGGTGGCACATGCCTGAAATCCCAGCTACTCGGG);
    #my $expected_reconstructed_alt_seq=qw(CACCAGCAACTACCTCAGCCAGTCAGCTCCGTTCTACCTCTGTCATCTCAGATGAGAAGAGCAGGCCAGTATCTCTGGCCTTACCTGAAATATCTTAAGGCCGTAATTTACATTTTAGGCATGAATGATTTTCTAAAACCCACGATCAGAGTTTCTCTGGGAATCGGCGTCTGGCTTAGGAACACATTCATTTGTTTGACAAATACCTTCCCAAAACTATTTTAAAACACAGCTGCTGGGCGGGACGCAGTGGCTCACACCTGTAATCCCAGGACTTTGGGTGGCCGAGGCGGGTGAATCACTTGAGGTCAGCAGTTCAAGACCAGCTTGGCCAACATAGTGAAATCCTGTCTCTACTAAAAATACAAAAATTAGCCGGGTGTGGCAGTGCATGCCTATAATCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATCGTTTGAACCTGGGAGGCGGAGGTTGCAGTGAGCCGAGATTGCACCACTGCACTCCAGCCTGGGCGACAGAACAAGACTCTGTCTCAAAAAAGTAAATAAATAAATAAATAAATAAAGCTTCATATCAGCATTTCCTTTTTGGGAACTATACTATTCATCTGAATTAGCATATATATATATGGGGCCGGACACAGCGGCTCACACCTGTAATCTCAAAACTTTGGAAGGCCAAAACAGGTGGTTCACCGGAGGTCAGGTGTTTTGAGACATGTCTGGCCAACGTGGTGAAACCCCATCTCTACTAAAAATACCAAAATTAGCCAGGCGTGGTGGTACGCCGCACCTGTAATCCCAGCTACTCAGGATGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCAAAGATTGCAGTGAGCCGAGATCACGCCATTGCACTCCAGCAGGGGTGACAGACTGAGACTCCATCTCAAAAAAGAAGTCTACCACATTTTACTCTGAGACAAGGAAATGTCCACAGGGAAGTGGCCACACACAGAAGTTAACCTAAAAGACAATGAATTCAGAGGACGGACATGAACAAATGTGCAATTTAAAACACAGGCCAGGTGCAGTGGCAACCCCTATAATCCCAGAGCTTTGGGAGGCCAAGGCGGGCTCATCACATGAGGTCAGGACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAGTACAAAAATTAGCCAGGCGTGGTGGCACATGCCTGAAATCCCAGCTACTCGGG);
    #my $resulted_ref_seq = Sanger::CGP::VcfCompare::VcfMergeAndPileup::_get_dna_segment($expected_obj->{'samplec'},$g_pu_new->{'chr'},$g_pu_new->{'pos_5p'},$g_pu_new->{'pos_3p'});
    
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
    my($info_tag_val,$vcf_file_obj)=$vcf_obj->getVcfHeaderData;
    my($bam_objects,$bas_files)=$vcf_obj->_get_bam_object;
    #lib size information only applicable for indels....
    my($bam_header_data,$lib_size)=$vcf_obj->_get_bam_header_data($bam_objects,$bas_files);
    #create variant object
    my $variant=Sanger::CGP::Vaf::Process::Variant->new( 
        'location'         => undef,
        'varLine'         => undef,
        'varType'         => $vcf_obj->{'_a'},
        'libSize'         => defined $lib_size?$lib_size:undef,
        'samples'         => $vcf_obj->{'allSamples'},
        'tumourName'    => $vcf_obj->getTumourName,
        'normalName'    => $vcf_obj->getNormalName,
        'vcfStatus'     => $vcf_obj->{'vcf'},
        'noVcf'            => defined $vcf_obj->{'noVcf'}?$vcf_obj->{'noVcf'}:undef,
        'outDir'            => $vcf_obj->getOutputDir,
        'passedOnly'  => $vcf_obj->{'_r'},
    'tabix_hdr'   => defined  $vcf_obj->{'_hdr'}?Bio::DB::HTS::Tabix->new(filename => $vcf_obj->{'_hdr'}):undef,
        );
       $variant->setLocation('1:16902712:T:C');
        $variant->setVarLine('samplea-UM;MN;MQ');
        my($g_pu)=$variant->formatVarinat();
        is_deeply($g_pu,$expected_g_pu,'ReadVcf_processMergedLocations:g_pu');
         $g_pu=$variant->populateHash($g_pu,'samplea',$bam_header_data);
         $g_pu=$variant->getPileup($bam_objects->{'samplea'},$g_pu);
};

subtest 'CleanTestResults' => sub {
      my $vcf_obj = Sanger::CGP::Vaf::Data::ReadVcf->new($options);
  my($cleaned1) = $vcf_obj->cleanTempdir($options->{'tmp'}); 
        is_deeply($cleaned1,$options->{'tmp'},'CleanTestResults:tmpdir');
        my($cleaned2) = $vcf_obj->cleanTempdir($test_output);
        is_deeply($cleaned2,$test_output,'CleanTestResults:testOutdir');
};

system("rm $test_genome");
done_testing();
