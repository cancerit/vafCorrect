#!/use/bin/perl

##########LICENCE############################################################
# Copyright (c) 2016-2017 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of cgpVAF.
#
# cgpVAF is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##############################################################

BEGIN {
  $SIG{__WARN__} = sub {warn $_[0] unless(($_[0] =~ m/^Use of uninitialized value \$buf/) || ($_[0] =~ m/gzip: stdout: Broken pipe/))};
};

use strict;
use Bio::DB::HTS::Tabix;
use File::Path qw(mkpath);
use FindBin qw($Bin);
use English qw( -no_match_vars );
use Pod::Usage qw(pod2usage);
use warnings FATAL => 'all';
use Carp;
use Const::Fast qw(const);
use Getopt::Long;
use Data::Dumper;
use Try::Tiny qw(try catch finally);
use Capture::Tiny qw(:all);

use Log::Log4perl;
use Log::Log4perl::Level;
Log::Log4perl->init("$Bin/../config/log4perl.vaf.conf");
use lib "$Bin/../lib";
use Sanger::CGP::Vaf::Data::ReadVcf;
use Sanger::CGP::Vaf::VafConstants;

my $log = Log::Log4perl->get_logger(__PACKAGE__);

my $tags=$Sanger::CGP::Vaf::VafConstants::SNP_TAGS;

try {
    my ($options) = option_builder();
    $log->debug("================Using Parameters===========================");
    $log->debug(Dumper($options));
    if ($options->{'a'} eq 'indel') {
        $tags=$Sanger::CGP::Vaf::VafConstants::INDEL_TAGS;
    }
    my $vcf_obj = Sanger::CGP::Vaf::Data::ReadVcf->new($options);

    # progress checked before the processing starts , speed ups concatenation step
    my ($chromosomes, $chr_count)=$vcf_obj->getChromosomes($options->{'chr'});
    $log->logcroak('It appears that no data has been detected for processing, please check command line params.') if($chr_count == 0);
    $chromosomes=$vcf_obj->getProgress($chromosomes);

    # this is called only once to add allSample names to vcf object
    $vcf_obj->getAllSampleNames;
    my($info_tag_val,$vcf_file_obj)=$vcf_obj->getVcfHeaderData();
    my($bam_objects,$bas_files)=$vcf_obj->get_bam_object();
    my $max_lib_size = $vcf_obj->get_lib_n_read_read_length($bam_objects,$bas_files);
    # create variant object
    my $variant=Sanger::CGP::Vaf::Process::Variant->new(
        'location'      => undef,
        'tmp'           =>$options->{'tmp'},
        'varLine'       => undef,
        'varType'       => $vcf_obj->{'_a'},
        'samples'       => $vcf_obj->{'allSamples'},
        'tumourName'    => $vcf_obj->getTumourName,
        'normalName'    => $vcf_obj->getNormalName,
        'vcfStatus'     => $vcf_obj->{'vcf'},
        'noVcf'         => defined $vcf_obj->{'noVcf'}?$vcf_obj->{'noVcf'}:undef,
        'lib_size'      => $max_lib_size,
        'outDir'        => $vcf_obj->getOutputDir,
        'passedOnly'    => $vcf_obj->{'_r'},
        'tabix_hdr'     => defined $vcf_obj->{'_hdr'}?Bio::DB::HTS::Tabix->new(filename => $vcf_obj->{'_hdr'}):undef,
        'mq'            => $vcf_obj->{'_mq'},
        'bq'            => $vcf_obj->{'_bq'},
        'exp'           => $vcf_obj->{'_exp'},
        'emb'           => $vcf_obj->{'_emb'},
    );

    foreach my $chr(keys %$chromosomes) {
        my($data_for_all_samples,$unique_locations)=$vcf_obj->getMergedLocations($chr, $vcf_file_obj);
        if(defined $options->{'b'} ){
            my ($bed_locations)=$vcf_obj->getBedHash($chr);
            if( $options->{'bo'} == 1 && (defined $data_for_all_samples) ) {
                ($data_for_all_samples,$unique_locations)=$vcf_obj->populateBedLocationsFromVCF($data_for_all_samples,$unique_locations,$bed_locations);
            }else{
                ($data_for_all_samples,$unique_locations)=$vcf_obj->populateBedLocations($data_for_all_samples,$unique_locations,$bed_locations);
            }
        }
        $vcf_obj->processMergedLocations($data_for_all_samples ,$unique_locations,$variant,$chromosomes->{$chr} ,$bam_objects,$chr,$tags,$info_tag_val);
    }# completed all chromosomes;

    # run following steps only if chromosome option is empty or user has selected option to concatenate files.
    if($options->{'ct'} || @{$options->{'chr'}} == 0 ) {
        if($options->{'m'}) {
            $log->info('Completed analysis for all PASS locations, writing non PASS variants');
            my($aug_vcf_fhs,$aug_vcf_names)=$vcf_obj->WriteAugmentedHeader();
            foreach my $sample (keys %$aug_vcf_names){
            if(-e $aug_vcf_names->{$sample}.'.gz') {
            $log->logcroak(sprintf 'Output file : %s.gz exists , skipping concatenation step', $aug_vcf_names->{$sample});
            }
            $vcf_obj->catFiles($options->{'tmp'},'vcf',$sample,$aug_vcf_names->{$sample});
            $log->info('Compressing and Validating augmented VCF file for sample: '.$sample);
            my($aug_gz,$aug_tabix)=$vcf_obj->gzipAndIndexVcf($aug_vcf_names->{$sample});
            }
        }

        my($outfile_name_no_ext)=$vcf_obj->writeFinalFileHeaders($info_tag_val,$tags);

        if($outfile_name_no_ext){
            $vcf_obj->catFiles($options->{'tmp'},'vcf',undef,$outfile_name_no_ext);
            $vcf_obj->catFiles($options->{'tmp'},'tsv',undef,$outfile_name_no_ext);
            $log->info("Compressing and Validating VCF file");
            my($outfile_gz,$outfile_tabix)=$vcf_obj->gzipAndIndexVcf("$outfile_name_no_ext.vcf");
            if ((-e $outfile_gz) && (-e $outfile_tabix) && !$options->{'dbg'} ) {
            my($cleaned)=$vcf_obj->check_and_cleanup_dir($options->{'tmp'});
            }
        }
    }
}
catch {
  croak "\n\n".$_."\n\n" if($_);
};

# get options from user

sub option_builder {
        my ($factory) = @_;
        my %options = ( 'emb' => 6000,
                        'exp' => 92,
                        'oe' => '.vaf.vcf',
                        's' => undef,
                        't' => 'VD,VW,VT,VC',
                        'be' => '.bam',
                        'finc' => $Sanger::CGP::Vaf::VafConstants::DEFAULT_READLEN_INCLUDE,
                        'fexc' => $Sanger::CGP::Vaf::VafConstants::DEFAULT_READS_EXCLUDE_PILEUP,
                        'bo' => 0,
         );
        &GetOptions (
                'h|help'    => \$options{'h'},
                't|infoTags=s' => \$options{'t'},
                'd|inputDir=s' => \$options{'d'},
                'b|bedIntervals=s' => \$options{'b'},
                'e|vcfExtension=s' => \$options{'e'},
                'be|bamExtension=s' => \$options{'be'},
                'g|genome=s' => \$options{'g'},
                'a|variant_type=s' => \$options{'a'},
                'r|restrict_flag=i' => \$options{'r'},
                'chr|chromosome=s{,}' => \@{$options{'chr'}},
                'ct|concat=i'  => \$options{'ct'},
                'o|outDir=s'  => \$options{'o'},
                'm|augment' => \$options{'m'},
                'mq|map_quality=i' => \$options{'mq'},
                'bq|base_quality=i' => \$options{'bq'},
                # provide at least 1 tumour sample name
                'tn|tumour_name=s{1,}' => \@{$options{'tn'}},
                'nn|normal_name=s' => \$options{'nn'},
                'tb|tumour_bam=s{1,}' => \@{$options{'tb'}},
                'nb|normal_bam=s' => \$options{'nb'},
                'bo|bed_only=i' => \$options{'bo'},
                'oe|output_vcfExtension=s' => \$options{'oe'},
                'tmp|tmpdir=s' => \$options{'tmp'},
                'dp|depth=s' => \$options{'dp'},
                'hdr|high_depth_bed=s' => \$options{'hdr'},
                'pid|id_int_project=s' => \$options{'pid'},
                'exp|exonerate_pct=i' => \$options{'exp'},
                'vcf|vcf_files=s{,}' => \@{$options{'vcf'}},
                'finc|filter_inc=i' => \$options{'finc'},
                'fexc|filter_exc=i' => \$options{'fexc'},
                'emb|exonerate_mb=i' => \$options{'emb'},
                'dbg|debug' => sub { $log->level($TRACE) },
                'v|version'  => \$options{'v'}
    );

  pod2usage(-message => Sanger::CGP::Vaf::license, -verbose => 1) if(defined $options{'h'});

    if(defined $options{'v'}){
        my $version = Sanger::CGP::Vaf->VERSION;
        print "$version\n";
        exit;
    }

    pod2usage(q{'-g' genome must be defined}) unless (defined $options{'g'});
    pod2usage(q{'-d' input directory path must be defined}) unless (defined $options{'d'});
    pod2usage(q{'-a' variant type must be defined}) unless (defined $options{'a'});
    pod2usage(q{'-tn' toumour sample name/s must be provided}) unless (defined $options{'tn'});
    pod2usage(q{'-nn' normal sample name/s must be provided}) unless (defined $options{'nn'});
    pod2usage(q{'-e' Input vcf file extension must be provided}) unless (defined $options{'bo'} || defined $options{'e'} || defined $options{'vcf'});
    pod2usage(q{'-b' bed file must be specified }) unless (defined $options{'b'} || defined $options{'e'} || defined $options{'vcf'});
    pod2usage(q{'-o' Output folder must be provided}) unless (defined $options{'o'});

    $options{'d'}=~s/\/$//g;
    mkpath($options{'o'});
    if(!defined $options{'tmp'}) {
        $options{'o'}=~s/\/\//\//g;
        my $tn_name=@{$options{'tn'}}[0];
        mkpath($options{'o'}.'/tmpvaf_'.$tn_name);
        $options{'tmp'}=$options{'o'}.'/tmpvaf_'.$tn_name;
    }
    if(defined $options{'a'} and ( (lc($options{'a'}) eq 'indel') || (lc($options{'a'}) eq 'snp') ) ) {
        $log->info(sprintf 'Analysing: %s', $options{'a'});
    }
    else{
        $log->logcroak('Not a valid variant type [should be either [snp or indel]');
        exit(0);
    }

    # use PASS flag if bed only is set
    if($options{'bo'}==1) {
        $options{'r'} = 0;
    }
    if(!defined $options{'r'}) {
        $options{'r'} = 1;
    }

    if($options{'a'} eq 'indel' && !defined $options{'dp'}) {
        $options{'dp'}='NR,PR';
    }

    if($options{'ct'}) {
        $log->debug('Concatenation option selected, chromosome option will be set to all chromosomes');
        $options{'chr'}=[];
    }

    if($options{'m'} && lc($options{'a'}) eq 'snp' ) {
        $log->logcroak('VCF augment option is only supported for indels');
    }
  if(!defined $options{'hdr'} && lc($options{'a'}) eq 'indel') {
     $log->warn('-hdr high depth reagions file not provided for indel analysis, high depth regions will take longer to run');
    }
 \%options;
}

__END__

=head1 NAME

cgpVaf.pl merge the variants in vcf files for a given Tumour - normal pairs in a project  and write pileup and exonerate output for merged locations

=head1 SYNOPSIS

cgpVaf.pl [-h] -d -a -g -tn -nn -e  -o [ -tb -nb -b -t -c -r -m -ao -mq -pid -bo -vcf -v]

  Required Options (inputDir and variant_type must be defined):

   --variant_type         (-a)     variant type (snp or indel) [default snp]
   --inputDir             (-d)     input directory path containing bam and vcf files
   --genome               (-g)     genome fasta file name (default genome.fa)
   --tumour_name          (-tn)    Tumour sample name(s), space separated [ co-located bam/cram, index and bas files expected, unless -tb specified]
   --normal_name          (-nn)    Normal sample name [ co-located bam/cram, index and bas files expected, unless -nb specified ]
   --outDir               (-o)     Output folder

  Optional
   --vcfExtension         (-e)     vcf file extension string after the sample name - INCLUDE preceding dot (default: .caveman_c.annot.vcf.gz)
                                   - optional if -bo 1
                                   - ignored  when -vcf defined
   --tumour_bam           (-tb)    tumour bam/cram file(s) space separated list [ optional if -tn is specified]
                                   - if not defined will be deduced from tumour_name
                                   - should be specified in same order as tumour sample names
   --normal_bam           (-nb)    normal bam/cram file [optional if -nn is specified]
                                   - if not defined will be deduced from --normal_name
   --infoTags             (-t)     comma separated list of tags to be included in the tsv output, vcf file by default includes all data
                                   (default: VD,VW,VT,VC for Vagrent annotations)
   --bedIntervals         (-b)     tab separated file containing list of intervals in the form of <chr><pos> <ref><alt> (e.g 1  14000  A  C)
   --restrict_flag        (-r)     restrict analysis (possible values 1 : PASS or 0 : ALL) [default 1 ]
   --chromosome           (-chr)   restrict analysis chromosomes [space separated]
   --concat               (-ct)    concat per chromosome results to a single vcf file
   --augment              (-m)     Augment original vcf file (valid for indels)
                                   - this will add additional fields[ MTR, WTR, AMB] to FORMAT column of NORMAL and TUMOUR samples ] (default FALSE)
   --output_vcfExtension  (-oe)    Extension for augmented VCF, see --augment [vaf.vcf]
   --map_quality          (-mq)    read mapping quality threshold
   --base_quality         (-bq)    base quality threshold for snp
   --exonerate_pct        (-exp)   report alignment over a percentage of the maximum score attainable by each query (exonerate/indel parameter) [default 92]
   --exonerate_mb         (-emb)   Max memory to alow Exonerate to use (indel) [default: 6000]
   --depth                (-dp)    comma separated list of field(s) as specified in FORMAT field representing total depth at given location
   --high_depth_bed       (-hdr)   High Depth Region(HDR) bed file (tabix indexed) to mask high depth regions in the genome
   --bed_only             (-bo)    Only analyse bed intervals in the file (default 0: analyse vcf and bed interval)
   --vcf                  (-vcf)   user defined input vcf file path(s) [ optional if -tn is specified ]
                                   - if not defined will be deduced from --tumour_name and --vcfExtension
                                   - should be specified in same order as tumour sample names
   --filter_inc           (-finc)  Sam flag values to include when checking reads for read length
   --filter_exc           (-fexc)  Sam flag values to exclude when checking reads for read length

    --id_int_project       (-pid)   Internal project id [WTSI only]
   --help                 (-h)     Display this help message
   --version              (-v)     provide version information for vaf

   Examples:
      #Merge vcf files to create single vcf containing union of all the variant sites and provides pileup output for each location
      perl cgpVaf.pl -d tmpvcfdir -o testout -a snp -g genome.fa -e .caveman_c.annot.vcf.gz -nn PD21369b -tn PD26296a PD26296c2
      #Merge vcf files to create single vcf containing union of all the variant sites and provides allele count for underlying indel location
      perl cgpVaf.pl -d tmpvcfdir -o testout -a indel -g genome.fa -e .caveman_c.annot.vcf.gz -nn PD21369b -tn PD26296a PD26296c2
      # Run per chromosome analysis
      perl cgpVaf.pl -d tmpvcfdir -o testout -a indel -g genome.fa -e .caveman_c.annot.vcf.gz -nn sampleb -tn samplea samplec -chr 1
      perl cgpVaf.pl -d tmpvcfdir -o testout -a indel -g genome.fa -e .caveman_c.annot.vcf.gz -nn sampleb -tn samplea samplec -chr 2
      # concatenate per chromosome output to single vcf
      perl cgpVaf.pl -d tmpvcfdir -o testout -a indel -g genome.fa -e .caveman_c.annot.vcf.gz -nn sampleb -tn samplea samplec -ct 1

=cut
