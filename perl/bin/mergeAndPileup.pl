#!/software/perl-5.16.3/bin/perl
BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) || ($_[0] =~ m/^Use of uninitialized value \$buf/) || ($_[0] =~ m/gzip: stdout: Broken pipe/))};
};

use strict;
$main::SQL_LIB_LOC = '.'; # this suppresses warnings about uninitialised values
use Sanger::CGP::Config::Config qw(vcfCommons);
use Sanger::CGP::VcfCompare::VcfMergeAndPileup;

use File::Path qw(mkpath);
use FindBin qw($Bin);
use English qw( -no_match_vars );
use Pod::Usage qw(pod2usage);
use warnings FATAL => 'all';
use Carp;
use Const::Fast qw(const);
use Getopt::Long;
use Try::Tiny qw(try catch finally);
use Capture::Tiny qw(:all);

my $log = Log::Log4perl->get_logger(__PACKAGE__);

const my $SEP => "/";

my $progress_fhw;
my $progress_fhr;
my $cfg;
try {
	my ($options) = option_builder();
	$options=Sanger::CGP::VcfCompare::VcfMergeAndPileup::get_options($options);
	# create progress file is it doesn't exists else open existing file in append mode
	print "\n >>>>>> To view analysis progress please check  .log file created in current directory >>>>>>>>>\n";	
	if (-e "$options->{'o'}/progress.out") {
		open($progress_fhw, '>>', "$options->{'o'}/progress.out");
	}
	else {
		open($progress_fhw, '>', "$options->{'o'}/progress.out");
	}
	open($progress_fhr, '<', "$options->{'o'}/progress.out");
	my @progress_data=<$progress_fhr>;
	close($progress_fhr);
	# create tmp ini
	if (!defined $options->{'i'} && $options->{'tn'} && $options->{'nn'}) {
		$log->debug("Config file not defined will be created from input paramaters");
		$options= Sanger::CGP::VcfCompare::VcfMergeAndPileup::create_tmp_ini($options);
	}  
	# create tabix object
	my $tabix_hdr = new Tabix(-data => "$Bin/hdr/seq.cov".$options->{'c'}.'.ONHG19_sorted.bed.gz');
	if($options->{'i'}) {
       		$cfg = Config::IniFiles->new( -file => $options->{'i'} );
	}
        else {
	$log->logcroak("No ini file defined or tumour normal samples specified on command line");
	}
  $log->debug("Using config file:".$options->{'i'});
	$log->debug("Running analysis for:".$options->{'a'});
	$log->debug("Using HDR [High Depth regions] cutoff:".$options->{'c'});
	$log->debug("Analysis progress log will be written in:".$options->{'o'}."/progress.out");
        
	if(!defined $options->{'b'}) {	
		$options->{'b'}=$cfg->val('UserData','bedfile')
	}
	foreach my $project_id ($cfg->Sections()) { # project_name
		if($project_id eq 'genome_build' || $project_id eq 'UserData') {next;}
		foreach my $sample_group($cfg->Parameters($project_id)) {  # sample group
			my $outfile_name=$options->{'o'}.$SEP.${project_id}."_".${sample_group}."_consolidated_".$options->{'a'};
			# check file name in progress file
			my($progress_flag)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::check_progress(\@progress_data,$outfile_name,$options);
			next if $progress_flag;
	  	my @tumour_sample_names = $cfg->val($project_id, $sample_group);
	  	my $subset_flag=0;
	  	if ( !$options->{'m'} && !$options->{'s'} && @tumour_sample_names < 2 ) { 
	  		$log->warn("Single sample in group nothing to merge : @tumour_sample_names"); 
				$log->warn("To run analysis on single sample please set option --single_sample=1"); 
	  		next;
	  	}
	  	if(defined $options->{'u'}) {
	  		$subset_flag=Sanger::CGP::VcfCompare::VcfMergeAndPileup::check_user_subset(\@tumour_sample_names,$options);
	  	}
	  	next if $subset_flag;
			$log->debug("Analysing data for =>project:".$project_id."Samples".@tumour_sample_names);
			#get unique locations hash for all samples...
			my ($union_of_locations_original, $data_for_all_samples,$info_tag_val,$vcf_flag,$normal_sample,$vcf_file_status,$augment_vcf)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::get_unique_locations(\@tumour_sample_names,$options,$sample_group);
			if ($vcf_flag == 0 && $options->{'a'} eq 'snp' && !defined $options->{'f'} ) {
				$log->warn("No caveman_c vcf files for these samples, try specifying option[ -f cave_java ]");
		  	next;
		  }
		  elsif($vcf_flag == 0 && (!$options->{'b'}|| $options->{'ao'} || $options->{'m'})) {
		  	$log->warn("No data to analyse for these samples, no VCF or bed file specified");
		  	next;
		  }  
		 Sanger::CGP::VcfCompare::VcfMergeAndPileup::run_and_consolidate_mpileup(\@tumour_sample_names,$union_of_locations_original,$data_for_all_samples,$info_tag_val,$normal_sample,$outfile_name,$options,$tabix_hdr,$vcf_file_status,$progress_fhw,$augment_vcf,$sample_group);
		}
	}	
}
	
catch {
  croak "\n\n".$_."\n\n" if($_);
}

finally {

close($progress_fhw);
};


# get options from user

sub option_builder {
        my ($factory) = @_;
        my %options;
        &GetOptions (
                'h|help'    => \$options{'h'},
                'i|configFile=s' => \$options{'i'},
                'd|inpuDir=s' => \$options{'d'},
                't|infoTags=s' => \$options{'t'},
                'b|bedIntervals=s' => \$options{'b'},
                'e|vcfExtension=s' => \$options{'e'},
                'c|hdr_cutoff=s' => \$options{'c'},
                'g|genome=s' => \$options{'g'},
                'a|variant_type=s' => \$options{'a'},
                'f|vcf_from=s' => \$options{'f'},
                'r|restrict_flag=s' => \$options{'r'},
                's|single_sample=s' => \$options{'s'},
                'o|outDir=s'  => \$options{'o'},
                'u|sample_names=s' => \$options{'u'},
                'm|augment=s' => \$options{'m'},
                'ao|augment_only=s' => \$options{'ao'},
                'tn|tumour_name=s' => \$options{'tn'},
                'nn|normal_name=s' => \$options{'nn'},
                'tb|tumour_bam=s' => \$options{'tb'},
                'nb|normal_bam=s' => \$options{'nb'},
                'vcf|input_vcf=s' => \$options{'vcf'},
                'pid|id_int_project=s' => \$options{'pid'},
                'bo|bed_only=s' => \$options{'bo'},
                'oe|output_vcfExtension=s' => \$options{'oe'},
                'p|depth=s' => \$options{'p'},
                'v|version'  => \$options{'v'}

	);
	
  pod2usage(-message => Sanger::CGP::VcfCompare::license, -verbose => 1) if(defined $options{'h'});
        
	if(defined $options{'v'}){
		my $version = Sanger::CGP::VcfCompare->VERSION;
		print "$version\n";
		exit;
	}
	#pod2usage(q{'-i' config file must be provided}) unless (defined $options{'i'});
	pod2usage(q{'-d' input directory must be provided}) unless (defined $options{'d'});
	pod2usage(q{'-a' variant type must be defined}) unless (defined $options{'a'});
	
\%options;
}

__END__

=head1 NAME

mergeAndPileup.pl merge the variants in vcf files for a given Tumour - normal pairs in a project  and write pileup and exonerate output for merged locations

=head1 SYNOPSIS

mergeAndPileup.pl [-h] -i -d -a [-t -b -e -c -r -g -f -o -v]

  Required Options (config file, inputDir and variant_type must be defined):

   --configFile    (-i)  path to config file listing tumour/normal samples names without 
                         file extension, this file can be generated  manually or using script (consolidate_results.pl)
   --inputDir      (-d)  path to directory containing links to vcf files to be merged, BAM files and reference genome.
   --variant_type  (-a)  specify variant data type (snp or indel) [default snp]
  Optional
   --infoTags      (-t)  comma separated list of tags to be included in the vcf INFO field 
                         (default: VD,VW,VT,VC for Vagrent annotations)
   --bedIntervals  (-b)  tab separated file containing list of intervals in the form of <chr><pos> <ref><alt> (e.g 1  14000  A  C)
                         bed file can be specified in the config file after the last sample pair,
                         if specified on command line then same bed file is used for all the tumour/normal pairs,
                         bed file name in config file overrides command line argument
   --vcfExtension  (-e)  specify the vcf file extension string after the sample name (default: caveman_c.annot.vcf.gz) 
   --hdr_cutoff    (-c)  specify High Depth Region(HDR) cutoff  value[ avoids extreme depth regions (default: 005 i.e top 0.05% )]
                         (possible values 001,005,01,05 and 1)
   --restrict_flag (-r)  restrict analysis on (possible values 1 : PASS or 0 : ALL) [default 1 ]   
   --single_sample (-s)  run analysis on single sample (default 0: skip single sample vcf )                    
   --genome        (-g)  specify genome fasta file name (default genome.fa)
   --vcf_from      (-f)  specify vcf source [ only used with varinat_type snp (default cave_c [specify-: cave_c or cave_java]  [ WTSI only ]
   --outDir        (-o)  Output folder [default to inputdir/outptut]
   --sample_names  (-u)  Comma separated list of samples within same project [ e.g., PD12345a,PD12345c]
   --augment       (-m)  Augment pindel file [ this will add additional fields[ MTR, WTR, AMB] to FORMAT column of NORMAL and TUMOUR samples ] (default 0: don not augment)
   --augment_only  (-ao) Augment pindel file [ do not  merge input files and add non passed varinats in output file ] (default 0: augment and merge )
   --tumour_name  (-tn)  Toumour sample name [ list of comma separated files ]
   --normal_name  (-nn)  Normal sample name [ only if ini file is not provide ]
   --tumour_name  (-tb)  Toumour sample bam files [ list of comma separated files ] [ only if ini file is not provide ]
   --normal_name  (-nb)  Normal sample bam file [ only if ini file is not provide ]
   --input_vcf    (-vcf) Input vcf file name [Comma separated list if more than one tumour sample provided]
   --depth         (-p)  comma separated list of field(s) as specified in FORMAT field representing total depth at given location
   --id_int_project  (-pid) Internal project id 
   --bed_only      (-bo) Only analyse bed intervals in the file (default 0: analyse vcf and bed interval)
   --vcfExtension  (-oe) Specify uncompressed output vcf file extension string after the sample name (default: .vaf.vcf) only for augmented vcf files. 
   --help          (-h)  This help
   --version       (-v)  provides version information of this code

   Examples:
      Merge vcf files to create single vcf containing union of variant sites and provide pileup output for each location
      perl mergeAndPileup.pl -i exampleConfig.ini -d tmpvcfdir -o testout -a snp
      Merge vcf files to create single vcf containing union of variant sites and provide exonerate output for each location
      perl mergeAndPileup.pl -i exampleConfig.ini -d tmpvcfdir -o testout -a indel
=cut

