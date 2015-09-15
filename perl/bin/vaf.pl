#!/software/perl-5.16.3/bin/perl
BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) || ($_[0] =~ m/^Use of uninitialized value \$buf/) || ($_[0] =~ m/gzip: stdout: Broken pipe/))};
};

use strict;
#$main::SQL_LIB_LOC = '.'; # this suppresses warnings about uninitialised values
#use Sanger::CGP::Config::Config qw(vcfCommons);
#use Sanger::CGP::VcfCompare::VcfMergeAndPileup;

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

use Sanger::CGP::Vaf::Data::ReadVcf;
use Sanger::CGP::Vaf::VafConstants;

my $log = Log::Log4perl->get_logger(__PACKAGE__);

const my $SEP => "/";

try {
	my ($options) = option_builder();
	my $vcf = Sanger::CGP::Vaf::Data::ReadVcf->new($options);
	my($data_for_all_samples,$unique_locations,$info_tag_val,$updated_info_tags)=$vcf->getUniqueLocations();
	$vcf->processLocations($data_for_all_samples,$unique_locations,$info_tag_val,$updated_info_tags);
	
		
}
	
catch {
  croak "\n\n".$_."\n\n" if($_);
};


# get options from user

sub option_builder {
        my ($factory) = @_;
        my %options;
        &GetOptions (
                'h|help'    => \$options{'h'},
                't|infoTags=s' => \$options{'t'},
                'b|bedIntervals=s' => \$options{'b'},
                'e|vcfExtension=s' => \$options{'e'},
                'c|hdr_cutoff=s' => \$options{'c'},
                'g|genome=s' => \$options{'g'},
                'a|variant_type=s' => \$options{'a'},
                'f|vcf_from=s' => \$options{'f'},
                'r|restrict_flag=s' => \$options{'r'},
                'o|outDir=s'  => \$options{'o'},
                'u|sample_names=s' => \$options{'u'},
                'm|augment=s' => \$options{'m'},
                'ao|augment_only=s' => \$options{'ao'},
                'tn|tumour_name=s' => \$options{'tn'},
                'nn|normal_name=s' => \$options{'nn'},
                'tb|tumour_bam=s' => \$options{'tb'},
                'nb|normal_bam=s' => \$options{'nb'},
                'vcf|input_vcf=s' => \$options{'vcf'},
                'bo|bed_only=s' => \$options{'bo'},
                'oe|output_vcfExtension=s' => \$options{'oe'},
                'ie|input_vcfExtension=s' => \$options{'oe'},
                'p|depth=s' => \$options{'p'},
                'v|version'  => \$options{'v'}
	);
	
  pod2usage(-message => Sanger::CGP::Vaf::license, -verbose => 1) if(defined $options{'h'});
        
	if(defined $options{'v'}){
		print $Sanger::CGP::Vaf::VafConstants::VERSION."\n";
		exit;
	}
	pod2usage(q{'-g' genome must be defined}) unless (defined $options{'g'});
	pod2usage(q{'-a' variant type must be defined}) unless (defined $options{'a'});
	pod2usage(q{'-tb' toumour sample bam files must be provided}) unless (defined $options{'tb'});
	pod2usage(q{'-nb' normal sample bam file must be provided}) unless (defined $options{'nb'});
	pod2usage(q{'-vcf' Input vcf file name must be provided}) unless (defined $options{'vcf'} || defined $options{'bo'});
	pod2usage(q{'-b' bed file must be specified }) unless (!defined $options{'vcf'} || !defined $options{'bo'});
  pod2usage(q{'-o' Output folder must be provided}) unless (defined $options{'o'});
	if(!defined $options{'bo'}) { $options{'bo'}=0;}
	mkpath($options{'o'});
	if(!defined $options{'e'}) { # variant extension
		if(defined $options{'a'} and lc($options{'a'}) eq 'indel'){
				$options{'e'}=".pindel.annot.vcf.gz";	
		}
		elsif(defined $options{'a'} and lc($options{'a'}) eq  'snp') {
			if(defined $options{'f'} and lc($options{'f'}) eq "cave_java") {
				$options{'e'}=".cave.annot.vcf.gz";
				}
			else{
				$options{'e'}=".caveman_c.annot.vcf.gz";
			}
		}	
		else{
			$log->logcroak("Not a valid variant type [should be either [snp or indel]");	
			exit(0); 
		}	
 	} 		 
	
	if(!defined $options{'t'}) { 
		$options{'t'}="VD,VW,VT,VC";
	}
	#use tabix file 
	if(!defined $options{'c'}) {
		$options{'c'}='005';
	}
	# use PASS flag
	if(!defined $options{'r'}) {
		$options{'r'}= 1;
	}
	if($options{'a'} eq 'indel' && !defined $options{'p'}) {
		$options{'p'}='NR,PR';
	}
	
	if(!defined $options{'s'}) {
		#analyse single sample no merge step 
		$options{'s'}=undef;
	}
	if(!defined $options{'ao'}) {
		# augment vcf no merging step
		$options{'ao'}=undef;
	}
	if(!defined $options{'oe'}) {
		# augment vcf extesnion
		$options{'oe'}='.vaf.vcf';
	}
		
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
   --outDir        (-o)  Output folder
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

