#!/software/perl-5.16.3/bin/perl
BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) || ($_[0] =~ m/^Use of uninitialized value \$buf/)|| ($_[0] =~ m/symlink exists/) || ($_[0] =~ m/gzip: stdout: Broken pipe/) )};
};

use strict;
$main::SQL_LIB_LOC = '.'; # this suppresses warnings about uninitialised values
use FindBin qw($Bin);

use Sanger::CGP::Config::Config qw(vcfCommons);
use Sanger::CGP::VcfCompare::VcfMergeAndPileup;

use Sanger::CGP::Database::Conn;
use Sanger::CGP::WholeGenome::Base;

use Config::IniFiles;
use English qw( -no_match_vars );
use Pod::Usage qw(pod2usage);
use Carp;
use File::Path qw(mkpath);
use Getopt::Long;
use Try::Tiny qw(try catch finally);
use Capture::Tiny qw(:all);
use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);
use Data::Dumper;

const my $ref_dir => '/nfs/cancer_ref01';
try {
  my ($options) = option_builder();
  my $project=$options->{'p'};
  my $base = Sanger::CGP::WholeGenome::Base->new;
  $base->db_type('live');
  my $conn = $base->connection;
  Sanger::CGP::VcfCompare::VcfMergeAndPileup::load_sql_config($conn);
  my $to_process = Sanger::CGP::VcfCompare::VcfMergeAndPileup::build_input_data($options, $conn);
  undef $conn;
  $base->clear_connection;
  my ($sample_group,$normal_samples, $species, $build,$symlinked_files)= Sanger::CGP::VcfCompare::VcfMergeAndPileup::generate_raw_output($options, $to_process);
  my($config_path,$output_dir)=Sanger::CGP::VcfCompare::VcfMergeAndPileup::write_config($options, $sample_group, $normal_samples, $project ,$species, $build,$ref_dir, $symlinked_files);
  if($config_path) {
		print "Would you like to merge VCF files? if YES please specify vcf file type: \n".
		 "[1] Pindel \n".
		 "[2] caveman_java \n".
		 "[3] caveman_c \n".
		 "[4] (1,3) \n".
		 "[5] (1,2) \n".
		 "[6] (2,3) \n".
		 "[7] (all) \n".
		 "[0] exit  \n ------:";
		 
				my $resp = <STDIN>;
				chomp $resp;
				if($resp=~/^(1|2|3|4|5|6|7)$/) {
					print "\n >>>>>> To view analysis progress please check .log file created in current directory >>>>>>>>>\n";	
					Sanger::CGP::VcfCompare::VcfMergeAndPileup::run_vcfmerge($config_path,$output_dir,$resp)
				}
				else{
					print "exiting...\n";
					exit(0);
				}
	}
}
catch {
 croak "\n\n".$_."\n\n" if($_);
}

finally {
#print "Script completed ----\n";
};



sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
					'h|help'    => \$opts{'h'},
					'p|project=s' => \$opts{'p'},
					'b|userBed=s' => \$opts{'b'},
					'o|outdir=s'  => \$opts{'o'},
					'n|normal=s'  => \$opts{'n'},
					'u|sample_names=s'  => \$opts{'u'},
					'v|version'  => \$opts{'v'},
	);

  pod2usage(-message => Sanger::CGP::VcfCompare::license, -verbose => 1) if(defined $opts{'h'});
	
	if(defined $opts{'v'}){
		my $version = Sanger::CGP::VcfCompare->VERSION;
		print "$version\n";
		exit;
	}
	
	pod2usage(q{'-p' project number must be defined.}) unless(defined $opts{'p'});
	pod2usage(q{'-o' output folder must be specified.}) unless(defined $opts{'o'}) ;
	
	unless(-e $opts{'o'}){
	  mkpath($opts{'o'});
	}
	# check if directory path has / at the end other wise add it
	$opts{'o'} .= '/' unless($opts{'o'} =~ m|/$|);

	return \%opts;
}

__END__

=head1 NAME

createConfig.pl - Create config file containing tumour normal sample pairs to be merged   

=head1 SYNOPSIS

createConfig.pl [-h] -p -o  [ -b -n - u -v ]

  Required Options (project must be defined):

    --project       (-p) project number [e.g 888]
    --outdir        (-o) outdir [ Path to output folder ]
  Optional

    --help          (-h)  This message and format of input file
     One or more of the following:
     
    --bedIntervals  (-b) tab separated file containing list of intervals in the form of <chr><pos> <ref><alt> (e.g 1  14000  A  C)
    --normal        (-n)  BOOLEAN   Only interrogate matched tumour normal pairs [Y/N default = N (all)].
                    Requires '-p'
    --sample_names  (-u)  comma separated list of samples within same project [ e.g., PD12345a,PD12345c]
    --version       (-v) displays version number of this software

  Examples:
    - Create config file containing tumour normal sample pairs to be merged 
      perl createConfig.pl -p 888 -o testdir
=cut

