#!/software/perl-5.16.3/bin/perl 

##########LICENCE############################################################
# Copyright (c) 2016 Genome Research Ltd.
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
  use Cwd qw(abs_path);
  use File::Basename;
  $ENV{POSIXLY_CORRECT}=1;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) || ($_[0] =~ m/^Use of uninitialized value \$buf/)|| ($_[0] =~ m/symlink exists/) || ($_[0] =~ m/gzip: stdout: Broken pipe/) )};
};

use strict;
use warnings FATAL => qw(all);
$main::SQL_LIB_LOC = '.'; # this suppresses warnings about uninitialised values
use FindBin qw($Bin);

#use lib("/software/CGP/projects/vcfCommons/perl/lib");
use Sanger::CGP::Config::Config qw(vcfCommons);
use Sanger::CGP::Vaf::Utility::CreateCommands;

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

try {
  my ($options) = option_builder();
  my $base = Sanger::CGP::WholeGenome::Base->new;
  $base->db_type($options->{'db'});
  print "Using database $options->{'db'}\n";
  my $conn = $base->connection;
  
  my $vaf=Sanger::CGP::Vaf::Utility::CreateCommands->new($options);
  $vaf->loadSql($conn);
  my $data_to_process = $vaf->buildInputData($conn);
  $base->clear_connection;
  my ($sample_group,$normal_samples, $species, $build,$symlinked_files)= $vaf->generateRawOutput($data_to_process);
  my($cfg)=$vaf->writeConfig($sample_group,$normal_samples, $species, $build,$symlinked_files,$conn);
  undef $conn;
  if(defined $cfg) {
		print "Please specify vcf file types to merge: \n".
		 "[1] Pindel \n".
		 "[2] caveman_java \n".
		 "[3] caveman_c \n".
		 "[4] (1,3) \n".
		 "[5] (1,2) \n".
		 "[6] (2,3) \n".
		 "[7] (all) \n".
		 "[0] exit  \n------:"; 
		my $resp = <STDIN>;
		chomp $resp;
		if($resp=~/^(1|2|3|4|5|6|7)$/) {
			$vaf->createVafCommands($resp,$cfg);
			print "\n >>>>>> To analyse the data please use commands from commands.txt file created in current directory >>>>>>>>>\n";	
		}
		else{
			warn "Please provide valid option... exiting...\n";
			exit(0);
		}
	}
}
catch {
 croak "\n\n".$_."\n\n" if($_);
};

sub option_builder {
	my ($factory) = @_;
	my %options;
	&GetOptions (
	'h|help'    => \$options{'h'},
	'pid|id_project_int=i' => \$options{'pid'},
	'b|bedIntervals=s' => \$options{'b'},
	'o|outdir=s'  => \$options{'o'},
	'db|databse_type=s'  => \$options{'db'},
	'u|sample_names=s{,}' => \@{$options{'u'}},
	#additional parameters to run VAF
	't|infoTags=s' => \$options{'t'},
	'c|hdr_cutoff=i' => \$options{'c'},
	'g|genome=s' => \$options{'g'},
	'r|restrict_flag=i' => \$options{'r'},
	'm|augment=i' => \$options{'m'},
	'ao|augment_only=i' => \$options{'ao'},
	# provide at least 1 tumour samples name
	'bo|bed_only=i' => \$options{'bo'},
	'oe|output_vcfExtension=s' => \$options{'oe'},
	'dp|depth=s' => \$options{'dp'},
	'vn|vcf_normal=i' => \$options{'vn'},
	'v|version'  => \$options{'v'},
	);

  pod2usage(-message => Sanger::CGP::VcfCompare::license, -verbose => 1) if(defined $options{'h'});
	
	if(defined $options{'v'}){
		my $version = Sanger::CGP::VcfCompare->VERSION;
		print "$version\n";
		exit;
	}
	pod2usage(q{'-pid' Project identifier must be provided.}) unless(defined $options{'pid'});
	pod2usage(q{'-o' Output folder must be specified.}) unless(defined $options{'o'}) ;
	$options{'db'}='live' unless(defined $options{'db'});
	unless(-e $options{'o'}){
	  mkpath($options{'o'});
	}
	# check if directory path has / at the end other wise add it
	$options{'o'} .= '/' unless($options{'o'} =~ m|/$|);
	return \%options;
}

__END__

=head1 NAME

createVafCmd.pl - Create commands to run for each sample group containing tumour normal sample pairs in a given project   

=head1 SYNOPSIS

createVafCmd.pl [-h -v] -pid -o  [-g -b -u -db -t -c -r -m -ao -p -bo -vn]

  Required Options (project and output directory must be defined):

    --project_id    (-pid)project id [e.g 888]
    --outdir        (-o)  outdir [ Path to output folder ]
  Optional
    --genome        (-g) 	genome fasta file name (default genome.fa)
    --bedIntervals  (-b) 	tab separated file containing list of intervals in the form of <chr><pos> <ref><alt> (e.g 1  14000  A  C)
    --sample_names  (-u) 	Restrict to the list of samples to analyse [ e.g., PD12345a PD12345c]
    --database_type (-db)	database type [live] e.g., test or live
    --infoTags      (-t) 	comma separated list of tags to be included in the vcf INFO field 
                         	(default: VD,VW,VT,VC for Vagrent annotations)
    --hdr_cutoff    (-c) 	High Depth Region(HDR) cutoff  value[ avoids extreme depth regions (default: 005 i.e top 0.05% )]
                         	(possible values 001,005,01,05 and 1)
    --restrict_flag (-r) 	restrict analysis on (possible values 1 : PASS or 0 : ALL) [default 1 ]   
    --augment       (-m) 	Augment pindel file [ this will add additional fields[ MTR, WTR, AMB] to FORMAT column of NORMAL and TUMOUR samples ] (default 0: don not augment)
    --augment_only  (-ao)	Only augment pindel VCF file (-m must be specified) [ do not  merge input files and add non passed varinats to output file ] (default 0: augment and merge )
    --bed_only      (-bo)	Only analyse bed intervals in the file (default 0: analyse vcf and bed interval)
    --vcf_normal    (-vn)	use normal sample defined in vcf header field [ default 1 ]
    --depth         (-dp)	comma separated list of field(s) as specified in FORMAT field representing total depth at given location
    --help          (-h) 	Display this help message
    --version       (-v) 	displays version number of this software

  Examples:
    - Create commands file for following project 
      perl createVafCmd.pl -pid 888 -o testdir
=cut

