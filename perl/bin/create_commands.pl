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

use File::Path qw(mkpath);
use FindBin qw($Bin);
use English qw( -no_match_vars );
use warnings FATAL => 'all';
use Carp;
use Try::Tiny qw(try catch finally);
use Capture::Tiny qw(:all);
use Config::IniFiles;

my $ref_dir = '/nfs/cancer_ref01';
my $ref_dir_x10 = '/lustre/scratch112/sanger/cgppipe/canpipe/live/ref';

my($indir,$ini_file,$ext)=@ARGV;

if(!$indir || !$ini_file ||!$ext) {
	croak("unable to find files");
}

my $cfg;
try {
	if(-e "$indir/$ini_file") {
       		$cfg = Config::IniFiles->new( -file => "$indir/$ini_file" );
	}
  else {
	croak("No ini file defined or tumour normal samples specified on command line");
	}
	foreach my $project_id ($cfg->Sections()) { # project_name
		if($project_id eq 'genome_build' || $project_id eq 'UserData') {next;}
		foreach my $sample_group($cfg->Parameters($project_id)) {  # sample group
	  	my @tumour_sample_names = $cfg->val($project_id, $sample_group);
			print "perl /nfs/users/nfs_s/sb43/scripts/cgpVAF/perl/bin/vaf.pl --outDir $indir/output/$sample_group\_out -g $indir/genome.fa --variant_type indel -nn $sample_group -tn @tumour_sample_names -e $ext \n";  
		}  
	}
}catch {
  croak "\n\n".$_."\n\n" if($_);
};

