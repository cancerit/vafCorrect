#Package- merge VCF files and run pileup for SNP and exonerate for indels

package Sanger::CGP::VcfCompare::CreateCommands; 
use Sanger::CGP::VcfCompare;
our $VERSION = Sanger::CGP::VcfCompare->VERSION;

BEGIN {
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) .
   ($_[0] =~ m/^Odd number of elements in hash assignment/) || ($_[0] =~m/^Use of uninitialized value \$gtype/) || ($_[0] =~ m/^Use of uninitialized value \$buf/)|| ($_[0] =~ m/symlink exists/) || ($_[0] =~ m/gzip: stdout: Broken pipe/) )};

};

$main::SQL_LIB_LOC = '.'; # this suppresses warnings about uninitialised values
use strict;
use Tabix;
use Config::IniFiles;

use Bio::DB::Sam::Constants;
use Data::Dumper;
use English;
use File::Path qw(mkpath);
use FindBin qw($Bin);
use List::Util qw(first reduce max min);
use Math::Round qw(round);
use Carp;
use Const::Fast qw(const);
use Capture::Tiny qw(:all);
use Try::Tiny qw(try catch finally);
use warnings FATAL => 'all';

use Log::Log4perl;
Log::Log4perl->init("$Bin/../config/log4perl.vaf.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

const my $SEP => "/";
const my $test_mode => 0;
const my $version => Sanger::CGP::VcfCompare->VERSION;
const my $ref_dir => '/nfs/cancer_ref01';
const my $ref_dir_x10 => '/lustre/scratch112/sanger/cgppipe/canpipe/live/ref';



use base qw(Sanger::CGP::Vaf::Data::AbstractVcf);


sub _localInit {
	my $self=shift;
} 


##### modules from create_config ######

=head2 load_sql_config
load sql statement in connection object
Inputs
=over 2
=item conn -sql connection object
=back
=cut

sub loadSql {
  my ($self,$conn) = @_;
	$conn->addQuery('nst::NL::getProjectBamAndVcf', q{
	select cs.id_sample cs_id_sample, cs.id_ind, ip.id_int_project, sip.sample_synonym, LOWER(s.species) SPECIES, ipat.build, ipat.design, ipat.sw, sipa.attr_value treat_as_tumour
	, max(decode(ar.result_type, 247, ar.result,decode(ar.result_type,7,ar.result))) BAM
	, max(decode(ar.result_type, 16, ar.result)) BAI
	, max(decode(ar.result_type, 250, ar.result)) BAS
	, max(decode(ar.result_type,140, ar.result)) CAVE
	, max(decode(ar.result_type,141, ar.result)) CAVE_IDX
	, max(decode(ar.result_type,183, ar.result)) CAVE_C
	, max(decode(ar.result_type,184, ar.result)) CAVE_C_IDX
	, max(decode(ar.result_type,132, ar.result)) PINDEL
	, max(decode(ar.result_type,133, ar.result)) PINDEL_IDX
	, max(decode(ar.result_type,27, ar.result)) PINDEL_BAM
	, max(decode(ar.result_type,28, ar.result)) PINDEL_BAI
	from (
		select id_int_project
		, max(decode(attr_type, 1, attr_value)) build
		, max(decode(attr_type, 3, attr_value)) design
		, nvl(max(decode(attr_type, 10, attr_value)),0) sw
		from internal_project_attributes
		group by id_int_project
	) ipat
	, internal_project ip
	, sample_int_project sip
	, sip_attributes sipa
	, sample s
	, analysis_results ar
	,cosi_summary cs
	where ip.team_name > 0
	and ip.id_int_project = ipat.id_int_project
	and ip.id_int_project = sip.id_int_project
	and sip.id_sample = s.id_sample
	and sip.id_sample_cosmic = sipa.id_sample_cosmic
	and sip.id_sample_cosmic = cs.id_sample
	and sipa.attr_type = 12
	and ip.id_int_project = ar.id_int_project
	and sip.id_sample_cosmic = ar.id_field
	and ar.is_current = 1
	and ar.result_type in (7,247,16,250,140,141,183,184,132,133,27,28)
	and ip.id_int_project = ?
	group by ip.id_int_project, sip.sample_synonym, s.species, ipat.build, ipat.design, ipat.sw, sipa.attr_value,cs.id_sample,cs.id_ind
	order by 1,2
		});
}

=head2 buildInputData
parse sql results 
Inputs
=over 2
=item options -user provide input parameters
=item conn -sql connection object
=back
=cut


sub buildInputData {
  my ($self,$conn) = @_;
  my $project_id=$self->{'_p'};
		my $all_data = $conn->executeArrHashRef('nst::NL::getProjectBamAndVcf',$project_id);
		my @retained_data;
		my $total_records = 0;
		for my $curr(@{$all_data}) {
			$total_records++;
			next if( (scalar @{$self->{'_u'}}) > 0 && (! first { $curr->{'SAMPLE_SYNONYM'} eq $_ } @{$self->{'_u'}}) && ($curr->{'TREAT_AS_TUMOUR'} eq 'Y') );
			next if(defined $self->{'_p'} && ! first { $curr->{'ID_INT_PROJECT'} == $_ } $self->{'_p'});
			if(@retained_data > 0) {
				croak "There are multiple species in your requested analysis.\n" if($curr->{'SPECIES'} ne $retained_data[-1]->{'SPECIES'});
				croak "There are multiple builds in your requested analysis.\n" if($curr->{'BUILD'} ne $retained_data[-1]->{'BUILD'});
			}
			push @retained_data, $curr;
		}
		print "Total samples  for this project: ".scalar @retained_data."\n";
	  return \@retained_data;
}
=head2 generateRawOutput
generate raw output using SQL results hash
Inputs
=over 2
=item options -user provide input parameters
=item to_process -results obtained by SQL query
=back
=cut

sub generateRawOutput {
  my ($self, $to_process) = @_;
	my (%sample_group, %normal_samples, $species, $build,$symlinked_files);
  for my $sample_row(@{$to_process}) {
  	if ($sample_row->{'TREAT_AS_TUMOUR'} eq 'N') {
    	$normal_samples{$sample_row->{'ID_INT_PROJECT'}}{$sample_row->{'ID_IND'}}=$sample_row->{'SAMPLE_SYNONYM'}; 
    	$symlinked_files=$self->_processFiles($sample_row,$symlinked_files);	
    	next;
    }
    if(!exists $sample_group{$sample_row->{'ID_INT_PROJECT'}}{$sample_row->{'ID_IND'}} ) {
    	$sample_group{$sample_row->{'ID_INT_PROJECT'}}{$sample_row->{'ID_IND'}}=$sample_row->{'SAMPLE_SYNONYM'};
		  $species = $sample_row->{'SPECIES'};
			$build = $sample_row->{'BUILD'};
		}
		else {
		  $sample_group{$sample_row->{'ID_INT_PROJECT'}}{$sample_row->{'ID_IND'}}.="\t$sample_row->{'SAMPLE_SYNONYM'}";
		}
   $symlinked_files=$self->_processFiles($sample_row,$symlinked_files);
	}
  
  return (\%sample_group, \%normal_samples, $species, $build, $symlinked_files );
}
=head2 _process_files
process the input data files
Inputs
=over 2
=item options -user provide input parameters
=item input_files -input files to process
=item symlinked_files -hash storing symlinked file path information

=back
=cut

sub _processFiles {
  my ($self, $input_files,$symlinked_files) = @_;
  my $root_path = $self->{'_o'};
  my $out_file = $root_path;
  my ($project, $sample);
  my ($file_types)=$self->_getFileTypes();
  my $file_ref = ref $input_files;
  $sample = $input_files->{'SAMPLE_SYNONYM'};
  if((defined $file_ref) && $file_ref eq 'HASH') {
    foreach my $file (keys %$file_types) {
      if(defined $input_files->{$file}) {
      	#if(defined ($input_files->{'CAVE_C'}) and (($file eq 'CAVE') or ($file eq 'CAVE_IDX'))) {next;} 
        my $sym_link=$root_path."/$sample.".$file_types->{$file};
        $self->_createSymlink($input_files->{$file}, $sym_link);
        push(@$symlinked_files,$sym_link);
      }	
    }  
  } 
 $symlinked_files;
}
=head2 _get_file_types
hash storing various files extension one gets from the database
Inputs
=over 2
=item
=back
=cut

sub _getFileTypes {
	my($self)=shift;
	my %file_types=('BAM' => 'bam',
	                'BAI' => 'bam.bai',
	                'BAS' => 'bam.bas',
	                'CAVE' => 'cave.annot.vcf.gz',
	                'CAVE_IDX' => 'cave.annot.vcf.gz.tbi',
	                'CAVE_C' => 'caveman_c.annot.vcf.gz',
	                'CAVE_C_IDX' => 'caveman_c.annot.vcf.gz.tbi',
	              	'PINDEL' => 'pindel.annot.vcf.gz',
	              	'PINDEL_IDX' => 'pindel.annot.vcf.gz.tbi',
	              	'PINDEL_BAM' => 'mt.pindel.bam',
                  'PINDEL_BAI' => 'mt.pindel.bam.bai',
	              	);         	
  return \%file_types;			
}

=head2 write_config
Write configuration file 
Inputs
=over 2
=item options -user provide input parameters
=item sample_group -hash storing sample data as a group per individual
=item normal_samples -hash storing normal samples
=item project -project name
=item species -species these samples belongs to
=build -reference sequence build
=ref_dir -path to reference build
=symlinked_files -list fo files for which symlink is created
=back
=cut


sub writeConfig {
  my ($self,$sample_group, $normal_samples, $species, $build, $symlinked_files)=@_;
	my ($group_name,$unmatched_normal_count, $data_in_config)='';
	my ($file_types)=$self->_getFileTypes();
	#my $data_in_config=0;
	my $single_sample_count=0;
	my $root_path = $self->{'_o'};
	my $project = $self->{'_p'};
  my $cfg=Config::IniFiles->new();
  my $config_path="$root_path/$project\_VcfMergeConfig.ini";
  $cfg->SetFileName($config_path);
  $cfg->AddSection($project); 
  foreach my $id_ind(keys %{$sample_group->{$project}}) {
		my @samples_to_analyse=split("\t",$sample_group->{$project}{$id_ind});
		if(@samples_to_analyse >0) {
			
			if(!exists $normal_samples->{$project}{$id_ind} && (defined $self->{'_n'} && $self->{'_n'} eq 'Y')) {
				next;
			} 
			if(!exists $normal_samples->{$project}{$id_ind}) {
				$unmatched_normal_count++;
				$group_name=join("_",@samples_to_analyse)."_$unmatched_normal_count";
			}
			else {
			  $group_name="$normal_samples->{$project}{$id_ind}";
			}
			
			my $subset_flag=0;
	  	if((scalar @{$self->{'_u'}}) > 0) {
	  		$subset_flag=$self->_checkSubset(\@samples_to_analyse);
	  	}
	  	if ($subset_flag) {
			 	next;
	  	}
	  	$data_in_config++;
			$cfg->newval($project,$group_name,@samples_to_analyse);
			
		}
	  else {
		  $single_sample_count++;
		}
	}
	
	if (!defined $data_in_config) {
  	$log->warn("No data to merge for this project exit.....");
  	exit(0);
  }
	$log->info("Number of tumour-normal pairs in config file: $data_in_config");
	if (defined $unmatched_normal_count) {
		$log->info("Nunmber of unmatched tumour-normal pairs: $unmatched_normal_count");
	}
	#add bed file to config file
	if(defined $self->{'_b'}) {
		$cfg->AddSection('UserData'); 
		$cfg->newval('UserData','bedfile',$self->{'_b'});
	}
	#create symlink for genome build
	$cfg=$self->_symlinkGenomeBuild($species,$build,$cfg,$root_path,$ref_dir,$ref_dir_x10);
  $cfg->RewriteConfig();
  
  $config_path,$root_path;
}


sub _checkSubset {
	my ($self,$tumour_sample_names)=@_;
	foreach my $user_sample(@{$self->{'_u'}}) {
		foreach my $cofig_sample(@$tumour_sample_names) {
			if($user_sample eq $cofig_sample) {
				return 0;
			}
		}
	}
	$log->debug("No user sample in config:".@$tumour_sample_names);
	return 1;
}


=head2 _symlink_genome_build
crete symlink for genome build
Inputs
=over 2
=item species -species these samples belongs to
=item build -reference sequence build
=item cfg -config file object
=item root_path -path to store symlinked files
=ref_dir -path to reference build
=back
=cut

sub _symlinkGenomeBuild {
	my ($self,$species,$build,$cfg,$root_path,$ref_dir,$ref_dir_x10)=@_;
	my $genome_fasta="$root_path/genome.fa";
	my $genome_index="$root_path/genome.fa.fai";
	my $ref_fasta="$ref_dir/$species/$build/genome.fa";
	my $ref_index="$ref_dir/$species/$build/genome.fa.fai";
	if(! -e $ref_fasta) {
		$ref_fasta="$ref_dir_x10/$species/$build/genome.fa";
		$ref_index="$ref_dir_x10/$species/$build/genome.fa.fai";
	}
	
	$self->_createSymlink($ref_fasta, $genome_fasta);
	$self->_createSymlink($ref_index, $genome_index);
	$cfg->AddSection('genome_build'); 
  $cfg->newval('genome_build', 'genome', $species);
  $cfg->newval('genome_build', 'build', $build);
  return $cfg;
}

=head2 run_vcfmerge
run vcf merge by calling the perl script
Inputs
=over 2
=item config_path -config file path
=item input_dir -input folder storing all the data
=item resp -user response to run specific algorithm
=back
=cut


sub createVafCommands {
	my($config_path,$input_dir,$resp)=@_;
	#--varinat_type  (-v)  specify variant data type (default SNV [ specify snp or indel])
  #--vcf_from      (-f)  specify vcf source [ only used with varinat_type snp (default cave_c [specify-: cave_c or cave_java]  [ WTSI only ]
	#--outDir        (-o)  Output folder [default to inputdir/outptut]
  #--restrict_flag (-r)  restrict analysis on PASS filter  [default 1 ](possible values 1 PASS or 0 ALL)                      
  
	my ($filter,$algorithm);
	if($resp=~/^(1|2|3|4|5|6|7)$/) {
	  print "\033[2J"; # clears the screen
	  print "\033[0;0H"; # jump to 0,0
		print "Restrict analysis on (possible values 1 : PASS varinats only  or 0 : ALL varinats) [default 1 ] =>:";
		$filter = <STDIN>;
		chomp $filter;
		if($filter!~/^(0|1)$/) { $filter=1;}
	}

	print "Output will be stored at:\n $input_dir/pindel/output/ \n $input_dir/caveman_c/output/ \n and/or $input_dir/caveman_java/output/ \n";
	
	if($resp == 1 || $resp == 4 || $resp == 5 || $resp == 7 ) {
		print "Merging pindel vcf files:\n";
		my $cmd="$Bin/mergeAndPileup.pl -i $config_path -d $input_dir --outDir $input_dir/pindel --variant_type indel  --vcf_from pindel -r $filter --augment 1";
		$log->debug("COMMAND: $cmd");
		system("$cmd");				
	}
	if($resp == 2 || $resp == 5 || $resp == 6 || $resp == 7 ) {
		print "Merging caveman_java vcf files:\n";
		my $cmd="$Bin/mergeAndPileup.pl -i $config_path -d $input_dir --outDir $input_dir/caveman_java --variant_type snp  --vcf_from cave_java -r $filter ";
		$log->debug("COMMAND: $cmd");
		system("$cmd");
	}
	if($resp == 3 || $resp == 4 || $resp == 6 || $resp == 7 ) {
	print "Merging caveman_c vcf files:\n";
		my $cmd="$Bin/mergeAndPileup.pl -i $config_path -d $input_dir --outDir $input_dir/caveman_c --variant_type snp  --vcf_from caveman_c -r $filter ";
		$log->debug("COMMAND: $cmd");
		system("$cmd");
		
	}

}


##### End of create config modules ########

######Common modules ######################

=head2 _create_symlink
Create symlink
Inputs
=over 2
=item file
destination file for which symlink will be created
=item symlink_path
path where symlink will be created
=back
=cut

sub _createSymlink {
	my ($self,$file, $symlink_path)=@_;
	if( -l $symlink_path) { $log->warn("symlink exists, skipped symlink creation"); return;} 
	symlink $file, $symlink_path;
}


###### End of common modules ##############

sub print_hash {
	my $hash=shift;
	#print "---------------Printing data---------------\n";
	foreach my $key (sort keys %$hash) {
		if(defined($hash->{$key}))
		{	
			print "$key:==>$hash->{$key}\n";	
		}
	}
}

1;
