#Package- merge VCF files and run pileup for SNP and exonerate for indels

package Sanger::CGP::VcfCompare::VcfMergeAndPileup; 
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
use Vcf;
use Bio::DB::Sam;
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
 
#open(my $new_interval,'>','test_modified_pos_new.tsv');

const my $LIB_MEAN_INS_SIZE => 'mean_insert_size';
const my $LIB_SD => 'insert_size_sd';
const my $EXECUTE_EXTERNAL => 1;
const my $EXONERATE_SCORE_MULTIPLIER => 5;
const my $EXONERATE_SCORE_FACTOR => 70;
const my $READ_LENGTH_CUTOFF => 2;
const my $SEP => "/";
const my $test_mode => 0;
const my $INSERT_SIZE_FACTOR => 1;
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
const my $version => Sanger::CGP::VcfCompare->VERSION;
const my $BASIC_COLUMN_TITLES => ['Normal', 'VariantID','Chrom','Pos','Ref','Alt','Qual','Filter','Gene','Transcript','RNA', 'CDS','Protein', 'Type', 'Effect'];
const my $BASIC_COLUMN_DESCS => ['Normal sample name', 'ID of the Variant','Chromosome','Position, the position of a sub or the position immediatly before an indel',
                        'Reference sequence of the Variant','Alternative sequence of the Variant','VCF Quality field','VCF Filter field',
                        'Gene affected, according to the VD tag','Transcript affected, according to the VD tag','Variant description at RNA level',
                        'Variant description at CDS level','Variant description at protein level','Variant type','Summary variants of overall effect'];

=head2 get_options
get the user input and adds file extensions based on input parameters
Inputs
=over 2
=item options - user provided options to get file extension, input dir path , tag values and bed file
=back
=cut

sub get_options {
	my ($options) = shift;
	if(!defined $options->{'o'}) { $options->{'o'}=$options->{'d'};}
	$options->{'o'}=$options->{'o'}.$SEP."output";
	mkpath($options->{'o'});
	if(!defined $options->{'bo'}) { $options->{'bo'}=0;}
	if(!defined $options->{'g'}) { $options->{'g'}=$options->{'d'}.$SEP."genome.fa"}
	if(!defined $options->{'e'}) { # variant extension
		if(defined $options->{'a'} and lc($options->{'a'}) eq 'indel'){
				$options->{'e'}=".pindel.annot.vcf.gz";	
		}
		elsif(defined $options->{'a'} and lc($options->{'a'}) eq  'snp') {
			if(defined $options->{'f'} and lc($options->{'f'}) eq "cave_java") {
				$options->{'e'}=".cave.annot.vcf.gz";
				}
			else{
				$options->{'e'}=".caveman_c.annot.vcf.gz";
			}
		}	
		else{
			$log->logcroak("Not a valid variant type [should be either [snp or indel]");	
			exit(0); 
		}	
 	} 		 
 	# use annotation tags
	if(!defined $options->{'t'}) { 
		$options->{'t'}="VD,VW,VT,VC";
	}
	#use tabix file 
	if(!defined $options->{'c'}) {
		$options->{'c'}='005';
	}
	# use PASS flag
	if(!defined $options->{'r'}) {
		$options->{'r'}='1';
	}
	if($options->{'a'} eq 'indel' && !defined $options->{'p'}) {
		$options->{'p'}='NR,PR';
	}
	
	if(!defined $options->{'s'}) {
		#analyse single sample no merge step 
		$options->{'s'}=0;
	}
	if(!defined $options->{'ao'}) {
		# augment vcf no merging step
		$options->{'ao'}=0;
	}
	if(!defined $options->{'oe'}) {
		# augment vcf extesnion
		$options->{'oe'}='.vaf.vcf';
	}
$options;
}


=head2 create_tmp_ini
create tmp_ini files and symlinks for for associated input files provided as input paramaters
Inputs
=over 2
=item options
=back
=cut


sub create_tmp_ini {
	my ($options)=@_;
	chomp $options->{'tn'};
  $log->debug("Tumour Name".$options->{'tn'});
	my @tn_samples=split(',',$options->{'tn'});
	my $cfg=Config::IniFiles->new();
	my $config_path="vafConfig.ini";
  $cfg->SetFileName($config_path);
  $cfg->AddSection($options->{'pid'}); 
  $cfg->newval($options->{'pid'},$options->{'nn'},@tn_samples);
  $log->debug("Paramaters in tmp ini files:".$options->{'pid'}.':'.$options->{'nn'}.':'.@tn_samples);
  $cfg->WriteConfig($config_path);	
  $options->{'i'}='vafConfig.ini';
  
  if($options->{'tb'} && $options->{'nb'}) {
  	_create_bam_symlinks($options);
  }
  if($options->{'vcf'}) {
  	_create_vcf_symlinks($options);
  }
  $options;
}

=head2 _create_vcf_symlinks
create symlink for bam files to be accessible in input folder
Inputs
=over 2
=item options
=back
=cut

sub _create_bam_symlinks {
	my ($options)=@_;
	
	my @tn_samples=split(',',$options->{'tn'});
	my @tn_bam=split(',',$options->{'tb'});
	if (scalar(@tn_samples) != scalar(@tn_bam)){
		warn "Unequal number of tumour sample names and bam files\n";
		exit;
	}
	# create tumor sample symlinks
	foreach my $i (0 .. $#tn_samples) {
	  if (-e $tn_bam[$i] ) {
			_create_symlink($tn_bam[$i], $options->{'d'}.'/'.$tn_samples[$i].'.bam');
			_create_symlink($tn_bam[$i].'.bai', $options->{'d'}.'/'.$tn_samples[$i].'.bam.bai');
		}
		if (-e $tn_bam[$i].'.bas') {
			_create_symlink($tn_bam[$i].'.bas', $options->{'d'}.'/'.$tn_samples[$i].'.bam.bas');
		}
	}
	# create normal sample symlinks
	if (-e $options->{'nb'} ) {
		_create_symlink($options->{'nb'}, $options->{'d'}.'/'.$options->{'nn'}.'.bam');
		_create_symlink($options->{'nb'}.'.bai', $options->{'d'}.'/'.$options->{'nn'}.'.bam.bai');
	}
	if (-e $options->{'nb'}.'.bas' ) {
		_create_symlink($options->{'nb'}.'.bas', $options->{'d'}.'/'.$options->{'nn'}.'.bam.bas');
	}
	
 return 1;
}

=head2 _create_vcf_symlinks
create symlink for vcf files to be accessible in input folder
Inputs
=over 2
=item options
=back
=cut

sub _create_vcf_symlinks {
	my ($options)=@_;
	my @tn_samples=split(',',$options->{'tn'});
	my @vcf=split(',',$options->{'vcf'});
	if (scalar(@tn_samples) != scalar(@vcf)){
		warn "Unequal number of tumour sample names and vcf files\n";
		exit;
	}
	# create vcf symlink
	foreach my $i (0 .. $#tn_samples) {
	  if (-e $vcf[$i] && $vcf[$i]=~m/vcf.gz$/) {
			_create_symlink($vcf[$i], $options->{'d'}.'/'.$tn_samples[$i].'.vcf.gz');
			_create_symlink($vcf[$i].'.tbi', $options->{'d'}.'/'.$tn_samples[$i].'vcf.gz.tbi');
		}
		else {
			warn "Provide bgzipped and tabix indexed vcf file\n";
			exit;
		}
	}
 return 1;
}

=head2 check_progress
check if analysis had been run for this sample group
Inputs
=over 2
=item $progress_data - array of file names for which analysis has been run
=item $outfile_name -- new out file name to be checked in progress file
=item options - user provided options
=back
=cut
sub check_progress {
	my ($progress_data,$outfile_name,$options)=@_;
	my $progress_flag=0;
	foreach my $file (@$progress_data) {
				chomp $file;
				if ($file eq "$outfile_name.vcf") {
					warn "WARNING!!! Outfile: $outfile_name.vcf exists skipping analysis for this sample\n";
					warn "To rerun analysis on all samples please remove $options->{'o'}/progress.out file\n";
					$progress_flag=1;
					return $progress_flag;
				}
	}
$progress_flag;
}

=head2 check_user_subset
restrict analysis for use defined sample group
Inputs
=over 2
=item tumour_sample_names - array of sample names
=item options - user provided options to get file extension, input dir path , tag values and bed file
=back
=cut

sub check_user_subset {
	my ($tumour_sample_names,$options)=@_;
	chomp $options->{'u'};
	my @user_subset=split(',',$options->{'u'});
	foreach my $user_sample(@user_subset) {
		foreach my $cofig_sample(@$tumour_sample_names) {
			if($user_sample eq $cofig_sample) {
				return 0;
			}
		}
	}
	$log->debug("No user sample in config:".@$tumour_sample_names);
	return 1;
}

=head2 get_unique_locations
create union of locations for a set of input files
Inputs
=over 2
=item sample_names - array of sample names
=item options - user provided options to get file extension, input dir path , tag values and bed file
=back
=cut
sub get_unique_locations {
	my ($sample_names,$options,$normal_name)=@_;
	my ($unique_locations, $data_for_all_samples, $info_data, $info_tag_val, $normal_sample, $vcf_file_status,$updated_info_tags,$augment_vcf);
	my $vcf_flag=0;
	my $read_depth=0;
	my $tumour_count=0;
	
	my @depth_tags;
	if(defined $options->{'p'}) {
		 @depth_tags=split(',',$options->{'p'});
	}
	foreach my $sample (@$sample_names) {
		my $vcf_file;
	  if($options->{'vcf'}) {
			$vcf_file=$options->{'vcf'};
			$log->debug("Analysing data for user provided vcf: $vcf_file");
		}		
	  elsif(!$options->{'vcf'} && $options->{'tn'} && $options->{'nn'}) {
			$vcf_file=$options->{'d'}.$SEP.$sample.'_vs_'.$normal_name.$options->{'e'};
			
		}
		elsif($vcf_file && -e $vcf_file && !$options->{'vcf'}) {
		  $log->debug("Analysing data for vcf file in input folder: $vcf_file");
		}
		else {
			$vcf_file=$options->{'d'}.$SEP.$sample.$options->{'e'};
			$log->debug("Analysing data for default vcf: $vcf_file");
		}
		
	  if(-e $vcf_file ) {
	  	$tumour_count++;
	  	my $vcf = Vcf->new(file => $vcf_file);
	  	 if ($options->{'m'} && @$sample_names) {
	  			$augment_vcf->{$sample} = $vcf_file;
	  	}
			$vcf->parse_header();	
			($info_tag_val,$normal_sample,$updated_info_tags)=_get_old_header($vcf,$info_tag_val,$normal_sample,$tumour_count,$sample,$updated_info_tags,$options);
			while ((my $x=$vcf->next_data_array()) && ($options->{'bo'}=~m/0/)) {
			  if (!defined $options->{'t'}) { undef $info_data;}
				else {$info_data=_parse_info($vcf,$$x[7],$updated_info_tags);}
				my $filter_flag=1;
				if($$x[6] eq "PASS") { $filter_flag=2; }
				#location key consists of CHR:POS:REF:ALT
				my $location_key="$$x[0]:$$x[1]:$$x[3]:$$x[4]";
				if(defined $options->{'p'}) {
					$read_depth=_get_read_depth($vcf,$x,$options);
				}
				$data_for_all_samples->{$sample}{$location_key}={'INFO'=>$info_data, 'FILTER'=>$$x[6],'RD'=>$read_depth };	
																																																																																		
				if(!exists $unique_locations->{$location_key}) {
					$unique_locations->{$location_key}="$sample-$$x[6]";
				}
				else{
					$unique_locations->{$location_key}.=':'."$sample-$$x[6]";
				}
			}	
			
			$vcf->close();
			$vcf_flag=1;
			$vcf_file_status->{$sample}=1;
		}
		else {
			$vcf_file_status->{$sample}=0;
		}		
	}
	# add bed locations to unique hash...
	if(defined $options->{'b'} ) {
		($unique_locations, $data_for_all_samples,$info_tag_val)=_populate_bed_locations($unique_locations,$data_for_all_samples,$sample_names,$info_tag_val,$options,$updated_info_tags,$vcf_flag);	
	}	
	if($vcf_flag == 0 && !$options->{'b'} && !$options->{'m'}) { return $unique_locations,$data_for_all_samples,$info_tag_val, $vcf_flag}
	$log->debug("Created unique VCF locations");
	return ($unique_locations,$data_for_all_samples,$info_tag_val,$vcf_flag,$normal_sample,$vcf_file_status,$augment_vcf);
}

=head2 _get_read_depth
get read depth for given variant site for underlying sample
line[8] is format field, get index of tag used for depth e.g., NR(reads on -ve strand) and  PR (reads on +ve strand)
once we had index get value for given tag in the tumour sample
Inputs
=over 2
=item vcf 	-- vcf object
=item line 	-- vcf line for give location
=item tags 	--user specified tags that account for total read depth
=back
=cut

sub	_get_read_depth {
	my ($vcf,$line,$options)=@_;
	my @tags=split(',',$options->{'p'});
	my $depth=0;
	for my $tag (@tags) {
		my $idx = $vcf->get_tag_index($$line[8],$tag,':'); 
		my $pl  = $vcf->get_field($$line[10],$idx) unless $idx==-1;
		$depth+=$pl;
	}
	
$depth;
} 

=head2 _get_old_header
parse VCF header and get data from old header lines
Inputs
=over 2
=item vcf - vcf object
=item info_tag_val - hash to store data as perl VCF header specifications 
declared in get_unique_locations so that samples names for each VCF were stored
=item info_tags -user defined annotation tags to be included in the INFO field [e.g Vagrent annotations]
=item norma_sample -name of the normal sample to be stored 
=item tumour_count -counter to increment tumour sample name in header
=back
=cut

sub _get_old_header {
	my ($vcf,$info_tag_val,$normal_sample,$tumour_count,$sample,$updated_info_tags,$options)=@_;
	#stores hash key data in an array for header tags...
	if ($tumour_count == 1){
		my ($vcf_filter,$vcf_info,$vcf_format)=_custom_header_lines($options);
		#add contig info
		my $contig_info=$vcf->get_header_line(key=>'contig');
		foreach my $contig_line ( @$contig_info) {
			foreach my $chr (sort keys %$contig_line) {
				push(@$info_tag_val,$contig_line->{$chr});
			}
		}
		
		# get user defined tags
		if(defined $options->{'t'}) {
			foreach my $tag (split(",",$options->{'t'})) {
				my $line=$vcf->get_header_line(key=>'INFO', ID=>$tag);
				if(@$line > 0) {
					push (@$updated_info_tags,$tag);
				}
				push(@$info_tag_val,@$line);
			}
		}
		#add custom info
		foreach my $key (sort keys %$vcf_info) {
			push(@$info_tag_val,$vcf_info->{$key});
		}
		
		#add old FILTER
		my $filter_info=$vcf->get_header_line(key=>'FILTER');
		#$filter_info->{(key='FILTER_OLD')}=delete $filter_info->{(key='FILTER')} 
		foreach my $filter_info ( @$filter_info) {
			foreach my $filter_id (sort keys %$filter_info) {
					foreach my $key (keys %{$filter_info->{$filter_id}}) {
						if($filter_info->{$filter_id}{$key} eq 'FILTER') {
							$filter_info->{$filter_id}{$key}='ORIGINAL_FILTER';
						}
					}
				push(@$info_tag_val,$filter_info->{$filter_id});
			}
		}
		#add custom filter
		foreach my $key (sort keys %$vcf_filter) {
			push(@$info_tag_val,$vcf_filter->{$key});
		}
		#add format
		#add custom filter
		foreach my $key (sort keys %$vcf_format) {
			push(@$info_tag_val,$vcf_format->{$key});
		}
		
	}
	#add samples from old header
	my $header_sample=$vcf->get_header_line(key=>'SAMPLE');
	foreach my $header_sample_line (@$header_sample) {
		foreach my $sample_data (sort keys %$header_sample_line) {
			if($sample_data eq "TUMOUR") {
				$header_sample_line->{'TUMOUR'}{'ID'}="TUMOUR$tumour_count";
			}	
		push(@$info_tag_val,$header_sample_line->{$sample_data});
		}
	}
	
	#add sample names..
	my $sample_info=$vcf->get_header_line(key=>'SAMPLE', ID=>'NORMAL');
	foreach my $sample_line ( @$sample_info) {
		foreach my $key (keys %$sample_line) {
			if(defined $sample_line->{$key} and $key eq "SampleName") {
				$normal_sample->{$sample_line->{$key}}.="|$sample";
			}
		}		
	}	

	$info_tag_val,$normal_sample,$updated_info_tags;
}

=head2 _custom_header_lines
stores custom VCF header information as per VCF specifications
Inputs
=over 2
=item vcf_filter - hash created to store FILTER field
=item vcf_info - hash created to store INFO field
=item vcf_foramt - hash created to store FORMAT field
=back
=cut

sub _custom_header_lines {
	my ($options)=shift;
	my ($vcf_filter,$vcf_info,$vcf_format);
  my ($log_key,$process_param)=_get_process_log($options);
  $vcf_format->{'process_log'}={(key=>$log_key,
			InputVCFSource => _trim_file_path($0),
			InputVCFVer => $version,
			InputVCFParam =>$process_param )};

  $vcf_filter->{'filter1'}={(key=>'FILTER',ID=>'1',Description=>"New filter status 1=not called in any")}; # use of 0 is not allowed in FILTER column
	$vcf_filter->{'filter2'}={(key=>'FILTER',ID=>'2',Description=>"New filter status 2=called in any")};
	$vcf_filter->{'filter3'}={(key=>'FILTER',ID=>'3',Description=>"New filter status 3=called + passed")};	
        
	$vcf_info->{'NS'}={(key=>'INFO',ID=>'NS',Number=>'1',Description=>"Number of samples analysed [Excludes designated normal sample]")}; 
	$vcf_info->{'NC'}={(key=>'INFO',ID=>'NC',Number=>'1',Description=>"Number of samples where variant originally Called [Excludes designated normal sample]")};
	$vcf_info->{'NP'}={(key=>'INFO',ID=>'NP',Number=>'1',Description=>"Number of samples where variant Passed the Filter [Excludes designated normal sample]")};
	$vcf_info->{'ND'}={(key=>'INFO',ID=>'ND',Number=>'1',Description=>"Number of samples where sequencing Depth[>0 reads] found [Excludes designated normal sample]")};
	$vcf_info->{'NA'}={(key=>'INFO',ID=>'NA',Number=>'1',Description=>"Number of samples where original algorithm is run [Excludes designated normal sample]")};
	$vcf_info->{'NVD'}={(key=>'INFO',ID=>'NVD',Number=>'1',Description=>"Number of samples where sequencing depth[>0 reads] found for variant reads[Excludes designated normal sample]")};
	
	$vcf_format->{'MTR'}={(key=>'FORMAT',ID=>'MTR', Number=>'1',Type=>'Integer',Description=>"Reads reporting the variant allele")};
	$vcf_format->{'WTR'}={(key=>'FORMAT',ID=>'WTR', Number=>'1',Type=>'Integer',Description=>"Reads reporting the reference allele")};
	$vcf_format->{'DEP'}={(key=>'FORMAT',ID=>'DEP', Number=>'1',Type=>'Integer',Description=>"Total reads covering this position (for subs del positions should be ignored)")};
	$vcf_format->{'MDR'}={(key=>'FORMAT',ID=>'MDR', Number=>'1',Type=>'Integer',Description=>"Variant allele read directions 0=no reads; 1=Forward; 2=Reverse; 3=Forward + Reverse")};
	$vcf_format->{'WDR'}={(key=>'FORMAT',ID=>'WDR', Number=>'1',Type=>'Integer',Description=>"Reference allele read directions 0=no reads; 1=Forward; 2=Reverse; 3=Forward + Reverse")};
	$vcf_format->{'OFS'}={(key=>'FORMAT',ID=>'OFS', Number=>'1',Type=>'String',Description=>"Original filter status as defined in old_FILTER field")};
	if ($options->{'a'} eq 'indel') {
		$vcf_format->{'AMB'}={(key=>'FORMAT',ID=>'AMB', Number=>'1',Type=>'String',Description=>"Reads mapping on both the alleles with same specificity")};
	}
	
	if ($options->{'a'} eq 'snp') {
		$vcf_format->{'FAZ'}={(key=>'FORMAT',ID=>'FAZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a A for this position, forward strand")};
		$vcf_format->{'FCZ'}={(key=>'FORMAT',ID=>'FCZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a C for this position, forward strand")};
		$vcf_format->{'FGZ'}={(key=>'FORMAT',ID=>'FGZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a G for this position, forward strand")};
		$vcf_format->{'FTZ'}={(key=>'FORMAT',ID=>'FTZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a T for this position, forward strand")};
		$vcf_format->{'RAZ'}={(key=>'FORMAT',ID=>'RAZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a A for this position, reverse strand")};
		$vcf_format->{'RCZ'}={(key=>'FORMAT',ID=>'RCZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a C for this position, reverse strand")};
		$vcf_format->{'RGZ'}={(key=>'FORMAT',ID=>'RGZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a G for this position, reverse strand")};
		$vcf_format->{'RTZ'}={(key=>'FORMAT',ID=>'RTZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a T for this position, reverse strand")};	
	}
	
	
	$vcf_filter,$vcf_info,$vcf_format;
}

=head2 _parse_info
parse VCF INFO line for given locations and returns values for respective tags
Inputs 
=over 2
=item vcf - vcf object
=item INFO - VCF INFO filed value
=item INFO - INFO field tags defined in header 
=back
=cut
sub _parse_info {
	my ($vcf, $INFO,$info_tags)=@_;
	my %info_data;
	#my @tags=split(',',$info_tags);
	foreach my $tag (@$info_tags) {
		my $info_val = $vcf->get_info_field($INFO,$tag);
		if($info_val) {
			$info_data{$tag}=$info_val;
		}
		else {
			$info_data{$tag}='-';
		}
	}
	return \%info_data;
}
=head2 _populate_bed_locations
populate additional bed locations for each vcf file
Inputs
=over 2
=item unique_locations -union of locations as created by merging vcf files  
=item data_for_all_samples - hash storing data for all samples as key-val pair for OFS, FILTER and INFO fields.
=item vcf_files -array of sample names in a tumour normal group
=item info_tag_val - INFO filed data for bed locations
=item bed_file -Bed file in tab separated format chr start stop alt and ref allele
=back
=cut
sub _populate_bed_locations {
	my ($unique_locations,$data_for_all_samples,$vcf_files,$info_tag_val,$options,$updated_info_tags,$vcf_flag)=@_;
	my $bed_file=$options->{'b'};
	open(my $bedFH, $options->{'b'})|| $log->logcroak("unable to open file $!");
	my $location_counter=0;
	my %info_tag;
	
	$info_tag{'Interval'}='BedFile';
  
	
	if(!$vcf_flag && !$options->{'ao'}) {
		my ($vcf_filter,$vcf_info,$vcf_format)=_custom_header_lines($options);
		foreach my $key (sort keys %$vcf_info) {
			push(@$info_tag_val,$vcf_info->{$key});
		}
		foreach my $key (sort keys %$vcf_format) {
			push(@$info_tag_val,$vcf_format->{$key});
		}
		foreach my $key (sort keys %$vcf_filter) {
			push(@$info_tag_val,$vcf_filter->{$key});
		}
	}
	
	my %bed_vcf_info=(key=>'FILTER',ID=>'BD', Description=>"Location from bed file");
	push(@$info_tag_val,\%bed_vcf_info);
	while(<$bedFH>) {
		chomp;
		next if $_=~/^#/;
		my $fileds=split("\t",$_);
		if ($fileds < 3) { 
			$log->debug("Not a valid bed location".$_); 
			next;
		}
		my $location_key=join(":",split("\t",$_));
		next if(exists $unique_locations->{$location_key});
		$unique_locations->{$location_key}="$bed_file-BEDFILE";
		#create custom tag hash
		my $temp_tag_val;
		
		foreach my $tag (@$updated_info_tags) { 
			$temp_tag_val->{$tag}='-';
		}
		#add location key to all samples
		foreach my $sample(@$vcf_files)	{
			$data_for_all_samples->{$sample}{$location_key}={'INFO'=>$temp_tag_val,'FILTER'=>'NA','RD'=>1 };
		}
		$location_counter++;
	}
	$log->debug(" Added additional ( $location_counter ) locations from bed file:$bed_file");
	return ($unique_locations,$data_for_all_samples,$info_tag_val);
}

=head2 run_and_consolidate_mpileup
Main wrapper method to call other methods to analyse and write the results
Inputs
=over 2
=item input_bam_files - sample names in a given group
=item unique_locations - union of locations as created by merging vcf files  
=item data_for_all_samples - hash storing data for all samples as key-val pair for OFS, FILTER and INFO fields.
=item vcf_files - array of sample names in a tumour normal group
=item info_tag_val - INFO filed data for bed locations
=item normal_sample -Normal sample name 
=item conn - database connection object
=item outfile_name - merged VCF outfile name
=item options - user provided options to get file extension, input dir path , tag values and bed file
=back
=cut

#create wrapper script for this method

sub run_and_consolidate_mpileup {
	my ($input_bam_files, $unique_locations, $data_for_all_samples,$info_tag_val,$normal_sample,$outfile_name,$options,$tabix_hdr, $vcf_file_status,$progress_fh,$augment_vcf,$sample_group)=@_;
	my $store_results;
	my $sample_counter=0;
	my ($bam_header_data,$lib_size,$aug_vcf_fh,$aug_vcf_name);
	my $no_vcf_counter=0;
	print $progress_fh "WARNING!!! more than one normal sample detected for this group \n".print_hash($normal_sample)."\n" if scalar keys %$normal_sample > 1;
	my ($normal,$normal_flag)=_check_normal_bam($options->{'d'},$normal_sample,$sample_group);
	if($normal_flag){
		unshift(@$input_bam_files,$normal);
		$vcf_file_status->{$normal}=0;
		foreach my $key(keys %$vcf_file_status) {
			if($vcf_file_status->{$key} eq '0' &&  $key ne $normal) {
				$no_vcf_counter++;
			}
		}
	}
	else {
		print $progress_fh "WARNING!! No normal sample BAM found in the input directory skipping analysis for this group @$input_bam_files\n";
		return;
	}
	my @tags=qw(FAZ FCZ FGZ FTZ RAZ RCZ RGZ RTZ MTR WTR DEP MDR WDR OFS);
	my($bam_objects,$bas_files)=_get_bam_object($input_bam_files,$options);
	
	if($options->{'a'} eq 'indel') {
	 	@tags=qw(MTR WTR DEP AMB MDR WDR OFS);
		($bam_header_data,$lib_size)=_get_bam_header_data($bam_objects,$bas_files);
		if($options->{'m'}) {
		 	($aug_vcf_fh,$aug_vcf_name)=_write_augmented_header($augment_vcf,$options,$input_bam_files,$normal);
		}
	}
	my($vcf,$WFH_VCF,$WFH_TSV)=_write_vcf_header($input_bam_files,$info_tag_val,$normal,$vcf_file_status,$options,\@tags,$normal_sample,$outfile_name);
	my $count=0;
	my $total_locations=keys %$unique_locations;
	my $no_pass=0;
	my $loc_counter=0;
	foreach my $location (sort keys %$unique_locations) {	
	  $loc_counter++;
	  #print "------$location---\n";
	  #next if $location !~/49445525/;
		my (%pileup_results,$ref_seq_file,$alt_seq_file,$ref_n_alt_seq_file,$ref_file,$alt_file);
		my($g_pu)=_get_region($location,$options->{'a'});
		my $add_no_pass=1;
		$no_pass++;
		$sample_counter=0;
			if ($options->{'r'} && $unique_locations->{$location}!~/PASS/ && $unique_locations->{$location}!~/BEDFILE/) {		 
						 if($options->{'ao'}) {
					      #write_no_pass_lines_to_vcf
					     foreach my $sample (@$input_bam_files) {
					      $sample_counter++;
					      if($data_for_all_samples->{$sample}{$location}) {
									if($sample_counter > 1) {
					     			$store_results = _store_results($store_results,$g_pu,$add_no_pass,$sample,$unique_locations->{$location},$location);	
					     		}
					     	}
					     }
    					 if($no_pass % 10000 == 0 ) {
								$log->debug("Added Non passed varinat positions:$no_pass"); 
							} 	  
					}
					next;
			}	
		$count++;
		$add_no_pass=0;
				
		my ($original_vcf_info, $NFS, $original_flag,$max_depth)=_get_original_results($input_bam_files, $data_for_all_samples,$location,$vcf_file_status,$no_vcf_counter,$options, $unique_locations->{$location});
		
		if($options->{'a'} eq 'indel') {
			$g_pu=_get_range($bam_header_data,$normal,$g_pu,$tabix_hdr,$max_depth);
			my $ref_seq = _get_dna_segment($bam_objects->{$normal},$g_pu->{'chr'},$g_pu->{'pos_5p'},$g_pu->{'pos_3p'});
			my $reconstructed_alt_seq = _get_alt_seq($bam_objects->{$normal},$g_pu);
			$ref_n_alt_seq_file="$options->{'o'}/temp.ref"; 
			open (my $ref_n_alt_FH,'>'.$ref_n_alt_seq_file);
			print $ref_n_alt_FH ">alt\n$reconstructed_alt_seq\n>ref\n$ref_seq\n";
			close($ref_n_alt_FH);
			
			#print $new_interval "Before\t".$g_pu->{'alt_seq'}."\t".$g_pu->{'region'}."\t".$normal."\t".$g_pu->{'ins_flag'}."\t".$g_pu->{'ref_pos_5p'}."\t".$g_pu->{'ref_pos_3p'}."\t".$g_pu->{'alt_pos_3p'}."\t";
	 
			($g_pu)=_get_ref_5p_pos($ref_seq,$reconstructed_alt_seq,$g_pu);
	
			#print $new_interval "After\t".$g_pu->{'ins_flag'}."\t".$g_pu->{'ref_pos_5p'}."\t".$g_pu->{'ref_pos_3p'}."\t".$g_pu->{'alt_pos_3p'}."\n";
			
		}
		
		$sample_counter=0;
		my $mutant_depth=0;
		my $depth=0;
		foreach my $sample (@$input_bam_files) {
			$sample_counter++;
			$g_pu=_populate_hash($g_pu,$sample,$bam_header_data); # to reset the counter for counts and lib size;
			if($options->{'a'} eq 'indel') {	
				my $temp_read_file="$options->{'o'}/temp.reads";
				open (my $Reads_FH, '>',$temp_read_file) || $log->logcroak("unable to open file $!");
				$g_pu=_fetch_features($bam_objects->{$sample},$g_pu,$Reads_FH,$options);
				$g_pu=_do_exonerate($ref_n_alt_seq_file,$temp_read_file,$g_pu,$alt_file,$ref_file);
				if($sample_counter eq '1' && $options->{'m'}) {
						$g_pu->{'normal_MTR'}=$g_pu->{'alt_p'} + $g_pu->{'alt_n'};
						$g_pu->{'normal_WTR'}=$g_pu->{'ref_p'} + $g_pu->{'ref_n'};
						$g_pu->{'normal_AMB'}=$g_pu->{'amb'};
					}
				if($options->{'m'} && $data_for_all_samples->{$sample}{$location}) {
					if($sample_counter > 1) {
								$store_results=_store_results($store_results,$g_pu,$add_no_pass,$g_pu->{'sample'},$unique_locations->{$location},$location);							
					}	
				}	
				#next if $options->{'ao'};
			}
			# don't follow further steps if augment only option chosen 						
			# for snp data...
			elsif ($options->{'a'} eq 'snp'){
				$g_pu=_get_pileup($bam_objects->{$sample},$g_pu);
				#print Dumper $g_pu;
			}
			if($options->{'ao'} == 0) {
				my($pileup_line)=_format_pileup_line($original_flag,$g_pu,$options);
				#counter where depth is found
				if($g_pu->{'sample'} ne "$normal" && $pileup_line->{'MTR'} > 0 ) {
						$mutant_depth++;
				}
				if($g_pu->{'sample'} ne "$normal" && $pileup_line->{'DEP'} > 0 ) {
						$depth++;
				}
				$pileup_results{$g_pu->{'sample'}}=$pileup_line;
			}
		}
		#get specific annotations from original VCF INFO field....
		if($options->{'ao'} == 0) {
			$original_vcf_info->{'ND'} =	$depth;
			$original_vcf_info->{'NVD'} =	$mutant_depth;	
			_write_output($input_bam_files,$original_vcf_info,$NFS,\%pileup_results,$vcf,\@tags,$WFH_VCF,$normal,$WFH_TSV,$g_pu,$options);
			$depth=0;
			$mutant_depth=0;
		 }
		  #if($count >0) {exit;}			
			if($count % 100 == 0 ) {
				if($options->{'r'}){
					$log->debug("Normal_id:$normal: completed:$count PASS variants of $total_locations total merged variant sites");
				}
				else {
					$log->debug("Normal_id:$normal: completed:$count out of $total_locations total merged variant sites");
				}
			}
	}# foreach location
		
		if($options->{'m'}) {
		 	_write_results($augment_vcf,$aug_vcf_fh,$store_results,$input_bam_files,$aug_vcf_name); 
		}
		
	if ($vcf) {
		$vcf->close();
		close($WFH_VCF);
		close($WFH_TSV);
		$log->debug("Validating VCF file");
		my($outfile_gz,$outfile_tabix)=_compress_vcf("$outfile_name.vcf",$options);
                if ((-e $outfile_gz) && (-e $outfile_tabix)) {
		  unlink "$outfile_name.vcf" or $log->warn("Could not unlink".$outfile_name.'.vcf:'.$!);
		}
		validate("$outfile_gz");
	}
	
	$log->debug("Completed analysis for: $count locations");
	print $progress_fh "$outfile_name.vcf\n";
	return 1;
}



sub _compress_vcf {
  my ($annot_vcf,$options)=@_;
  my $annot_gz = $annot_vcf.'.gz';
  my $tmp_sorted_vcf=$options->{'o'}.'/tmp_sorted_vcf';
  my $command= 'vcf-sort '.$annot_vcf.' >'.$tmp_sorted_vcf;
  _run_external($command, 'vcf-sort', undef, 1, 1); # croakable, quiet, no data
  $command = 'bgzip -c '.$tmp_sorted_vcf.' > '.$annot_gz;
  _run_external($command, 'bgzip', undef, 1, 1); # croakable, quiet, no data

  $command = 'tabix -p vcf ';
  $command .= $annot_gz;
  _run_external($command, 'tabix', undef, 1, 1); # croakable, quiet, no data

  my $annot_tabix = "$annot_gz.tbi";
  croak "Tabix index does not appear to exist" unless(-e $annot_tabix);

  return ($annot_gz, $annot_tabix);
}

sub _compress_augmented_vcf {
  my ($annot_vcf)=@_;
  my $annot_gz = $annot_vcf.'.gz';
 
  my $command = 'bgzip -c '.$annot_vcf.' > '.$annot_gz;
  _run_external($command, 'bgzip', undef, 1, 1); # croakable, quiet, no data

  $command = 'tabix -p vcf ';
  $command .= $annot_gz;
  _run_external($command, 'tabix', undef, 1, 1); # croakable, quiet, no data

  my $annot_tabix = "$annot_gz.tbi";
  croak "Tabix index does not appear to exist" unless(-e $annot_tabix);

  return ($annot_gz, $annot_tabix);
}


sub _run_external {
        my ($command, $ext_prog, $no_croak, $quiet, $no_data, $FH, $binary) = @_;
        croak "Filehandle must be defined for binary output." if($binary && !$FH);
        return _run_external_core(q{-|}, $command, $ext_prog, $no_croak, $quiet, $no_data, $FH, $binary);
}

sub _run_external_core {
  my ($open_type, $command, $ext_prog, $no_croak, $quiet, $no_data, $FH, $binary) = @_;
 
  # ensure that commands containing pipes give appropriate errors
  $ENV{SHELL} = '/bin/bash'; # have to ensure bash is in use
  my $command_prefix = q{};
  $command_prefix = 'set -o pipefail; ' if($command =~ m/[|]/); # add pipefail to command if pipes are detected

  my (@prog_data, $tmp_fn);
        my $error_to_check;
  $log->warn(">>>>> $command");
  $command = $command_prefix.$command;
        if($EXECUTE_EXTERNAL == 1) {
          try {
      if(defined $FH && ref \$FH eq 'SCALAR') {
        $tmp_fn = $FH;
        undef $FH;
        open $FH, '>', $tmp_fn or croak "Failed to create $tmp_fn: $OS_ERROR";
      }
      my ($pid, $process);
      if($open_type eq q{2>&1}) {
        $pid = open $process, $command.' 2>&1 |' or croak 'Could not fork: '.$OS_ERROR;
      }
      else {
        $pid = open $process, q{-|}, $command or croak 'Could not fork: '.$OS_ERROR;
      }
      croak 'Failed to fork for: '.$command unless($pid);
      if($binary) {
        binmode $FH;
        binmode $process;
        my $buffer;
        my $buffer_max = 64*1024; # 64k
        while(read ($process, $buffer, $buffer_max) ){
          print $FH $buffer or croak "filehandle write failed $OS_ERROR";
        }
      }
      else {
        while (my $tmp = <$process>) {
          print $FH $tmp or croak "filehandle write failed $OS_ERROR" if($FH);
          $log->debug("<$ext_prog>\t$tmp") or $log->logcroak("logging write failed $OS_ERROR") unless($quiet);
          unless($no_data) {
            chomp $tmp;
            push @prog_data, $tmp ;
          }
        }
      }
      close $process or croak "Error closing pipe from $ext_prog";
    } catch {
      unless($no_croak) {
        my $message = q{This external process failed: }.$command;
      $message .= "\nERROR: $_";
        croak $message;
      } elsif($no_croak eq 'warn') {
        my $message = q{This external process failed: }.$command;
        $message .= "\nWARNING: $_";
        $log->warn($message);
      } elsif($no_croak eq 'check') {
        $error_to_check = $_;
      }
      # else just ignore the error
    }
    finally {
      close $FH or croak "Failed to close $tmp_fn: $OS_ERROR" if(defined $tmp_fn);
    };
        }
        $log->debug( "<<<<<") unless($quiet);

        if (defined $error_to_check) {
                $log->logcroak("_run_external should be called in list context if called with no_croak=check, returning an error value") unless wantarray;
                return (\@prog_data, $error_to_check);
        } else {
                return (\@prog_data);
        }
}


=head2 _write_results
Write augmented vcf file
Inputs
=over 2
=item input_dir
user provided input dir
=item normal_sample -normal sample hash
=item conn -database connection object
=back
=cut

sub _write_results {
 my ($augment_vcf,$aug_vcf_fh, $store_results,$input_bam_files,$aug_vcf_name)= @_;
		my $sample_counter=0;
			foreach my $sample (@$input_bam_files) {
				$sample_counter++;
				if($sample_counter > 1) {
					_write_final_vcf($augment_vcf->{$sample},$aug_vcf_fh,$sample,$store_results,$aug_vcf_name);
					$aug_vcf_fh->{$sample}->close();
				}
			}
 $log->debug("Completed writing VCF file");
}
 
=head2 _create_symlink_for_normal
create symlink for normal sample BAM file
Inputs
=over 2
=item input_dir
user provided input dir
=item normal_sample -normal sample hash
=item conn -database connection object
=back
=cut

sub _check_normal_bam {
	my ($input_dir,$normal_sample,$sample_group)=@_;
	my ($normal,$flag);

	my ($file_types)=_get_file_types();
	if(%$normal_sample) {
		foreach my $key (keys %$normal_sample) {
			$normal=$key;
		}
	}
	else {
		$normal=$sample_group;
	}
	if( -e "$input_dir/$normal.bam") {
		$log->debug("Using $normal as normal sample");
		$flag=1; 
		return $normal,$flag;
	}
	else {
	  $log->debug("Unable to access Normal BAM:$input_dir/$normal.bam"); 
	}
	
	return $normal,$flag;
}

=head2 _get_bam_object
create bam object using Bio::DB::Sam
Inputs
=over 2
=item bam_files sample names
=item input_dir
used defined input folder
=item genome
path of the genome file

=back
=cut

sub _get_bam_object {
	my ($bam_files,$options)=@_;
	my (%bam_objects,%bas_files);
	foreach my $file (@$bam_files) {
		my $sam = Bio::DB::Sam->new(-bam => "$options->{'d'}/$file.bam",
															-fasta =>$options->{'g'},
															-expand_flags => 1);
		$sam->max_pileup_cnt($MAX_PILEUP_DEPTH);
		$bam_objects{$file}=$sam;
		$bas_files{$file}="$options->{'d'}/$file.bam.bas";
	}
	\%bam_objects,\%bas_files;
}

=head2 _get_bam_header_data
get_bam_header_data -insert size and chr length
Inputs
=over 2
=item bam_objects - Bio::DB sam object
=back
=cut

sub _get_bam_header_data {
	my ($bam_objects,$bas_files)=@_;
  my ($chr_len,%bam_header_data);
  my $lib_size=0;
  foreach my $key (keys %$bam_objects) {	
  	my($mapped_length)=_get_read_length($bam_objects->{$key});
  	my $read_len=reduce{ $mapped_length->{$a} > $mapped_length->{$b} ? $a : $b } keys %$mapped_length;
   	my $header=$bam_objects->{$key}->header;
   	my $n_targets = $header->n_targets;
   	my $chr_names=$header->target_name;
   	foreach my $chr_no ((0..$n_targets)) {
	 		my $chr_name = $header->target_name->[$chr_no];
	 		my $chr_length = $header->target_len->[$chr_no];
	 		if(defined $chr_name && defined $chr_length)	{
	 			$bam_header_data{$key}{$chr_name}=$chr_length;
	 		}
  	}
  	#get insert size
  	# taken from brass code...
  	my $max=0;
  	my %sample_names = ();
  	foreach (split /\n/, $header->text) {
                next unless /^\@RG.*/;
                if(/.*\tMI:(?:Z:)?(\d+).*/){
                        $max = $1 if $1 > $max;
                        #added to store lib size for all samples..
                        $lib_size=$max if $max > $lib_size;
                }
                $sample_names{$1}++    if(/.*\tSM:([^\t]+).*/);
    }
        $log->warn("Multiple sample names detected ") if scalar keys %sample_names > 1;

        my @names = keys %sample_names;
    # case where MI median inser size tag is absent in BAM file [ use the insert size from bas file ]
    
  	if ($lib_size == 0 ) {$lib_size=$read_len*2};
  	
  	if( -e $bas_files->{$key} ) {
  		 $lib_size=_get_lib_size_from_bas($bas_files->{$key});
  	}
  	elsif ($lib_size == 0 )  {
  	  $lib_size=$read_len*2;
  	  $log->debug("No insert size tag in header or no BAS file found for $key using lib size (read_len x 2):$lib_size");
  	}
		$bam_header_data{$key}{'lib_size'}=$lib_size;
		if(defined $read_len) {
			$bam_header_data{$key}{'read_length'}=$read_len;
		}
	}
	
	\%bam_header_data,$lib_size;
}

=head2 _get_lib_size_from_bas
get library size from BAS file
Inputs
=over 2
=item bas file name
=back
=cut

sub _get_lib_size_from_bas {
  my ($bas_file)=@_;
	open(my $bas,'<',$bas_file) || $log->logcroak("unable to open file $!");
	my $index_mi=undef;
	my $index_sd=undef;
	my $lib_size=0;
	while (<$bas>) {
					my @array=split(/\t/, $_);
					if (defined $index_mi && $index_sd ) {
					 $lib_size=$array[$index_mi]+ ($array[$index_sd] * 2);
					 last;
					}
					if (!$index_mi && !$index_sd ) {
							$index_mi = first { $array[$_] eq $LIB_MEAN_INS_SIZE } 0 .. $#array;
							$index_sd = first { $array[$_] eq $LIB_SD } 0 .. $#array;
					}
	}
	close($bas);
	$lib_size;
}

=head2 _get_read_length
get read length per sample
Inputs
=over 2
=item sam - Bio::DB::Sam objects
=back
=cut

sub _get_read_length {
	my ($sam)=shift;
	my ($mapped_length,$read_counter);
	my $iterator  = $sam->features(-iterator=>1);
		while (my $a = $iterator->next_seq) { 
 			my $qseq=$a->qseq;
 			if(defined $qseq) {
				$read_counter++;
				my $len=length($qseq);
				$mapped_length->{$len}=$read_counter;
				if($read_counter > 100){
					return $mapped_length;
				}
			}
 		}
}

=head2 _write_augmented_header
write VCF header data
Inputs
=over 2
=item input_bam_files -analysed sample names 

=back
=cut

sub _write_augmented_header {
	my($augment_vcf,$options,$input_bam_files,$normal)=@_;
	my ($aug_vcf_fh,$aug_vcf_name);
	if($augment_vcf) {
		my ($vcf_filter,$vcf_info,$vcf_format)=_custom_header_lines($options);
		foreach my $sample (keys %$augment_vcf) { 
			my ($tmp_file)=(split '/', $augment_vcf->{$sample})[-1];
			$tmp_file=~s/(\.vcf|\.gz)//g;
			my $aug_vcf=$options->{'o'}.'/'.$tmp_file.$options->{'oe'};
			$aug_vcf_name->{$sample}=$aug_vcf;
			open(my $tmp_vcf,'>',$aug_vcf)|| $log->logcroak("unable to open file $!");
			$log->debug("Augmented vcf file:".$options->{'o'}."/$tmp_file.".$options->{'oe'});
			my $vcf_aug = Vcf->new(file => $augment_vcf->{$sample});
			$vcf_aug->parse_header();
			#to get all the FORMAT fileds in one group 
			my $tmp_header=$vcf_aug->get_header_line(key=>'FORMAT');
			$vcf_aug->remove_header_line(key=>'FORMAT');
			foreach my $line($tmp_header) {
				foreach my $format_line ( @$line) {
					foreach my $format_field (sort keys %$format_line) {
						$vcf_aug->add_header_line($format_line->{$format_field});
					}
				}
			}	
			$vcf_aug->add_header_line($vcf_format->{'MTR'});
			$vcf_aug->add_header_line($vcf_format->{'WTR'});
			$vcf_aug->add_header_line($vcf_format->{'AMB'});
			$vcf_aug->add_header_line($vcf_format->{'process_log'});

			print $tmp_vcf $vcf_aug->format_header();
			$vcf_aug->close();
			$aug_vcf_fh->{$sample}=$tmp_vcf;
		}
	}
    
    if($options->{'b'}) {
    	my @bed_header=qw(chr pos ref alt MTR_NORMAL WTR_NORMAL AMB_NORMAL MTR_TUMOUR WTR_TUMOUR AMB_TUMOUR);
    	foreach my $sample (@$input_bam_files) {
    		if ($sample ne $normal) {
    			open(my $tmp_bed,'>',"$options->{'o'}/$sample.augmented.bed");
    			$log->debug("Augmented bed file:".$options->{'o'}."/$sample.augmented.bed");
    			print $tmp_bed join("\t",@bed_header)."\n";
    			$aug_vcf_fh->{"$sample\_bed"}=$tmp_bed;
    			
    		}
    	}
    }
    $aug_vcf_fh,$aug_vcf_name;
}

# 
sub _trim_file_path{
	my ( $string ) = @_;
	my @bits = (split("/", $string));
	return pop @bits;
}



=head2 _write_vcf_header
write VCF header data
Inputs
=over 2
=item input_bam_files -analysed sample names 
=item info_tag_val -VCF INFO tag values
=item WFH -output file handler
=item genome -name of the genome file
=back
=cut


sub _write_vcf_header {
	my($input_bam_files,$info_tag_val,$normal,$vcf_file_status,$options,$tags,$normal_sample,$outfile_name)=@_;
	my $i=0;
	
	my ($vcf,$WFH_VCF,$WFH_TSV);
	
	# return if no VCF file found or augment only option is provided for indel data ...
	if((!%$normal_sample && !$options->{'b'}) || ($options->{'ao'})) {
			return $vcf,$WFH_VCF,$WFH_TSV;
	}
		
	$log->debug("VCF outfile:$outfile_name.vcf");
	$log->debug("TSV outfile:$outfile_name.tsv");
	open($WFH_VCF, '>',"$outfile_name.vcf");
	open($WFH_TSV, '>',"$outfile_name.tsv");

	$vcf=Vcf->new();
	my $genome_name=(split '/' ,$options->{'g'})[-1];
	my $script_name=(split '/' ,$0)[-1];
	$vcf->add_header_line({key=>'reference', value=>$genome_name}); 
	$vcf->add_header_line({key=>'source', value=>$script_name}); 
	$vcf->add_header_line({key=>'script_version', value=>$version});
	$vcf->add_header_line({key=>'Date', value=>scalar(localtime)});
	if($info_tag_val) {
		foreach my $hash_val(@$info_tag_val) {
			$vcf->add_header_line($hash_val);	
		}
	}
	else {
		$log->debug("No data in info filed");
	}
	
	foreach my $key (keys %$vcf_file_status) {
		
		if(($key ne $normal) && ($vcf_file_status->{$key} ne '1')) {
			$i++;
			my %temp=(key=>"SAMPLE",ID=>"NO_VCF_TUMOUR_$i", Description=>"NO_VCF_DATA_$i", SampleName=>$key);
			$vcf->add_header_line(\%temp);
		}
		
		if (!%$normal_sample && $key eq $normal) {
				my %temp=(key=>"SAMPLE",ID=>"NO_VCF_NORMAL", Description=>"NO_VCF_DATA", SampleName=>$key);
				$vcf->add_header_line(\%temp);
			}
	}
	
	$vcf->add_columns(@$input_bam_files);
	#To do add process log 
	print $WFH_VCF $vcf->format_header();	
	
	# writing results in tab separated format
	my $tmp_file="$options->{'o'}/temp.vcf";
	open(my $tmp_vcf,'>',$tmp_file);
	print $tmp_vcf $vcf->format_header();
	close($tmp_vcf);

	my($col_names,$header,$format_col)=_get_tab_sep_header($tmp_file);
	my $temp_cols=$col_names->{'cols'};
	
	foreach my $sample(@$input_bam_files){
		foreach my $tag_name(@$tags){
			push ($temp_cols,"$sample\_$tag_name");
		}
	}   
	print $WFH_TSV "@$header\n".join("\t",@$temp_cols)."\n";
	
		
	$vcf,$WFH_VCF,$WFH_TSV;
}
=head2 _get_original_results
get sample specific original values for INFO and FILTER fields
Inputs
=over 2
=item sample_names -sample names
=item data_for_all_samples -original VCF fileds for given sample
=item key_location -variant position
=item normal -normal sample name

=back
=cut

sub _get_original_results	{
	my ($sample_names,$data_for_all_samples,$key_location,$vcf_file_status,$no_vcf_counter,$options,$location)=@_;
	my (%original_vcf_info, %flag_val,$NFS,$DNFS,$old_info_val);
	#if($options->{'b'} && $key_location=~/BEDFILE/) {return;}
	#remove normal sample added at the beginning of the array
	my $max_rd=0;
	my $max_depth=0;
	foreach my $sample (@$sample_names) {
		next if $vcf_file_status->{$sample} eq '0';
		if(exists $data_for_all_samples->{$sample}{$key_location} ) {
			my $info_line=$data_for_all_samples->{$sample}{$key_location}->{'INFO'};
			my $filter_val=$data_for_all_samples->{$sample}{$key_location}->{'FILTER'};
			$max_rd=$data_for_all_samples->{$sample}{$key_location}->{'RD'};
			$max_depth=$max_rd if $max_rd > $max_depth;
			$flag_val{$sample}=$filter_val;
			$original_vcf_info{$sample}=$info_line;
		}
		else {
			$original_vcf_info{$sample}="NA";
			$flag_val{$sample}="NA";
		}
	}
	
		($NFS,$DNFS)=_get_NFS(\%flag_val,$location);
		$old_info_val =_get_INFO(\%original_vcf_info,$key_location,$DNFS,$no_vcf_counter);
	$old_info_val,$NFS,\%flag_val,$max_depth;
}
=head2 _get_NFS
get new filter status based on original filter status
Inputs
=over 2
=item flag_val -Original flag values in filter column
=item normal -normal sample name

=back
=cut


sub _get_NFS {
	my ($flag_val,$location)=@_;
	my $called=0;
	my $passed=0;
	my $not_called=0;
	my $called_n_passed=0;
	my $samples=0;
	my $NFS=0;
	
	foreach my $key (keys %$flag_val) {	
		if ($flag_val->{$key} eq "PASS") 	{$passed++;}
		elsif ($flag_val->{$key} eq "NA")	{$not_called++;}
		else{$called++;}
		$samples++;
		
	}
	$called_n_passed=$called+$passed;
		if ($location=~m/BEDFILE/) {
		return 'BD',"$samples:$called_n_passed:$passed";
	}
	
	# format as per the custom header
	#$vcf_filter{'filter1'}={(key=>'FILTER',ID=>'1',Description=>"New filter status 1=not called in any")}; # use of 0 is not allowed in FILTER column
	#$vcf_filter{'filter2'}={(key=>'FILTER',ID=>'2',Description=>"New filter status 2=called in any")};
	#$vcf_filter{'filter3'}={(key=>'FILTER',ID=>'3',Description=>"New filter status 3=passed  in any")};	
	#$vcf_filter{'filter4'}={(key=>'FILTER',ID=>'BD',Description=>"New filter status BD= location from bed file)};
	#0 - not called; 1 - called ;2 - passed ; 3 - called + passed
	#passed in all
	if(($passed == $samples) && $passed>0 ) { $NFS='PASS'; } #=2
	#called in all
	elsif(($called == $samples) && $passed>0) { $NFS=2; } #= 1
	#not called in all
	elsif(($not_called == $samples) && $not_called>0) { $NFS=1; } #=0
	# passed in all called samples
	elsif(($not_called + $passed == $samples) && $passed > 0 ) { $NFS='PASS'; } #=2
	#called + passed + not called
	#14 20924792 PD9179a PD9179b PD9179c PD9179d [ sample PD9179d doesn't contain the merged location]
	elsif( ( ($passed + $called + $not_called) == $samples ) && ($passed >0 && $called >0)) { $NFS=3; } #=3
	#called but not passed
	elsif( ( ($not_called + $called) == $samples ) && $called > 0) { $NFS=2; } #=1
	# if location is not present in any of the samples -- only valid for bed file., 
	else {$NFS='BD';} #=0
	
	#descriptive flag status added in the INFO field excludes normal sample
	my $DNFS="$samples:$called_n_passed:$passed";
	$NFS,$DNFS;
	
}
=head2 _get_INFO
parse info filed and 
Inputs
=over 2
=item original_vcf_info -data from INFO field
=item key_location -variant position
=item DNFS -Descriptive fiter status [ref header info]
=back
=cut


sub _get_INFO {
	my ($original_vcf_info,$key_location,$DNFS,$no_vcf_counter)=@_;
	my %old_info_val;
	foreach my $sample (keys %$original_vcf_info) {
		next if !defined $original_vcf_info->{$sample} || $original_vcf_info->{$sample} eq "NA";
		my $info=$original_vcf_info->{$sample};
		foreach my $tag (keys %$info) {
			$old_info_val{$tag}=$info->{$tag};
		}
	}
	$old_info_val{'NS'}=((split ':', $DNFS)[0]) + $no_vcf_counter;
	$old_info_val{'NC'}=(split ':', $DNFS)[1];
	$old_info_val{'NP'}=(split ':', $DNFS)[2];
	$old_info_val{'NA'}=(split ':', $DNFS)[0];
	
\%old_info_val;
		
}	
=head2 _get_region
get formatted position and incase of pindel data flag values for insertion del
Inputs
=over 2
=item location -variant position
=item variant_type - variant type [indel or SNP ]

=back
=cut

sub _get_region {
  my ($location,$variant_type) = @_;
  my ($chr,$start,$end,$ref,$alt,$insertion_flag,$del_flag,$g_pu);
  my (@data)=split(':',$location);
  $chr=$data[0];
  $start=$data[1];
  $end=$data[1];
  $ref=$data[2];
  $alt=$data[3];
  #different for indels...
  if($variant_type eq 'indel') {
		# if this is an insertion event no change in ref segment
		# if deletion event add length of del bases to start position to get end of deletion 
		# data for this analysis is not taken from info field to accomodate the custom bed locations...
		$insertion_flag=1;
		$del_flag=0;
		if(defined $ref && length($ref) > 1 ) { # if deletion 
			$end=$start+(length($ref)); # no need to do +1 ref base already included in the deletion seq
			$insertion_flag=0;
			$del_flag=1;
		}
		#case of ID --> complex indel
		if(defined $ref && length($alt) > 1) { $insertion_flag=1;}
		
  }

$g_pu={'chr'=>$chr, 'start' => $start, 'end' => $end, 'ref_seq' => $ref, 'alt_seq' => $alt,
				'region' => "$chr:$start-$end", 'ins_flag' => $insertion_flag, 'del_flag' => $del_flag, 
			};
$g_pu

}
=head2 _populate_hash
populate has for chr start and end
Inputs
=over 2
=item g_pu - has containing sample specific info
=item sample - sample name
=item bam_header_data - sample specific information on read length and lib size
=back
=cut

sub _populate_hash {
  my ($g_pu,$sample,$bam_header_data) = @_;
  
  $g_pu->{'ref_p'} = 0;
	$g_pu->{'ref_n'} = 0;
	$g_pu->{'alt_p'} = 0;
	$g_pu->{'alt_n'} = 0;
	
	$g_pu->{'FAZ'} = 0;
	$g_pu->{'FCZ'} = 0;
	$g_pu->{'FGZ'} = 0;
	$g_pu->{'FTZ'} = 0;
	$g_pu->{'RAZ'} = 0;
	$g_pu->{'RCZ'} = 0;
	$g_pu->{'RGZ'} = 0;
	$g_pu->{'RTZ'} = 0;
	 
	$g_pu->{'amb'} = 0;  
	$g_pu->{'sample'} = $sample; 
	if(exists $bam_header_data->{$sample}{'read_length'}){
		$g_pu->{'read_length'}=$bam_header_data->{$sample}{'read_length'};
		
	}
	else {
		$g_pu->{'read_length'}=100;
	}
	if(exists $bam_header_data->{$sample}{'lib_size'}){
		$g_pu->{'lib_size'}=$bam_header_data->{$sample}{'lib_size'};
	}
	else {
		$g_pu->{'lib_size'}=100;
	} 
	# exonerate score is 5 per base , we allow max 4 mismatches * 9 = 36 score units, 4 bases * 5 for readlength = 20 score units to be safe side added another 14 score units
	$g_pu->{'exonerate_score_cutoff'} = (($g_pu->{'read_length'}) * $EXONERATE_SCORE_MULTIPLIER) - $EXONERATE_SCORE_FACTOR;	                 
  $g_pu;
}

=head2 _get_range
get left and right span from  indel position
Inputs
=over 2
=item bam_header_data - sample specific information on read length, chr length and lib size 
=item sample -sample name
=item g_pu - has containing sample specific info
=item tabix_hdr - tabix object created from UCSC high depth region bed file
=back
=cut

sub _get_range {
  my($bam_header,$sample,$g_pu,$tabix_hdr,$max_depth)=@_;
  my ($left_pos,$right_pos,$chr_len,$spanned_region);
  my $lib_size=$bam_header->{$sample}{'lib_size'};
  my ($hdr_flag)=check_hdr_overlap($g_pu->{'chr'},$g_pu->{'start'},$g_pu->{'end'},$tabix_hdr);
  my $spanning_seq_denom=2;
  #if location is in high depth region and has depth >1000 then hdr_flag is true
  if($hdr_flag && $max_depth > 1000){$spanning_seq_denom=4;}
  else {$hdr_flag=0;}
  $chr_len=$bam_header->{$sample}{$g_pu->{'chr'}};
  if(defined $lib_size && defined $chr_len) {
  	$spanned_region = round(($lib_size *  $INSERT_SIZE_FACTOR )/$spanning_seq_denom);
  	#$spanned_region=round($spanned_region);
  	# 
		if(($spanned_region < $g_pu->{'start'}) && (($chr_len - $g_pu->{'end'}) > $spanned_region)) {
			$left_pos=$g_pu->{'start'} - $spanned_region;
			$right_pos=$g_pu->{'end'} + $spanned_region;
		}
	}

	$g_pu->{'pos_5p'}=$left_pos;
	$g_pu->{'pos_3p'}=$right_pos;
	$g_pu->{'hdr'}=$hdr_flag;
	# to get exact location of variant base and avoid borderline matching of a  
	# added padding to reference posotion 
	$g_pu->{'ref_pos_5p'}=($g_pu->{'start'} - $g_pu->{'pos_5p'}) + 1 ;
	if($g_pu->{'ins_flag'} && !$g_pu->{'del_flag'}){
		$g_pu->{'ref_pos_3p'}=$g_pu->{'ref_pos_5p'} + ( $g_pu->{'end'} - $g_pu->{'start'} );
	}
	else{
		$g_pu->{'ref_pos_3p'}=$g_pu->{'ref_pos_5p'} + ( $g_pu->{'end'} - $g_pu->{'start'} ) - 1;
	}
	$g_pu->{'alt_pos_3p'}=$g_pu->{'ref_pos_5p'} + length( $g_pu->{'alt_seq'}) -1;

	$g_pu;
}
=head2 _get_alt_seq
get alternative reference seq
Inputs
=over 2
=item bam_objects - Bio::DB sam object
=item g_pu - has containing sample specific info
=back
=cut
sub _get_alt_seq {
	my ($bam_objects,$g_pu)=@_;
	#<---1.5*insert_size---[alt_left]|indel|[alt_right]------1.5*insert_size--->	
	my $tmp_start=$g_pu->{'start'};
	my $tmp_end=$g_pu->{'end'};
	#insertion			
	
	if($g_pu->{'ins_flag'} == 1  && $g_pu->{'del_flag'} == 0) {
		$tmp_start= $g_pu->{'start'} - 1;
		$tmp_end= $g_pu->{'end'} + 1;
	}
	#complex indel
	if($g_pu->{'del_flag'} == 1  && $g_pu->{'ins_flag'} == 1) {
		$tmp_start=$g_pu->{'start'} - 1;
	}
	
	my ($alt_left_seq)=_get_dna_segment($bam_objects,$g_pu->{'chr'},$g_pu->{'pos_5p'},$tmp_start);
	my ($alt_right_seq)=_get_dna_segment($bam_objects,$g_pu->{'chr'},$tmp_end,$g_pu->{'pos_3p'});
	#indel[insertion] add alt between two segments...
	my $reconstructed_alt_seq;
	
	if(defined $g_pu->{'ins_flag'} && $g_pu->{'ins_flag'} == 1) {	
		$reconstructed_alt_seq=$alt_left_seq.$g_pu->{'alt_seq'}.$alt_right_seq;			
	}
	#indel[deletion] join two segments...
	else {
		$reconstructed_alt_seq=$alt_left_seq.$alt_right_seq;
	}
	
	$reconstructed_alt_seq;

}
=head2 _get_dna_segment
get reference dna 
Inputs
=over 2
=item sam_object - Bio::DB sam object
=item chr - chromosome number
=item start - start position
=item end - end position
=back
=cut

sub _get_dna_segment {
	my ($bam_object,$chr,$start,$end)=@_;
	my ($segment)=$bam_object->segment($chr,$start,$end);
	$segment->dna;
}
=head2 _get_pileup
get pileup output for given location
Inputs
=over 2
=item sam_object - Bio::DB sam object
=item g_pu - has containing sample specific info
=back
=cut


sub _get_pileup {
	my ($bam_object,$g_pu)=@_;
		$bam_object->fast_pileup($g_pu->{'region'}, sub {
									my ($seqid, $pos, $pu) = @_;
									return if($pos != $g_pu->{'start'});
									my $refbase = $bam_object->segment($seqid,$pos,$pos)->dna;
									foreach my $p (@{$pu}) {
										next if($p->is_del || $p->is_refskip);
										my $a = $p->alignment;
										
										my $flags = $a->flag;
										# & bitwise comparison
										## Ignore reads if they match the following flags:
										#Brass/ReadSelection.pm
		
										next if $flags & $NOT_PRIMARY_ALIGN;
										next if $flags & $VENDER_FAIL;
										next if $flags & $UNMAPPED;
									
										#next if($g_pu_mapq && $a->qual < $g_pu_mapq);
										#if($g_pu_baseq) {
										#	my $fa = Bio::DB::Bam::AlignWrapper->new($a, $bam_object);
										#	my $qual = ($fa->qscore)[$p->qpos];
											#next if($qual <= $g_pu_baseq);
										#}
										# get the base at this pos
										#my $refbase = $bam_object->segment($seqid,$pos,$pos)->dna;
										my $qbase  = substr($a->qseq, $p->qpos, 1);
										my $strand = $a->strand;
										next if $qbase =~ /[nN]/; #in case of insertion ....
										#$g_pu->{'depth'}++; # commented as for paired end it is calculated twice
										my $key;
										if($refbase eq $qbase && $strand > 0) {
											$g_pu->{'ref_p'}++;
											$key='F'.$qbase.'Z';
										}
										elsif($refbase eq $qbase && $strand < 0) {
											$g_pu->{'ref_n'}++;
											$key='R'.$qbase.'Z';
										}
										elsif(($g_pu->{'alt_seq'} eq $qbase) && $strand > 0) {
											$g_pu->{'alt_p'}++;
											$key='F'.$qbase.'Z';
										}
										elsif(($g_pu->{'alt_seq'} eq $qbase) && $strand < 0) {
											$g_pu->{'alt_n'}++;
											$key='R'.$qbase.'Z';
										}
										elsif ($strand > 0 ) {
											$key='F'.$qbase.'Z';
										}
										elsif ($strand < 0 ) {
											$key='R'.$qbase.'Z';
										}
										
										$g_pu->{$key}++;

									}
		});
	
	$g_pu;	
}
=head2 _fetch_features
get reads in an given location using Bio:Db:Sam fetch method
Inputs
=over 2
=item sam_object - Bio::DB sam object
=item g_pu - hash containing sample specific info
=item Reads_FH - file handler to store reads for a given location
=back
=cut


sub _fetch_features {
	my ($sam_object,$g_pu,$Reads_FH,$options)=@_;
	
	if(($g_pu->{'end'} - $g_pu->{'start'}) < $g_pu->{'lib_size'}){
		_fetch_reads($sam_object, "$g_pu->{'region'}",$Reads_FH);
		$g_pu->{'long_indel'}=0;
	}
	else {
		_fetch_reads($sam_object, "$g_pu->{'chr'}:$g_pu->{'start'}-$g_pu->{'start'}",$Reads_FH);
		_fetch_reads($sam_object, "$g_pu->{'chr'}:$g_pu->{'end'}-$g_pu->{'end'}",$Reads_FH);
		#indel longer than library size set flag on 
		$g_pu->{'long_indel'}=1;
		}
	#only get unmapped reads if not in HDR region and while analysing only PASSED varinat
	if(!$g_pu->{'hdr'} && $options->{'r'}){
		_fetch_unmapped_reads($sam_object,"$g_pu->{'chr'}:$g_pu->{'pos_5p'}-$g_pu->{'pos_3p'}",$Reads_FH);
	}
close($Reads_FH);
$g_pu;
}



=head2 _fetch_reads
fetch reads

Inputs
=over 2
=item sam_object - Bio::DB sam object
=item region - chr:start-stop format region info to get reads
=item Reads_FH - temp file handler to store reads
=back
=cut

sub _fetch_reads {
my ($sam_object,$region,$Reads_FH)=@_;
my ($mate_info,%mapped_length);
my $read_counter=0;
		$sam_object->fetch($region, sub {
		my $a = shift;
		my $paired=0;
		my $flags = $a->flag;
		# & bitwise comparison
		## Ignore reads if they match the following flags:
		#Brass/ReadSelection.pm
		
		return if $flags & $NOT_PRIMARY_ALIGN;
		return if $flags & $VENDER_FAIL;
		return if $flags & $UNMAPPED;
	  #return if $flags & $DUP_READ;
		return if $flags & $SUPP_ALIGNMENT;
		
		#if $flags & $READ_PAIRED;
		#my $cigar  = $a->cigar_str;;
		my $mseqid = $a->mate_seq_id;
		my $seqid = $a->seq_id;
		#target gives read seq as it comes from sequencing machine i.e softclipped bases included
		my $qseq = $a->target->dna();
		return if $qseq =~m/[nN]/;
		my $name=$a->display_name;
		#my $strand = $a->strand;
		my $mstart = $a->mate_start;
		my $start = $a->start;
		$read_counter++;
		print  $Reads_FH ">$name\_$read_counter\n$qseq\n";
	# fetch mate only if on another chromosome
		if(defined $mseqid and defined $seqid and ($seqid ne $mseqid)) {
		#if(defined $mseqid and defined $seqid ) {
			$mate_info->{$name}="$mseqid:$mstart-$mstart";
		}
				
	});
	
	#added separately as it was pulling only one read 
	foreach my $key (keys %$mate_info)
	{
		_fetch_mate_seq($sam_object,$mate_info->{$key},$key, $Reads_FH);
	}
}
=head2 _fetch_mate_seq
get mate sequence
Inputs
=over 2
=item sam_object - Bio::DB sam object
=item region - postion to get reads
=item readname - mate readname
=item Reads_FH - file handler to store read sequence
=back
=cut

sub _fetch_mate_seq {
	my ($sam_object,$region,$readname,$Reads_FH)=@_;
	my ($read,$mate_seq);
	my $callback= sub {
		my $a = shift; 
		my $flags = $a->flag;
		return if $flags & $NOT_PRIMARY_ALIGN;
		return if $flags & $VENDER_FAIL;
	  #return if $flags & $DUP_READ;
	  return if $flags & $UNMAPPED;
		return if $flags & $SUPP_ALIGNMENT;
		if ($readname eq $a->display_name) {
			my $tmp_seq=$a->target->dna();
			return if $tmp_seq=~m/[nN]/;
			$read=$a->display_name;
			$mate_seq=$tmp_seq;
			return;
		}
	};
	
	$sam_object->fetch($region,$callback);
	if($read){
		print  $Reads_FH ">$read\_0\n$mate_seq\n";
	}
}

=head2 _fetch_mate_unmapped_reads
fetch reads whose mate is unmapped
Inputs
=over 2
=item sam_object - Bio::DB sam object
=item region - chr:start-stop format region info to get reads
=item Reads_FH - temp file handler to store reads
=back
=cut

sub _fetch_unmapped_reads {
my ($sam_object,$region,$Reads_FH)=@_;
my ($mate_info,%mapped_length);
my $read_counter=0;
		$sam_object->fetch($region, sub {
		my $a = shift;
		my $paired=0;
		my $flags = $a->flag;
		# & bitwise comparison
		## Ignore reads if they match the following flags:
		#Brass/ReadSelection.pm
		return if $flags & $NOT_PRIMARY_ALIGN;
		return if $flags & $VENDER_FAIL;
		#return if $flags & $DUP_READ;
		return if $flags & $SUPP_ALIGNMENT;
		# only consider reads from wider range where mate is unmapped 
		if ($flags & $UNMAPPED) {
			my $qseq = $a->target->dna();
			return if $qseq=~m/[nN]/;
			my $mseqid = $a->mate_seq_id;
			my $seqid = $a->seq_id;
			#target gives read seq as it comes from sequencing machine
			
			my $name=$a->display_name;
			#my $strand = $a->strand;
			my $mstart = $a->mate_start;
			my $start = $a->start;
			$read_counter++;
			print  $Reads_FH ">$name\_$read_counter\n$qseq\n";
			# fetch mate only if on another chromosome
			if(defined $mseqid and defined $seqid and ($seqid ne $mseqid)) {
					$mate_info->{$name}="$mseqid:$mstart-$mstart";
			}
	}
				
	});
	
	#added separately as it was pulling only one read 
	foreach my $key (keys %$mate_info)
	{
		_fetch_mate_seq($sam_object,$mate_info->{$key},$key, $Reads_FH);
	}
}

=head2 _do_exonerate
parse exonerate output 
Inputs
=over 2
=item ref_seq —reference sequence and alt sequences in fasta format
=item temp_read_file  —temp read file in fasta format
=item g_pu -- stores relative positions of variant region  
=back
=cut

sub _do_exonerate {
	my($ref_seq_file,$temp_read_file,$g_pu,$alt_file,$ref_file)=@_;
	my ($ref_count_p, $ref_count_n, $alt_count_p, $alt_count_n, $read_track_alt, $read_track_ref, $amb_reads);
	
	#print $new_interval "Before\t".$g_pu->{'alt_seq'}."\t".$g_pu->{'region'}."\t".$g_pu->{'sample'}."\t".$g_pu->{'ins_flag'}."\t".$g_pu->{'ref_pos_5p'}."\t".$g_pu->{'ref_pos_3p'}."\t".$g_pu->{'alt_pos_3p'}."\t";
	 
	#($g_pu)=_get_ref_5p_pos($alt_file,$ref_file,$g_pu);
	
	#print $new_interval "After\t".$g_pu->{'ins_flag'}."\t".$g_pu->{'ref_pos_5p'}."\t".$g_pu->{'ref_pos_3p'}."\t".$g_pu->{'alt_pos_3p'}."\n";

	
	# -E | --exhaustive <boolean>
  #Specify whether or not exhaustive alignment should be used.  By default, this is FALSE, and alignment heuristics will be used.  If it is set to TRUE, an exhaus‐
  #tive alignment will be calculated.  This requires quadratic time, and will be much, much slower, but will provide the optimal result for the given model. 
  #-S | --subopt <boolean>
  # using exhaustive OFF as it is fast and gives identical answer 
  
  my $cmd="exonerate -E 0 -S 0".
	" --score $g_pu->{'exonerate_score_cutoff'} --percent 95 --fsmmemory 12000 --verbose 0 --showalignment no  --wordjump 3".
	" --querytype dna --targettype dna --query $temp_read_file  --target $ref_seq_file".
	" --showvulgar 0 --bestn 1 --ryo '%qi %ti %qal %tS %tab %tae %qS\n' ";
	#for testing only
	#if($test_mode)
	#{
	#my $cmd2="exonerate -E 0 -S 0".
	#	" --score $g_pu->{'exonerate_score_cutoff'} --percent 95 --fsmmemory 12000 --verbose 0 --showalignment yes --wordjump 3".
	#	" --querytype dna --targettype dna --query $temp_read_file   --target $ref_seq_file".
	#	" --showvulgar 0 --bestn 1 --ryo '%qi %ti %qal %tS %tab %tae %qS\n' ";
	#	my ($exonerate_output1, $stderr1, $exit1) = capture {system("$cmd2")};
	#open (my $tfh1, '>',"exonerate_results_Alignment.out");
	#print $tfh1 $exonerate_output1;
	#my ($exonerate_output2, $stderr2, $exit2) = capture {system("$cmd")};
	
 # }
	
	my ($exonerate_output, $stderr, $exit) = capture {system("$cmd")};
	if ($exit) { $log->logcroak("exonerate log: EXIT:$exit EROOR:$stderr CMD:$cmd"); }
	#----- parse exonerate output ------
	foreach my $line((split("\n", $exonerate_output))) {
		my ($read,$target,$match_len,$t_strand,$t_start,$t_end,$q_strand)=(split ' ', $line);
		if ($match_len < ($g_pu->{'read_length'} - $READ_LENGTH_CUTOFF)) {
		 next;
		}
		my $strand=$t_strand;
		#<--5p--|*******|--3p-->
		my $temp_start=$t_start;
		my $org_read=$read;
		$read=~s/_\d+$//g;
		if($strand eq '-') { $t_start=$t_end ; $t_end=$temp_start;}
		if( $target eq 'ref') {	
			# ref_pos stores the varinat interval relative to subset created using gnomic seq	
			if( ($t_start < $g_pu->{'ref_pos_5p'} &&  $t_end >$g_pu->{'ref_pos_5p'}) || ($t_start < $g_pu->{'ref_pos_3p'} &&  $t_end >$g_pu->{'ref_pos_3p'}) ) 
			{
				$read_track_ref->{$org_read}++;
				if($strand eq '+') {
				#store diff to check the distance of variant pos from either end of the read
					$ref_count_p->{$read} = abs( ($g_pu->{'ref_pos_5p'} - $t_start) - ($t_end - $g_pu->{'ref_pos_3p'}) )
				} 
				else {
					$ref_count_n->{$read} = abs( ($g_pu->{'ref_pos_5p'} - $t_start) - ($t_end - $g_pu->{'ref_pos_3p'}) );
				}
			}			
		}
		# checks overlap...
		elsif( ($t_start < $g_pu->{'ref_pos_5p'} &&  $t_end >$g_pu->{'ref_pos_5p'} ) || ($t_start < $g_pu->{'alt_pos_3p'} &&  $t_end >$g_pu->{'alt_pos_3p'}) )  
		{
			$read_track_alt->{$org_read}++;
			if($strand eq '+') {
				$alt_count_p->{$read} = abs( ($g_pu->{'ref_pos_5p'} - $t_start) - ($t_end - $g_pu->{'alt_pos_3p'}) );
			} 
			else {
				$alt_count_n->{$read} = abs( ($g_pu->{'ref_pos_5p'} - $t_start) - ($t_end - $g_pu->{'alt_pos_3p'}) );
			}
		}
	}	
	
$g_pu=_cleanup_read_ambiguities($g_pu,$read_track_alt,$read_track_ref, $alt_count_p,$alt_count_n,$ref_count_p,$ref_count_n); 


#print_hash($g_pu);

return $g_pu;


}


=head2 _get_ref_5p_pos
get updated 5p positions
Inputs
=over 2
=item g_pu -- stores relative positions of variant region  
=item alt_file —alt seq
=item ref_file -ref seq
=back
=cut

sub _get_ref_5p_pos {
	my ($ref_seq,$reconstructed_alt_seq,$g_pu) = @_;		
			my $new_pos;
			my $exclusive_OR=$ref_seq^$reconstructed_alt_seq;
			
			if($exclusive_OR =~ /[^\0]/g) {
				$new_pos=$-[0]; #gives offset of the beginning of last successful match
				$g_pu->{'new_pos'}=$new_pos;
			}	
		if( ($g_pu->{'ins_flag'}) && ($new_pos != $g_pu->{'ref_pos_5p'}) ){
				my $insert_length = ($g_pu->{'alt_pos_3p'} - $g_pu->{'ref_pos_5p'});
				# get run over after insert string due to match with reference bases
				$g_pu->{'insert_mod'}=($new_pos - $g_pu->{'ref_pos_5p'}) % $insert_length;
				$g_pu->{'alt_pos_3p'}=$new_pos + ($insert_length) - $g_pu->{'insert_mod'};
				$g_pu->{'ref_pos_5p'}=$new_pos;
				$g_pu->{'ref_pos_3p'}=$new_pos;
		}
		$g_pu;
}

=head2 _cleaup_read_ambiguities
cleaup of ambiguous reads and read of a same read pair mapping on +ve and/or -ve strand at same location
Inputs
=over 2
=item g_pu -- stores relative positions of variant region  
=item read_track_alt —hash stroring reads mapping on alt base
=item read_track_ref  —hash stroring reads mapping  on ref base
=item  alt_count_p  —hash stroring +ve reads mapping  on alt base
=item  alt_count_n  —hash stroring -ve reads mapping  on alt base
=item  ref_count_p  —hash stroring +ve reads mapping  on ref base
=item  ref_count_n  —hash stroring -ve reads mapping  on ref base
=back
=cut

sub _cleanup_read_ambiguities {
	my ($g_pu,$read_track_alt,$read_track_ref,$alt_count_p,$alt_count_n,$ref_count_p,$ref_count_n)=@_;
	my $amb_reads;
	foreach my $read (sort keys %$read_track_alt) {
			if(exists $read_track_ref->{$read}) {
				
				$read=~s/_\d+$//g;
				delete $ref_count_n->{$read} if $ref_count_n->{$read};
				delete $ref_count_p->{$read} if $ref_count_p->{$read};
				delete $alt_count_n->{$read} if $alt_count_n->{$read};
				delete $alt_count_p->{$read} if $alt_count_p->{$read};
				$amb_reads->{$read}++;
			}
		}	
	
	# check if read pairs maps at same location , consider only one read from pair which is aligned properly to variant site 
	foreach my $read (sort keys %$ref_count_p) {
		# if variant position on +ve strand is towards middle of the read then remove read from -ve strand 
		if ($ref_count_n->{$read}) {
			if($ref_count_p->{$read} < $ref_count_n->{$read}) {
				delete $ref_count_n->{$read};
			}
			else {
				delete $ref_count_p->{$read};
			}
		}
	}
	foreach my $read (sort keys %$alt_count_p) {
		# if variant position on +ve strand is towards middle of the read then remove read from -ve strand 
		if ($alt_count_n->{$read}) {
			if($alt_count_p->{$read} < $alt_count_n->{$read}) {
				delete $alt_count_n->{$read};
			}
			else {
				delete $alt_count_p->{$read};
			}
		}
	}

		if($ref_count_p) { $g_pu -> {'ref_p'}=keys %$ref_count_p; }
		if($ref_count_n) { $g_pu -> {'ref_n'}=keys %$ref_count_n; }
		if($alt_count_p) { $g_pu -> {'alt_p'}=keys %$alt_count_p; }
		if($alt_count_n) { $g_pu -> {'alt_n'}=keys %$alt_count_n; }
		if($amb_reads)	 { $g_pu -> {'amb'}=keys %$amb_reads; }
		
		#print Dumper $amb_reads;
		
return $g_pu;

}


=head2 _format_pileup_line
format pileup/ exonerate results as per VCF specifications
=over 2
=item original_flag -orignal flag values
=item g_pu - hash containing results and sample specific info for give location
=back
=cut

sub _format_pileup_line {
	my ($original_flag,$g_pu,$options)=@_;
	my $VCF_OFS;	
	my $pileup_results;
	if(defined $original_flag->{$g_pu->{'sample'}}) {
		$VCF_OFS=$original_flag->{$g_pu->{'sample'}};
	}
	else {
		$VCF_OFS='NA';
	}
	
	my $MTR = $g_pu->{'alt_p'} + $g_pu->{'alt_n'};
	my $WTR = $g_pu->{'ref_p'} + $g_pu->{'ref_n'};
	my $DEP = $MTR + $WTR + $g_pu->{'amb'};
	
	## determine read direction
	my $MDR =0;
	# only +ve reads 
	if($g_pu->{'alt_p'} > 0 && $g_pu->{'alt_n'} == 0 ) 
	{ $MDR=1; }
	# only -ve reads
	elsif($g_pu->{'alt_p'} == 0 && $g_pu->{'alt_n'} > 0 ) 
	{ $MDR=2; }
	# +ve & -ve
	elsif($g_pu->{'alt_p'} > 0 && $g_pu->{'alt_n'} > 0 )
	 { $MDR=3; }
 
	my $WDR=0;
	# only +ve
	if($g_pu->{'ref_p'} > 0 && $g_pu->{'ref_n'} == 0 )
	{ $WDR=1; }
	# only -ve
	elsif($g_pu->{'ref_p'} == 0 && $g_pu->{'ref_n'} > 0 ) 
	{ $WDR=2; }
	# +ve & -ve
	elsif($g_pu->{'ref_p'} > 0 && $g_pu->{'ref_n'} > 0 ) 
	{ $WDR=3; }
	
	if($options->{'a'} ne 'indel') {
	$pileup_results={ 
										'FAZ' =>$g_pu->{'FAZ'},
										'FCZ' =>$g_pu->{'FCZ'},
										'FGZ' =>$g_pu->{'FGZ'},
										'FTZ' =>$g_pu->{'FTZ'},
										'RAZ' =>$g_pu->{'RAZ'},
										'RCZ' =>$g_pu->{'RCZ'},
										'RGZ' =>$g_pu->{'RGZ'},
										'RTZ' =>$g_pu->{'RTZ'},
										'MTR'	=>$MTR,
										'WTR'	=>$WTR,
										'DEP'	=>$DEP,
										'MDR'	=>$MDR,
										'WDR'	=>$WDR,
										'OFS'	=>$VCF_OFS,
										};
	}
	
	else {
	$pileup_results={ 'MTR'=>$MTR,
										'WTR'=>$WTR,
										'DEP'=>$DEP,
										'MDR'=>$MDR,
										'WDR'=>$WDR,
										'OFS'=>$VCF_OFS,
										'AMB'=>$g_pu->{'amb'}
										};
	}
		
$pileup_results;	
}

=head2 _augment_vcf_file
Create a VCF object for a single vcf line by restricting it to a variant region
augment the resultant line with exonerate based variant read fraction
Inputs
=over 2
=item vcf_file --vcf file to augment
=item g_pu --hash storing information for a variant location
=item aug_vcf_fh  --hash stroring filehandler for given sample
=back
=cut

sub _augment_vcf_file {
	my ($vcf_file,$g_pu,$aug_vcf_fh, $location_val,$add_no_pass,$sample)=@_;
    
	  if ($add_no_pass) {
	  my $vcf = Vcf->new(file => $vcf_file, region => "$g_pu->{'chr'}:$g_pu->{'start'}-$g_pu->{'start'}" );
		$vcf->parse_header();	
		my $x = $vcf->next_data_hash();
		$vcf->add_format_field($x,'MTR');
		$vcf->add_format_field($x,'WTR');
		$vcf->add_format_field($x,'AMB');
		
		$$x{gtypes}{'TUMOUR'}{'MTR'}=".";
		$$x{gtypes}{'TUMOUR'}{'WTR'}=".";
		$$x{gtypes}{'TUMOUR'}{'AMB'}=".";
		
		$$x{gtypes}{'NORMAL'}{'MTR'}=".";
		$$x{gtypes}{'NORMAL'}{'WTR'}=".";
		$$x{gtypes}{'NORMAL'}{'AMB'}=".";
	  $aug_vcf_fh->{$sample}->print($vcf->format_line($x));
		$vcf->close();
		return;
	  }
	  
		my $MTR = $g_pu->{'alt_p'} + $g_pu->{'alt_n'};
		my $WTR = $g_pu->{'ref_p'} + $g_pu->{'ref_n'};
	
		if ($location_val=~/BEDFILE/) {
			my $bed_line="$g_pu->{'chr'}\t$g_pu->{'start'}\t$g_pu->{'ref_seq'}\t$g_pu->{'alt_seq'}".
						"\t$g_pu->{'normal_MTR'}\t$g_pu->{'normal_WTR'}\t$g_pu->{'normal_AMB'}\t".
						$MTR."\t".$WTR."\t".$g_pu->{'amb'}."\n";			
			$aug_vcf_fh->{"$g_pu->{'sample'}\_bed"}->print($bed_line);
			return;
		}
	  my $vcf = Vcf->new(file => $vcf_file, region => "$g_pu->{'chr'}:$g_pu->{'start'}-$g_pu->{'start'}" );
		$vcf->parse_header();	
		my $x = $vcf->next_data_hash();
		$vcf->add_format_field($x,'MTR');
		$vcf->add_format_field($x,'WTR');
		$vcf->add_format_field($x,'AMB');
		
		$$x{gtypes}{'TUMOUR'}{'MTR'}=$MTR;
		$$x{gtypes}{'TUMOUR'}{'WTR'}=$WTR;
		$$x{gtypes}{'TUMOUR'}{'AMB'}=$g_pu->{'amb'};
		
		$$x{gtypes}{'NORMAL'}{'MTR'}=$g_pu->{'normal_MTR'};
		$$x{gtypes}{'NORMAL'}{'WTR'}=$g_pu->{'normal_WTR'};
		$$x{gtypes}{'NORMAL'}{'AMB'}=$g_pu->{'normal_AMB'};
		
		# print data line to respective filehandler : following notation of accessing the filehandler 
		# from array or hash is called 'indirect object' notation which can be done by two methods
		# Method 1: print {$aug_vcf_fh->{$g_pu->{'sample'}}} $vcf->format_line($x);
		# Method 2: $aug_vcf_fh->{$g_pu->{'sample'}}->print($vcf->format_line($x));
		$aug_vcf_fh->{$g_pu->{'sample'}}->print($vcf->format_line($x));
		$vcf->close();
}

sub _store_results {
	my ($store_results,$g_pu,$add_no_pass,$sample,$location_val,$location)=@_;
	my $results = {'tMTR'=> '.', 
							 'tWTR'=> '.',
							 'tAMB'=> '.',
							 'nMTR'=> '.',
							 'nWTR'=> '.',
							 'nAMB'=> '.',								  
		};	
	if ($add_no_pass) {
		$store_results->{$sample}{$location}=$results;
		return $store_results;
	}
	
	
  my $MTR = $g_pu->{'alt_p'} + $g_pu->{'alt_n'};
	my $WTR = $g_pu->{'ref_p'} + $g_pu->{'ref_n'};

	if ($location_val=~/BEDFILE/) {
	    my $bed_line="$g_pu->{'chr'}\t$g_pu->{'start'}\t$g_pu->{'ref_seq'}\t$g_pu->{'alt_seq'}".
						"\t$g_pu->{'normal_MTR'}\t$g_pu->{'normal_WTR'}\t$g_pu->{'normal_AMB'}\t".
						$MTR."\t".$WTR."\t".$g_pu->{'amb'}."\n";
			$store_results->{"$sample\_bed"}{$location}=$bed_line;		
			return $store_results;
		}
	$results->{'tMTR'}=$MTR;
	$results->{'tWTR'}=$WTR;
	$results->{'tAMB'}=$g_pu->{'amb'};
	
	$results->{'nMTR'}=$g_pu->{'normal_MTR'};
	$results->{'nWTR'}=$g_pu->{'normal_WTR'};
	$results->{'nAMB'}=$g_pu->{'normal_AMB'};
	
  $store_results->{$sample}{$location}=$results;
  return $store_results; 
 
}

sub _write_final_vcf {
	  my ($vcf_file,$aug_vcf_fh,$sample,$store_results,$aug_vcf_name)=@_;
	  
	  if ($aug_vcf_fh->{"$sample\_bed"}) {
	      foreach my $key (keys %{$store_results->{"$sample\_bed"}}) {
	   			my $bed_line=$store_results->{"$sample\_bed"}{$key};		
					$aug_vcf_fh->{"$sample\_bed"}->print($bed_line);
				}
	  }
	  my $vcf = Vcf->new(file => $vcf_file);
		$vcf->parse_header();	
		while (my $x = $vcf->next_data_hash()) {
			my $location="$$x{'CHROM'}:$$x{'POS'}:$$x{'REF'}:@{$$x{'ALT'}}[0]";
			my $result_line=$store_results->{$sample}{$location};
			if ($store_results->{$sample}{$location}) { 
				$vcf->add_format_field($x,'MTR');
				$vcf->add_format_field($x,'WTR');
				$vcf->add_format_field($x,'AMB');
		
				$$x{gtypes}{'TUMOUR'}{'MTR'}=$result_line->{'tMTR'};
				$$x{gtypes}{'TUMOUR'}{'WTR'}=$result_line->{'tWTR'};
				$$x{gtypes}{'TUMOUR'}{'AMB'}=$result_line->{'tAMB'};
		
				$$x{gtypes}{'NORMAL'}{'MTR'}=$result_line->{'nMTR'};
				$$x{gtypes}{'NORMAL'}{'WTR'}=$result_line->{'nWTR'};
				$$x{gtypes}{'NORMAL'}{'AMB'}=$result_line->{'nAMB'};
				$aug_vcf_fh->{$sample}->print($vcf->format_line($x));
			}
	}
		$vcf->close();
		$aug_vcf_fh->{$sample}->close();
		my ($aug_gz,$aug_tabix)=_compress_augmented_vcf($aug_vcf_name->{$sample});
                 if ((-e $aug_gz) && (-e $aug_tabix)) {
                  unlink $aug_vcf_name->{$sample} or $log->warn("Could not unlink".$aug_vcf_name->{$sample}.':'.$!);
                }
		return;
}

=head2 _write_output
Write output to file
Inputs
=over 2
=item location -variant position
=item input_bam_files -type array :stores sample names
=item original_vcf_info - original VCF data
=item NFS - New Filter status
=item new_pileup_results --pileup/exonerate results for given position
=item vcf -VCF object
=item tags -custom tags
=item WFH -Write file handler
=back
=cut

sub _write_output {
	my ($input_bam_files,$original_vcf_info, $NFS,$new_pileup_results,$vcf,$tags,$WFH_VCF,$normal,$WFH_TSV,$g_pu,$options)=@_;
  if ((!$vcf && !$options->{'b'})|| $options->{'ao'}) {
		return 1;
	}
  
  if (!defined $original_vcf_info) {$original_vcf_info=['.'];}
  my $out;
	$out->{CHROM}  = $g_pu->{'chr'};
	$out->{POS}    = $g_pu->{'start'};
	$out->{ID}     = '.';
	$out->{ALT}    = [$g_pu->{'alt_seq'}];
	$out->{REF}    = $g_pu->{'ref_seq'};
	$out->{QUAL}   = '.';
	$out->{FILTER} = [$NFS];
	## check if field is in header if not use add_info_field...
	$out->{INFO}   = $original_vcf_info;
	$out->{FORMAT} = [@$tags];
	foreach my $sample (@$input_bam_files) {
	  $out->{gtypes}{$sample} = $new_pileup_results->{$sample};
	}
	
	$vcf->format_genotype_strings($out);
	
	# write to VCF at very end.....
	
	print $WFH_VCF $vcf->format_line($out);
	
	my ($line)=_parse_info_line($vcf->format_line($out),$original_vcf_info,$input_bam_files);
	
	# write to TSV at very end.....
	
	print $WFH_TSV "$normal\t".join("\t",@$line)."\n";
	
return 1;
}

=head2 
checks overlap with high depth regions
Inputs
=over 2
=item chr -chromosome number
=item zero_start start coordiante
=item one_end -end coordinate
=tabix - tabix object of bed file containing UCSC high depth regions
=back
=cut
sub check_hdr_overlap {
	my ($chr, $zero_start, $one_end, $tabix ) = @_;
	$chr=~s/chr//g;
	###
	# Querying is ALWAYS half open regardless of the underlying file type
	###
		my $i=0;
		my $res = $tabix->query("chr$chr", $zero_start, $one_end);
		if(defined $res->get) {
			#if uncomment if want to loop over the hdr regions
			#while(my $record = $tabix->read($res)){
				# the returned data is the ORIGINAL string from the file
				#print "$i : $record\n";
			#}
			$i=1;
			return $i;
		}
	$i;	
}

##### modules from create_config ######

=head2 load_sql_config
load sql statement in connection object
Inputs
=over 2
=item conn -sql connection object
=back
=cut

sub load_sql_config {
  my ($conn) = @_;
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




=head2 build_input_data
parse sql results 
Inputs
=over 2
=item options -user provide input parameters
=item conn -sql connection object
=back
=cut


sub build_input_data {
  my ($options, $conn) = @_;
  my $project_id=$options->{'p'};
  	my $all_data = $conn->executeArrHashRef('nst::NL::getProjectBamAndVcf',$project_id);
  			
		my @retained_data;
		my $total_records = 0;
		for my $curr(@{$all_data}) {
			$total_records++;
			next if(defined $options->{'u'} && (! first { $curr->{'SAMPLE_SYNONYM'} eq $_ } split(',',$options->{'u'})) && ($curr->{'TREAT_AS_TUMOUR'} eq 'Y') );
			next if(defined $options->{'p'} && ! first { $curr->{'ID_INT_PROJECT'} == $_ } $options->{'p'} );
			if(@retained_data > 0) {
				croak "There are multiple species in your requested analysis.\n" if($curr->{'SPECIES'} ne $retained_data[-1]->{'SPECIES'});
				croak "There are multiple builds in your requested analysis.\n" if($curr->{'BUILD'} ne $retained_data[-1]->{'BUILD'});
			}
			push @retained_data, $curr;
		}
		my $retainted_count = @retained_data;
		unless ( $options->{'f'} ) {
			print "Total samples for this project: $retainted_count\n";
			print "\nAre you sure you want to continue? y/n:";
			my $resp = <STDIN>;
			chomp $resp;
			if((lc $resp) ne 'y') {
				$log->warn("exiting...");
				exit(0);
			}
		}
  return \@retained_data;
}
=head2 generate_raw_output
generate raw output using SQL results hash
Inputs
=over 2
=item options -user provide input parameters
=item to_process -results obtained by SQL query
=back
=cut

sub generate_raw_output {
  my ($options, $to_process) = @_;
	my (%sample_group, %normal_samples, $species, $build,$symlinked_files);
  for my $sample_row(@{$to_process}) {
  	if ($sample_row->{'TREAT_AS_TUMOUR'} eq 'N') {
    	$normal_samples{$sample_row->{'ID_INT_PROJECT'}}{$sample_row->{'ID_IND'}}=$sample_row->{'SAMPLE_SYNONYM'}; 
    	$symlinked_files=_process_files($options,$sample_row,$symlinked_files);	
    	next;
    }
    if(!exists $sample_group{$sample_row->{'ID_INT_PROJECT'}}{$sample_row->{'ID_IND'}}) {
    	$sample_group{$sample_row->{'ID_INT_PROJECT'}}{$sample_row->{'ID_IND'}}=$sample_row->{'SAMPLE_SYNONYM'};
		  $species = $sample_row->{'SPECIES'};
			$build = $sample_row->{'BUILD'};
		}
		else {
		  $sample_group{$sample_row->{'ID_INT_PROJECT'}}{$sample_row->{'ID_IND'}}.="\t$sample_row->{'SAMPLE_SYNONYM'}";
		}
   $symlinked_files=_process_files($options,$sample_row,$symlinked_files);
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


sub _process_files {
  my ($options, $input_files,$symlinked_files) = @_;
  my $root_path = $options->{'o'};
  my $out_file = $root_path;
  my ($project, $sample);
  my ($file_types)=_get_file_types();
  my $file_ref = ref $input_files;
  $sample = $input_files->{'SAMPLE_SYNONYM'};
  if((defined $file_ref) && $file_ref eq 'HASH') {
    foreach my $file (keys %$file_types) {
      if(defined $input_files->{$file}) {
      	#if(defined ($input_files->{'CAVE_C'}) and (($file eq 'CAVE') or ($file eq 'CAVE_IDX'))) {next;} 
        my $sym_link=$root_path."/$sample.".$file_types->{$file};
        _create_symlink($input_files->{$file}, $sym_link);
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

sub _get_file_types {
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
                  'PINDEL_BAI' => 'mt.pindel.bam.bai'
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


sub write_config {
  my ($options, $sample_group, $normal_samples, $project, $species, $build, $ref_dir,$ref_dir_x10,$symlinked_files)=@_;
	my ($group_name,$unmatched_normal_count, $data_in_config)='';
	my ($file_types)=_get_file_types();
	#my $data_in_config=0;
	my $single_sample_count=0;
	my $root_path = $options->{'o'};
  my $cfg=Config::IniFiles->new();
  my $config_path="$root_path/$project\_VcfMergeConfig.ini";
  $cfg->SetFileName($config_path);
  $cfg->AddSection($project); 
  foreach my $id_ind(keys %{$sample_group->{$project}}) {
		my @samples_to_analyse=split("\t",$sample_group->{$project}{$id_ind});
		if(@samples_to_analyse >0) {
			
			if(!exists $normal_samples->{$project}{$id_ind} && (defined $options->{'n'} && $options->{'n'} eq 'Y')) {
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
	  	if(defined $options->{'u'}) {
	  		$subset_flag=Sanger::CGP::VcfCompare::VcfMergeAndPileup::check_user_subset(\@samples_to_analyse,$options);
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
	
	if($single_sample_count < 1) {
		$log->warn("Single tumour sample present nothing to merge for: $single_sample_count sample pairs");
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
	if(defined $options->{'b'}) {
		$cfg->AddSection('UserData'); 
		$cfg->newval('UserData','bedfile',$options->{'b'});
	}
	#create symlink for genome build
	$cfg=_symlink_genome_build($species,$build,$cfg,$root_path,$ref_dir,$ref_dir_x10);
  $cfg->RewriteConfig();
  
  $config_path,$root_path;
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

sub _symlink_genome_build {
	my ($species,$build,$cfg,$root_path,$ref_dir,$ref_dir_x10)=@_;
	my $genome_fasta="$root_path/genome.fa";
	my $genome_index="$root_path/genome.fa.fai";
	my $ref_fasta="$ref_dir/$species/$build/genome.fa";
	my $ref_index="$ref_dir/$species/$build/genome.fa.fai";
	if(! -e $ref_fasta) {
		$ref_fasta="$ref_dir_x10/$species/$build/genome.fa";
		$ref_index="$ref_dir_x10/$species/$build/genome.fa.fai";
	}
	
	_create_symlink($ref_fasta, $genome_fasta);
	_create_symlink($ref_index, $genome_index);
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


sub run_vcfmerge {
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

sub _create_symlink {
	my ($file, $symlink_path)=@_;
	if( -l $symlink_path) { $log->warn("symlink exists, skipped symlink creation"); return;} 
	symlink $file, $symlink_path;
}


###### End of common modules ##############

=head2 _parse_info_line
parse info line to generate tab sep outfile
Inputs
=over 2
=item vcf_line vcf --line to be parsed
=item original_vcf_info --VCF info field
=item samples --samples considered for this analysis group
=back
=cut


sub _parse_info_line {
	my ($vcf_line,$original_vcf_info,$samples)=@_;
   # info field prep
		chomp $vcf_line;
    my @data = split "\t", $vcf_line ;

    my (@record,@info_data);

    # basic variant data

    push(@record,$data[2],$data[0],$data[1],$data[3],$data[4],$data[5],$data[6]);

    foreach my $e (split ';',$data[7]){
      push @info_data, [split '=', $e];
    }

    # annotation fields
    my $vdseen = 0;
    foreach my $d (@info_data){
      if($d->[0] eq 'VD'){
        $vdseen = 1;
        my @anno_data = split('\|',$d->[1]);

        foreach my $i(0..4){
          if(defined $anno_data[$i]){
            push(@record,$anno_data[$i]);
          } else {
            push(@record,'-');
          }
        }
      }
    }

    if($vdseen == 0){
      push(@record,'-','-','-','-','-');
    }

    my $vtseen = 0;

		if($vtseen == 0){
      my $ref = $data[3];
      my $alt = $data[4];
      if(length($ref) > 0 && length($alt) > 0){
        if(length($ref) == 1){
          if(length($alt) > 1){
            push(@record,'Ins');
          } else {
            push(@record,'Sub');
          }
        } else {
          if(length($alt) == 1){
            push(@record,'Del');
          } else {
            push(@record,'Complex');
          }
        }
      } else {
        push(@record,'-');
      }

    }

    if($vdseen == 1){
      my $vcseen = 0;
      foreach my $d (@info_data){
        if($d->[0] eq 'VC'){
          $vcseen = 1;
          push(@record,$d->[1]);
        }
      }
      if($vcseen == 0){
        push(@record,'-');
      }
    } else {
      push(@record,'-');
    }

foreach my $info_key (sort keys %$original_vcf_info){
				next if ($info_key eq 'VT' || $info_key eq 'VC');
				if(defined $original_vcf_info->{$info_key}){
         	push(@record,$original_vcf_info->{$info_key});
        }
        else{
        	push(@record,'-');
        }
    }
# print FORMAT field values for each sample

my $i=9; # format field starts from 9
foreach my $sample(@$samples) {
	push(@record,split(':',$data[$i]));
	$i++;
}    

\@record;

}
=head2 _get_tab_sep_header
parse header data to generate header data and columns names
Inputs
=over 2
=item vcf_header
=back
=cut

sub _get_tab_sep_header {
	my ($vcf_header)=shift;
	my $vcf = Vcf->new(file => $vcf_header);
		 $vcf->parse_header();
	my $out;
	my @header;
	push @{$out->{'cols'}}, @$BASIC_COLUMN_TITLES;
	for my $i(0..(scalar @$BASIC_COLUMN_TITLES)-1){
    push @header, '#'.$BASIC_COLUMN_TITLES->[$i] . "\t" . $BASIC_COLUMN_DESCS->[$i]."\n";
  }
	
	my $line_info=$vcf->get_header_line(key=>'INFO');
	foreach my $info_data (@$line_info) {
		foreach my $key (sort keys %$info_data){
			next if ($key eq 'VT' || $key eq 'VC');
			push(@{$out->{'cols'}},$key);
			push @header, '#'.$key. "\t" .$info_data->{$key}{'Description'}."\n";
		}			
	}
	my ($format)=_get_header_lines($vcf->get_header_line(key=>'FORMAT'),'Description','FORMAT');
	push(@header,@$format);
	my ($filter)=_get_header_lines($vcf->get_header_line(key=>'FILTER'),'Description','FILTER');
	push(@header,@$filter);
	my ($org_filter)=_get_header_lines($vcf->get_header_line(key=>'ORIGINAL_FILTER'),'Description','ORIGINAL_FILTER');
	if($org_filter) {
	push(@header,@$org_filter);
	}
	my ($sample)=_get_header_lines($vcf->get_header_line(key=>'SAMPLE'),'SampleName','SAMPLE');
	push(@header,@$sample);
	
$out,\@header;

}

=head2 _get_header_lines
parse header data to generate tab file header and columns names
Inputs
=over 2
=item header_line -- specific header line to parse
=item val -- parse specific value from header
=back
=cut

sub _get_header_lines {
	my($header_line,$val,$prefix)=@_;
	my $header;
	foreach my $header_data (@$header_line) {
		foreach my $key (sort keys %$header_data){
			push @$header, '#'.$prefix.'-:'.$key."\t" .$header_data->{$key}{$val}."\n";
		}	
	}	

return $header;

}

sub _get_process_log {
	my $hash=shift;
	my ($date, $stderr, $exit) = capture {system("date +\"%Y%m%d\"")};
        chomp $date;
        my $log_key=join("_",'vcfProcessLog',$date);
        my $process_log;
	foreach my $key (sort keys %$hash) {
		if(defined($hash->{$key}))
		{
		   chomp $hash->{$key};
                   my $val=$hash->{$key};
		   $process_log->{$key}=_trim_file_path($val);	
		}
	}
return ($log_key,$process_log);
}
# generic sub used for testing...

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
