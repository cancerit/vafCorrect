#Package- merge VCF files and run pileup for SNP and exonerate for indels

package Sanger::CGP::Vaf::Data::ReadVcf; 

BEGIN {
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) .
   ($_[0] =~ m/^Odd number of elements in hash assignment/) || ($_[0] =~m/^Use of uninitialized value \$gtype/) || ($_[0] =~ m/^Use of uninitialized value \$buf/)|| ($_[0] =~ m/symlink exists/) || ($_[0] =~ m/gzip: stdout: Broken pipe/) )};

};

$main::SQL_LIB_LOC = '.'; # this suppresses warnings about uninitialised values
use strict;
use Tabix;
use Vcf;
use Data::Dumper;
use English;
use FindBin qw($Bin);
use List::Util qw(first reduce max min);
use warnings FATAL => 'all';
use Capture::Tiny qw(:all);

use Bio::DB::Sam;
use Bio::DB::Sam::Constants;
use Sanger::CGP::Vaf::VafConstants;
use Sanger::CGP::Vaf::Process::Variant;
use Log::Log4perl;
Log::Log4perl->init("$Bin/../config/log4perl.vaf.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

use base qw(Sanger::CGP::Vaf::Data::AbstractVcf);

1; 


sub _localInit {
	my $self=shift;
	
} 

sub getUniqueLocations {
	my ($self)=@_;
	$self->_checkData($self->getTumourBam,$self->getVcfFile);
	$self->_getData();
	my($data_for_all_samples,$unique_locations,$info_tag_val,$normal_sample,$updated_info_tags)=$self->_mergeVcf();
	
	if(defined $self->{'_b'} ) {
		($unique_locations, $data_for_all_samples,$info_tag_val)=$self->_populate_bed_locations($unique_locations,$data_for_all_samples,$info_tag_val,$updated_info_tags);		
	}	
	if(!defined $self->{'_bo'} ) { 
		print "WARNING!!! more than one normal sample detected for this group \n".print_hash($normal_sample)."\n" if scalar keys %$normal_sample > 1;
		$self->_setNormal($normal_sample);	
	}	

	return($data_for_all_samples,$unique_locations,$info_tag_val,$normal_sample,$updated_info_tags);
}

sub _checkData {
	my ($self,$a,$b)=@_;
	if( ((scalar @$a ) ne (scalar @$b)) && !defined $self->{'_bo'}) {
		$log->logcroak("Unequal number of elemnts in input array, please provide matched number of elements in input array");
	}
}


sub _getData {
	my $self=shift;
	my $count=0;
	foreach my $sample (@{$self->getTumourName}) {
		if( -e ${$self->getTumourBam}[$count]) {
			$self->{'bam'}{$sample}=${$self->getTumourBam}[$count];
		}
		else {
			$log->logcroak("Unble to find TUMOUR BAM file for sample: ".$sample);
		}
		if(defined ${$self->getVcfFile}[$count] &&  (-e ${$self->getVcfFile}[$count]) ) {
			$self->{'vcf'}{$sample}=${$self->getVcfFile}[$count];		
		}
		elsif(!defined $self->{'_bo'}) {
			$log->logcroak("Unble to find VCF file for sample: ".$sample);
		}	
		$count++;	
	}
	if(!-e $self->getNormalBam) {
		$log->logcroak("Unble to find NORMAL BAM file for sample");
	}
	$self->{'bam'}{$self->getNormalName}=$self->getNormalBam;
}

sub _mergeVcf {
	my($self)=@_;
	
	my $vcf_flag=0;
	my $read_depth=0;
	my $tumour_count=0;
	my $info_tag_val=undef;
	my $updated_info_tags=undef;
	my $data_for_all_samples=undef;
	my $vcf_normal_sample=undef;
	my $info_data=undef;
	my $unique_locations=undef;
	
	foreach my $sample (keys %{$self->{'vcf'}}) {
		$tumour_count++;
		if(-e $self->{'vcf'}{$sample}) {	
		
			my $vcf = Vcf->new(file => $self->{'vcf'}{$sample});
			$vcf->parse_header();	
			$vcf->recalc_ac_an(0);
			($info_tag_val,$vcf_normal_sample,$updated_info_tags)=$self->_getOriginalHeader($vcf,$info_tag_val,$vcf_normal_sample,$tumour_count,$sample,$updated_info_tags);
			my $count=0;
			while (my $x=$vcf->next_data_array()) {
				# test
				next if $count > 100;
				$count++;
				##
				if(defined $self->{'_t'}){
					$info_data=$self->_getInfo($vcf,$$x[7],$updated_info_tags);
				}
				#location key consists of CHR:POS:REF:ALT
				my $location_key="$$x[0]:$$x[1]:$$x[3]:$$x[4]";
				if(defined $self->{'_p'}) {
					$read_depth=$self->_getReadDepth($vcf,$x);
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
		}
		else {
		$self->{'noVcf'}{$sample}=0;
		}
 }
	# do bed or bed only analysis
	return ($data_for_all_samples,$unique_locations,$info_tag_val,$updated_info_tags);
}

=head2 _getOriginalHeader
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

sub _getOriginalHeader {
	my ($self,$vcf,$info_tag_val,$normal_sample,$tumour_count,$sample,$updated_info_tags)=@_;
	#stores hash key data in an array for header tags...
	if ($tumour_count == 1){
		my ($vcf_filter,$vcf_info,$vcf_format)=$self->_getCustomHeader();
		#add contig info
		my $contig_info=$vcf->get_header_line(key=>'contig');
		foreach my $contig_line ( @$contig_info) {
			foreach my $chr (sort keys %$contig_line) {
				push(@$info_tag_val,$contig_line->{$chr});
			}
		}
		# get user defined tags
		if(defined $self->{'_t'}) {
			foreach my $tag (split(",",$self->{'_t'})) {
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

return($info_tag_val,$normal_sample,$updated_info_tags);
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

sub _getCustomHeader {
	my ($self)=shift;
	my ($vcf_filter,$vcf_info,$vcf_format);
  my ($log_key,$process_param)=$self->_getProcessLog();
  $vcf_format->{'process_log'}={(key=>$log_key,
			InputVCFSource => _trim_file_path($0),
			InputVCFVer => $Sanger::CGP::Vaf::VafConstants::VERSION,
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
	if ($self->{'_a'} eq 'indel') {
		$vcf_format->{'AMB'}={(key=>'FORMAT',ID=>'AMB', Number=>'1',Type=>'String',Description=>"Reads mapping on both the alleles with same specificity")};
	}
	
	if ($self->{'_a'} eq 'snp') {
		$vcf_format->{'FAZ'}={(key=>'FORMAT',ID=>'FAZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a A for this position, forward strand")};
		$vcf_format->{'FCZ'}={(key=>'FORMAT',ID=>'FCZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a C for this position, forward strand")};
		$vcf_format->{'FGZ'}={(key=>'FORMAT',ID=>'FGZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a G for this position, forward strand")};
		$vcf_format->{'FTZ'}={(key=>'FORMAT',ID=>'FTZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a T for this position, forward strand")};
		$vcf_format->{'RAZ'}={(key=>'FORMAT',ID=>'RAZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a A for this position, reverse strand")};
		$vcf_format->{'RCZ'}={(key=>'FORMAT',ID=>'RCZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a C for this position, reverse strand")};
		$vcf_format->{'RGZ'}={(key=>'FORMAT',ID=>'RGZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a G for this position, reverse strand")};
		$vcf_format->{'RTZ'}={(key=>'FORMAT',ID=>'RTZ', Number=>'1',Type=>'Integer',Description=>"Reads presenting a T for this position, reverse strand")};	
	}
	
	return ($vcf_filter,$vcf_info,$vcf_format);
}


sub _getProcessLog {
	my ($self)=@_;
	my $hash=$self->{'options'};
	my ($date, $stderr, $exit) = capture {system("date +\"%Y%m%d\"")};
        chomp $date;
        my $log_key=join("_",'vcfProcessLog',$date);
        my $process_log;
	foreach my $key (sort keys %$hash) {
		if(defined $hash->{$key})
		{   
		   chomp $hash->{$key};
       my $val=$hash->{$key};
		   $process_log->{$key}=$self->_trim_file_path($val);	
		}
	}
return ($log_key,$process_log);
}


sub _trim_file_path{
	my ($self,$string ) = @_;
	
	return 'NA' unless (defined $string);
	
	my @bits = (split("/", $string));
	return pop @bits;
}

=head2 _getInfo
parse VCF INFO line for given locations and returns values for respective tags
Inputs 
=over 2
=item vcf - vcf object
=item INFO - VCF INFO filed value
=item INFO - INFO field tags defined in header 
=back
=cut
sub _getInfo {
	my ($self,$vcf,$INFO,$info_tags)=@_;
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


=head2 _getReadDepth
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

sub	_getReadDepth {
	my ($self,$vcf,$line)=@_;
	my @tags=split(',',$self->{'_p'});
	my $depth=0;
	for my $tag (@tags) {
		my $idx = $vcf->get_tag_index($$line[8],$tag,':'); 
		my $pl  = $vcf->get_field($$line[10],$idx) unless $idx==-1;
		$depth+=$pl;
	}
 return $depth;
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
	my ($self,$unique_locations,$data_for_all_samples,$info_tag_val,$updated_info_tags)=@_;
	my $bed_file=$self->{'_b'};
	my $bed_name=$self->_trim_file_path($bed_file);
	open(my $bedFH, $self->{'_b'})|| $log->logcroak("unable to open file $!");
	my $location_counter=0;
	my %info_tag;
	
	$info_tag{'Interval'}='BedFile';
	if( !defined $self->{'_vcf'} && !defined $self->{'_ao'} ) {
		my ($vcf_filter,$vcf_info,$vcf_format)=$self->_getCustomHeader();
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
		
		$unique_locations->{$location_key}="$bed_name-BEDFILE";
		#create custom tag hash
		my $temp_tag_val;
		
		foreach my $tag (@$updated_info_tags) { 
			$temp_tag_val->{$tag}='-';
		}
		#add location key to all samples
		foreach my $sample(keys %{$self->{'vcf'}})	{
			$data_for_all_samples->{$sample}{$location_key}={'INFO'=>$temp_tag_val,'FILTER'=>'NA','RD'=>1 };
		}
		$location_counter++;
	}
	$log->debug(" Added additional ( $location_counter ) locations from bed file:$bed_file");
	return ($unique_locations,$data_for_all_samples,$info_tag_val);
}

sub _setNormal {
	my($self,$vcf_normal)=@_;
	foreach my $key (keys %$vcf_normal) {
		if($key eq $self->getNormalName) {
			$log->debug(" User provided normal sample matches with VCF normal");
		}
		elsif(defined $self->{'_vn'}){
			$self->setNormal($key);
			$log->debug("Setting normal sample as specified in VCF header : make sure normal sample ($key) bam file is present");
		}
	}
}


sub processLocations {
	my($self,$data_for_all_samples,$unique_locations,$info_tag_val,$updated_info_tags)=@_;
		
	my($bam_objects,$bas_files)=$self->_get_bam_object();
	my($bam_header_data,$lib_size)=$self->_get_bam_header_data($bam_objects,$bas_files);
	my($aug_vcf_fh,$aug_vcf_name)=$self->_write_augmented_header();
	my($vcf,$WFH_VCF,$WFH_TSV)=$self->_write_vcf_header($info_tag_val);
	
	my $variant=Sanger::CGP::Vaf::Process::Variant->new( 
		'location' 		=> undef,
		'varLine' 		=> undef,
		'varType' 		=> $self->{'_a'},
		'libsize' 		=> $lib_size,
		'samples' 		=> $self->getAllSampleNames,
		'tumourName'	=> $self->getTumourName,
		'normalName'	=> $self->getNormalName,
		'vcfStatus' 	=> $self->{'_vcf'},
		'noVcf'    		=> $self->{'noVcf'},
		'outDir'			=> $self->getOutputDir,
		'passedOnly'  => $self->{'_r'},
		'tabix_hdr' 		=> new Tabix(-data => "$Bin/hdr/seq.cov".$self->{'_c'}.'.ONHG19_sorted.bed.gz')
		);
		
	my $store_results=undef;
	
	foreach my $location (sort keys %$unique_locations) {
		$variant->setLocation($location);
		$variant->setVarLine($unique_locations->{$location});
		my($g_pu)=$variant->formatVarinat();
		#process only passed varinats		
		if($self->{'_r'} && $variant->getVarLine!~/PASS/ && $variant->getVarLine!~/BEDFILE/) {
			if(defined $self->{'_ao'} || defined $self->{'_m'}) {
				foreach my $sample (@{$self->getTumourName}) {
					$store_results=$variant->storeResults($store_results,$g_pu,$sample);
				}
				next;
			}
		}
		
    my ($old_info_val,$NFS,$original_flag,$max_depth)=$variant->getVcfFields($data_for_all_samples);
   
    if ($self->{'_a'} eq 'indel') {
			$g_pu=$variant->createExonerateInput($bam_objects->{$self->getNormalName},$bam_header_data,$max_depth,$g_pu);
    }
   
  	foreach my $sample (@{$self->getAllSampleNames}) {
   		print "------$sample\n";
			#$sample_counter++;
			$g_pu=$variant->populateHash($g_pu,$sample,$bam_header_data); # reset the counter for counts and lib size;
			if($self->{'_a'} eq 'indel') {	
				$g_pu=$variant->getIndelResults($bam_objects->{$sample},$g_pu);
				if( ($sample eq $self->getNormalName) && (defined $self->{'_m'}) ) {
					$g_pu=$variant->addNormalCount($g_pu);
				}
				elsif($self->{'_m'} && $data_for_all_samples->{$sample}{$location}) {
					$store_results=$variant->_storeResults($store_results,$g_pu,$sample);
				}
			}
			else{
				$g_pu=$variant->getPileup($bam_objects->{$sample},$g_pu);
				print Dumper $g_pu;
			}
			if(!defined $options->{'ao'} ) {
				$variant->formatResults($original_flag,$g_pu);
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
		
			
			
			
  	}# sample	...	
   
   exit;
   
    
	}# location ...
	
		
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
	my ($self)=@_;
	my (%bam_objects,%bas_files);
	my $files=$self->{'bam'};
	foreach my $sample (keys %$files) {
		my $sam = Bio::DB::Sam->new(-bam => $files->{$sample},
															-fasta =>$self->getGenome,
															-expand_flags => 1);
		$sam->max_pileup_cnt($Sanger::CGP::Vaf::VafConstants::MAX_PILEUP_DEPTH);
		$bam_objects{$sample}=$sam;
		$bas_files{$sample}=$files->{$sample}.'.bas';
	}
	return(\%bam_objects,\%bas_files);
}


=head2 _get_bam_header_data
get_bam_header_data -insert size and chr length
Inputs
=over 2
=item bam_objects - Bio::DB sam object
=back
=cut

sub _get_bam_header_data {
	my ($self,$bam_objects,$bas_files)=@_;
  my ($chr_len,%bam_header_data);
  
  return if $self->{'_a'} ne 'indel';
  
  my $lib_size=0;
  foreach my $key (keys %$bam_objects) {	
  	my($mapped_length)=$self->_get_read_length($bam_objects->{$key});
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
  		 $lib_size=$self->_get_lib_size_from_bas($bas_files->{$key});
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
	
	return(\%bam_header_data,$lib_size);
}

=head2 _get_lib_size_from_bas
get library size from BAS file
Inputs
=over 2
=item bas file name
=back
=cut

sub _get_lib_size_from_bas {
  my ($self,$bas_file)=@_;
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
							$index_mi = first { $array[$_] eq $Sanger::CGP::Vaf::VafConstants::LIB_MEAN_INS_SIZE } 0 .. $#array;
							$index_sd = first { $array[$_] eq $Sanger::CGP::Vaf::VafConstants::LIB_SD } 0 .. $#array;
					}
	}
	close($bas);
	return $lib_size;
}

=head2 _get_read_length
get read length per sample
Inputs
=over 2
=item sam - Bio::DB::Sam objects
=back
=cut

sub _get_read_length {
	my ($self,$sam)=@_;
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
	my($self)=@_;
	my $aug_vcf_fh=undef;
	my $aug_vcf_name=undef;
	
	return if(!defined $self->{'_vcf'}); 
	
	if(defined $self->{'_m'}) {
		my $augment_vcf=$self->{'vcf'};
		my ($vcf_filter,$vcf_info,$vcf_format)=$self->_getCustomHeader();
		foreach my $sample (keys %$augment_vcf) { 
			my ($tmp_file)=(split '/', $augment_vcf->{$sample})[-1];
			$tmp_file=~s/(\.vcf|\.gz)//g;
			my $aug_vcf=$self->getOutputDir.'/'.$tmp_file.$self->{'_oe'};
			$aug_vcf_name->{$sample}=$aug_vcf;
			open(my $tmp_vcf,'>',$aug_vcf)|| $log->logcroak("unable to open file $!");
			$log->debug("Augmenting vcf file:".$self->getOutputDir."/$tmp_file.".$self->{'_oe'});
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
    
	if(defined $self->{'_b'} && defined $self->{'_m'}) {
		my $input_bam_files=$self->{'bam'};
		my @bed_header=qw(chr pos ref alt FAZ FCZ FGZ FTZ RAZ RCZ RGZ RTZ MTR_NORMAL WTR_NORMAL AMB_NORMAL FAZ FCZ FGZ FTZ RAZ RCZ RGZ RTZ MTR_TUMOUR WTR_TUMOUR AMB_TUMOUR);
		if($self->{'_a'} eq 'indel') {
			@bed_header=qw(chr pos ref alt  MTR_NORMAL WTR_NORMAL AMB_NORMAL MTR_TUMOUR WTR_TUMOUR AMB_TUMOUR);
		}
		foreach my $sample (keys %$input_bam_files) {
			if ($sample ne $self->getNormalName) {
				open(my $tmp_bed,'>',$self->getOutputDir."/$sample.augmented.bed");
				$log->debug("Augmented bed file:".$self->getOutputDir."/$sample.augmented.bed");
				print $tmp_bed join("\t",@bed_header)."\n";
				$aug_vcf_fh->{"$sample\_bed"}=$tmp_bed;	
			}
		}
	}
  return($aug_vcf_fh,$aug_vcf_name);
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
	my($self,$info_tag_val)=@_;
	
	my $vcf=undef;
	my $WFH_VCF=undef;
	my $WFH_TSV=undef;
	
	# return if no VCF file found or augment only option is provided for indel data ...
	if(( !defined $self->{'vcf'} && !defined $self->{'_b'}) || (defined $self->{'_ao'})) {
			return $vcf,$WFH_VCF,$WFH_TSV;
	}
	
	my $outfile_name=$self->getOutputDir.'/'.$self->getNormalName.'_'.@{$self->getTumourName}[0]."_consolidated_".$self->{'_a'};	
	$log->debug("VCF outfile:$outfile_name.vcf");
	$log->debug("TSV outfile:$outfile_name.tsv");
	open($WFH_VCF, '>',"$outfile_name.vcf");
	open($WFH_TSV, '>',"$outfile_name.tsv");

	$vcf=Vcf->new();
	my $genome_name=$self->_trim_file_path($self->getGenome);
	my $script_name=$self->_trim_file_path($0);
	$vcf->add_header_line({key=>'reference', value=>$genome_name}); 
	$vcf->add_header_line({key=>'source', value=>$script_name}); 
	$vcf->add_header_line({key=>'script_version', value=>$Sanger::CGP::Vaf::VafConstants::VERSION});
	$vcf->add_header_line({key=>'Date', value=>scalar(localtime)});
	if(defined $info_tag_val) {
		foreach my $hash_val(@$info_tag_val) {
			$vcf->add_header_line($hash_val);	
		}
	}
	else {
		$log->debug("No data in info filed");
	}
	
	if(!defined $self->{'_vcf'}) {
			my $i=0;
			foreach my $sample (@{$self->getAllSampleNames}) {
				if ($sample eq $self->getNormalName) {
					my %temp=(key=>"SAMPLE",ID=>"NO_VCF_NORMAL", Description=>"NO_VCF_DATA", SampleName=>$sample);
					$vcf->add_header_line(\%temp);
				}
				else {
					$i++;
					my %temp=(key=>"SAMPLE",ID=>"NO_VCF_TUMOUR_$i", Description=>"NO_VCF_DATA_$i", SampleName=>$sample);
					$vcf->add_header_line(\%temp);
				}
		}
	}
	
	$vcf->add_columns(@{$self->getAllSampleNames});
	print $WFH_VCF $vcf->format_header();	

  # writing results in tab separated format
	my $tmp_file=$self->{'_o'}.'/temp.vcf';
	open(my $tmp_vcf,'>',$tmp_file);
	print $tmp_vcf $vcf->format_header();
	close($tmp_vcf);
	my($col_names,$header,$format_col)=$self->_get_tab_sep_header($tmp_file);
	my $temp_cols=$col_names->{'cols'};
	
	my $tags=$Sanger::CGP::Vaf::VafConstants::SNP_TAGS;
	if($self->{'_a'} eq 'indel') {
		$tags=$Sanger::CGP::Vaf::VafConstants::INDEL_TAGS
	}
	
	foreach my $sample(@{$self->getAllSampleNames}){
		foreach my $tag_name(@$tags){
			push ($temp_cols,"$sample\_$tag_name");
		}
	}   
	print $WFH_TSV "@$header\n".join("\t",@$temp_cols)."\n";
			
	return ($vcf,$WFH_VCF,$WFH_TSV);
}


=head2 _get_tab_sep_header
parse header data to generate header data and columns names
Inputs
=over 2
=item vcf_header
=back
=cut

sub _get_tab_sep_header {
	my ($self,$vcf_header)=@_;
	my $vcf = Vcf->new(file => $vcf_header);
	$vcf->parse_header();
	my $out;
	my @header;
	push @{$out->{'cols'}}, @$Sanger::CGP::Vaf::VafConstants::BASIC_COLUMN_TITLES;
	for my $i(0..(scalar @$Sanger::CGP::Vaf::VafConstants::BASIC_COLUMN_TITLES) -1){
    push @header, '#'.$Sanger::CGP::Vaf::VafConstants::BASIC_COLUMN_TITLES->[$i] . "\t" . $Sanger::CGP::Vaf::VafConstants::BASIC_COLUMN_DESCS->[$i]."\n";
  }
	my $line_info=$vcf->get_header_line(key=>'INFO');
	foreach my $info_data (@$line_info) {
		foreach my $key (sort keys %$info_data){
			next if ($key eq 'VT' || $key eq 'VC');
			push(@{$out->{'cols'}},$key);
			push @header, '#'.$key. "\t" .$info_data->{$key}{'Description'}."\n";
		}			
	}
	my ($format)=$self->_get_header_lines($vcf->get_header_line(key=>'FORMAT'),'Description','FORMAT');
	push(@header,@$format);
	my ($filter)=$self->_get_header_lines($vcf->get_header_line(key=>'FILTER'),'Description','FILTER');
	push(@header,@$filter);
	my ($org_filter)=$self->_get_header_lines($vcf->get_header_line(key=>'ORIGINAL_FILTER'),'Description','ORIGINAL_FILTER');
	if($org_filter) {
	push(@header,@$org_filter);
	}
	my ($sample)=$self->_get_header_lines($vcf->get_header_line(key=>'SAMPLE'),'SampleName','SAMPLE');
	push(@header,@$sample);
	
return ($out,\@header);

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
	my($self,$header_line,$val,$prefix)=@_;
	my $header;
	foreach my $header_data (@$header_line) {
		foreach my $key (sort keys %$header_data){
			push @$header, '#'.$prefix.'-:'.$key."\t" .$header_data->{$key}{$val}."\n";
		}	
	}	

return $header;

}









