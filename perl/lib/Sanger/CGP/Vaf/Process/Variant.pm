#Package- merge VCF files and run pileup for SNP and exonerate for indels

package Sanger::CGP::Vaf::Process::Variant; 

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

use Log::Log4perl;
Log::Log4perl->init("$Bin/../config/log4perl.vaf.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

use base qw(Sanger::CGP::Vaf::Process::AbstractVariant);

1; 


sub _localInit {
	my $self=shift;
	
} 

=head2 formatVarinat
get formatted position, for pindel - flag values for insertion del
Inputs
=over 2
=item location -variant position
=item variant_type - variant type [indel or SNP ]

=back
=cut

sub formatVarinat {
  my ($self) = @_;
  my $insertion_flag=undef;
  my $del_flag=undef;
  my $g_pu=undef; 
  my ($chr,$start,$ref,$alt)=(split ':', $self->getLocation)[0,1,2,3];
  my $end=$start;
  #different for indels...
  if($self->getVarType eq 'indel') {
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
		#case of --> complex indel
		if(defined $ref && length($alt) > 1) { $insertion_flag=1;}
  }

	$g_pu={'chr'=>$chr, 
				'start' => $start, 
				'end' => $end, 
				'ref_seq' => $ref,
				'alt_seq' => $alt,
				'region' => $chr.':'.$start.'-'.$end, 
				'ins_flag' => $insertion_flag,
				'del_flag' => $del_flag, 
			 };

return $g_pu;

}



sub storeResults {
	my ($self,$store_results,$g_pu,$sample)=@_;
	my $results = {'tMTR'=> '.', 
							   'tWTR'=> '.',
							 	 'tAMB'=> '.',
							 	 'nMTR'=> '.',
							 	 'nWTR'=> '.',
							 	 'nAMB'=> '.',								  
		};	
	if (!exists $g_pu->{'alt_p'}) {
		$store_results->{$sample}{$self->getLocation}=$results;
		return $store_results;
	}
	
  my $MTR = $g_pu->{'alt_p'} + $g_pu->{'alt_n'};
	my $WTR = $g_pu->{'ref_p'} + $g_pu->{'ref_n'};

	if ($self->getVarLine=~/BEDFILE/) {
	    my $bed_line=$g_pu->{'chr'}."\t".
	    						 $g_pu->{'start'}."\t".
	    						 $g_pu->{'ref_seq'}."\t".
	    						 $g_pu->{'alt_seq'}."\t".
									 $g_pu->{'normal_MTR'}."\t".
									 $g_pu->{'normal_WTR'}."\t".
									 $g_pu->{'normal_AMB'}."\t".
									 $MTR."\t".
									 $WTR."\t".
									 $g_pu->{'amb'}."\n";
									 
			$store_results->{"$sample\_bed"}{$self->getLocation}=$bed_line;		
			return $store_results;
		}
	$results->{'tMTR'}=$MTR;
	$results->{'tWTR'}=$WTR;
	$results->{'tAMB'}=$g_pu->{'amb'};
	
	$results->{'nMTR'}=$g_pu->{'normal_MTR'};
	$results->{'nWTR'}=$g_pu->{'normal_WTR'};
	$results->{'nAMB'}=$g_pu->{'normal_AMB'};
	
  $store_results->{$sample}{$self->getLocation}=$results;
  
  return $store_results; 
 
}


=head2 getVcfFields
get sample specific original values for INFO and FILTER fields
Inputs
=over 2
=item sample_names -sample names
=item data_for_all_samples -original VCF fileds for given sample
=item key_location -variant position
=item normal -normal sample name

=back
=cut

sub getVcfFields {
	my($self,$data_for_all_samples)=@_;
	
	return if(!defined $self->getVcfStatus && $self->getVarLine!~/BEDFILE/);
	my $max_depth=0;
	my $flag_val=undef;
	my $original_vcf_info=undef;
	
  foreach my $sample (@{$self->{'_tumourName'}}) {
		if(exists $data_for_all_samples->{$sample}{$self->getLocation} ) {
			my $info_line=$data_for_all_samples->{$sample}{$self->getLocation}->{'INFO'};
			my $filter_val=$data_for_all_samples->{$sample}{$self->getLocation}->{'FILTER'};
			my $max_rd=$data_for_all_samples->{$sample}{$self->getLocation}->{'RD'};
			$max_depth=$max_rd if $max_rd > $max_depth;
			$flag_val->{$sample}=$filter_val;
			$original_vcf_info->{$sample}=$info_line;
		}
		else {
			$original_vcf_info->{$sample}='NA';
			$flag_val->{$sample}='NA';
		}
	}
  
  my ($NFS,$DNFS)=$self->_getFlagStatus($flag_val);
	my ($old_info_val)=$self->_getINFO($original_vcf_info,$DNFS);
	
	return ($old_info_val,$NFS,$flag_val,$max_depth);
}



=head2 _getFlagStatus
get new filter status based on original filter status
Inputs
=over 2
=item flag_val -Original flag values in filter column
=item normal -normal sample name

=back
=cut


sub _getFlagStatus {
	my ($self,$flag_val)=@_;
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
	if ($self->getVarLine=~m/BEDFILE/) {
		return ('BD',"$samples:$called_n_passed:$passed");
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
=item DNFS -Descriptive flag status [ref header info]
=back
=cut

sub _getINFO {
	my ($self,$original_vcf_info,$DNFS)=@_;
	my $old_info_val;
	foreach my $sample (keys %$original_vcf_info) {
		next if !defined $original_vcf_info->{$sample} || $original_vcf_info->{$sample} eq 'NA';
		my $info=$original_vcf_info->{$sample};
		foreach my $tag (keys %$info) {
			$old_info_val->{$tag}=$info->{$tag};
		}
	}
	$old_info_val->{'NS'}=((split ':', $DNFS)[0]) + scalar keys %{$self->{'_novcf'}};
	$old_info_val->{'NC'}=(split ':', $DNFS)[1];
	$old_info_val->{'NP'}=(split ':', $DNFS)[2];
	$old_info_val->{'NA'}=(split ':', $DNFS)[0];

	return $old_info_val;
}	


=head2 getRange
get left and right span from  indel position
Inputs
=over 2
=item bam_header_data - sample specific information on read length, chr length and lib size 
=item sample -sample name
=item g_pu - has containing sample specific info
=item tabix_hdr - tabix object created from UCSC high depth region bed file
=back
=cut

sub getRange {
  my($self,$bam_header,$g_pu,$max_depth)=@_;
  
  return if $self->getVarType ne 'indel';
  
  my ($left_pos,$right_pos,$chr_len,$spanned_region);
  
  my $lib_size=$bam_header->{$self->{'_normalSample'}}{'lib_size'};
  
  my ($hdr_flag)=$self->check_hdr_overlap($g_pu->{'chr'},$g_pu->{'start'},$g_pu->{'end'},$self->{'_tabix_hdr'});
  my $spanning_seq_denom=$Sanger::CGP::Vaf::VafConstants::SPANNING_SEQ_DENOMINATOR;
  #if location is in high depth region and has depth >1000 then hdr_flag is true
  if($hdr_flag && $max_depth > 1000){$spanning_seq_denom=4;}
  else {$hdr_flag=0;}
  $chr_len=$bam_header->{$self->{'_normalSample'}}{$g_pu->{'chr'}};
  if(defined $lib_size && defined $chr_len) {
  	$spanned_region = round(($lib_size *  $Sanger::CGP::Vaf::VafConstants::INSERT_SIZE_FACTOR )/$spanning_seq_denom);
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
	# to get exact location of variant base and avoid borderline matching   
	# added padding to reference posotion 
	$g_pu->{'ref_pos_5p'}=($g_pu->{'start'} - $g_pu->{'pos_5p'}) + 1 ;
	if($g_pu->{'ins_flag'} && !$g_pu->{'del_flag'}){
		$g_pu->{'ref_pos_3p'}=$g_pu->{'ref_pos_5p'} + ( $g_pu->{'end'} - $g_pu->{'start'} );
	}else{
		$g_pu->{'ref_pos_3p'}=$g_pu->{'ref_pos_5p'} + ( $g_pu->{'end'} - $g_pu->{'start'} ) - 1;
	}
	
	$g_pu->{'alt_pos_3p'}=$g_pu->{'ref_pos_5p'} + length( $g_pu->{'alt_seq'}) -1;

	return $g_pu;
}

























