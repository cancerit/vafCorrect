package Sanger::CGP::Vaf::Process::AbstractVariant;
																					 
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


use strict;

use Log::Log4perl;
use POSIX qw(ceil);
use Data::Dumper;
use Attribute::Abstract;

use Sanger::CGP::Vaf::VafConstants;
my $log = Log::Log4perl->get_logger(__PACKAGE__);

1;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
	my $self = {};
	bless($self, $class);
	$self->_init(@_);
	$self->_localInit(@_);
	return $self;
}


=head2 _init
get the user input
Inputs
=over 2
=item options - user provided options to get file extension, input dir path , tag values and bed file
=back
=cut
 
sub _init {
	my $self = shift;
	my (%vars) = @_;
	foreach my $key (keys %vars) {
		$self->{"_$key"}=$vars{$key};
	}
}

sub _localInit: Abstract;

sub getLocation {
	return shift->{'_location'};
}

sub setLocation {
	my($self,$loc)=@_;
	$self->{'_location'}=$loc;
}

sub getVarLine{
	return shift->{'_varLine'};
}

sub setVarLine {
	my($self,$var)=@_;
	$self->{'_varLine'}=$var;
}

sub getVarType {
	return shift->{'_varType'};
}

sub getLibSize {
	return shift->{'_libsize'};
}

sub getVcfStatus {
	return shift->{'_vcfStatus'};
}

sub getBamObj {
	return shift->{'_bamObj'};
}

sub getSample {
	return shift->{'_samples'};
}
#-----Legacy
sub addMessage {
	my ($self,$msg) = @_;
	push(@{$self->{_msg}},ref($self).": ".$msg);
}

sub _debug {
	my $self = shift;
	if(exists($self->{_debug}) && defined($self->{_debug}) && $self->{_debug}){
		return 1;
	} else {
		return 0;
	}
}

sub getMessages {
	my $self = shift;
	return @{$self->{_msg}} if defined($self->{_msg});
	return undef;
}

sub _clearMessages {
	shift->{_msg} = undef;
}

__END__

=head1 NAME

Sanger::CGP::Vaf::Process::AbstractVariant - Abstract base class for the variant allele fraction analysis 

=head1 DESCRIPTION

This is an abstract template class for the VAF, it provides a lot of shared behind the scenes functionality.  All

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $source = Sanger::CGP::Vaf::Process::AbstractVariant->new(%params);

=item Function :

Builds a new Sanger::CGP::Vaf::Process::AbstractVariant inheriting object

=item Returns :

Sanger::CGP::Vaf::Process::AbstractVariant object initialized with parameter values

=item Params :

Hash of parameter values
user input parameters

=back
