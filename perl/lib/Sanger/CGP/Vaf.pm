package Sanger::CGP::Vaf;
use strict;
use Const::Fast qw(const);

our $VERSION = '5.7.0';

const my $LICENSE =>
"#################
# Sanger::CGP::Vaf version %s, Copyright (C) 2014-2020 Wellcome Trust Cancer Genome Project (CGP)
# Sanger::CGP::Vaf comes with ABSOLUTELY NO WARRANTY
# See LICENSE for full details.
#################";

sub license {
  return sprintf $LICENSE, $VERSION;
}

1;
