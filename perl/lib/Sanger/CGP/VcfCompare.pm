package Sanger::CGP::VcfCompare;
use strict;
use Const::Fast qw(const);

our $VERSION = '3.1.0';

const my $LICENSE =>
"#################
# Sanger::CGP::VcfCompare version %s, Copyright (C) 2014 Wellcome Trust Cancer Genome Project (CGP)
# Sanger::CGP::VcfCompare comes with ABSOLUTELY NO WARRANTY
# See LICENSE for full details.
#################";

sub license {
  return sprintf $LICENSE, $VERSION;
}

1;

