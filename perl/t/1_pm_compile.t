##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of VAGrENT.
#
# VAGrENT is free software: you can redistribute it and/or modify it under
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
##########LICENCE##########

# this is a catch all to ensure all modules do compile

use strict;
use Data::Dumper;
use Test::More;
use List::Util qw(first);
use File::Find;
use Cwd;
use Try::Tiny qw(try finally);
use Capture::Tiny qw(capture);
use File::Spec;
use Const::Fast qw(const);

use FindBin qw($Bin);
my $lib_path = "$Bin/../lib";

# Add modules here that cannot be instantiated (should be extended and have no 'new')
# or need a set of inputs - these should be tested in own test script
const my @USE_SKIP => qw(  );
const my @NEW_SKIP => qw( Sanger::CGP::Vaf::Data::AbstractVcf
                          Sanger::CGP::Vaf::Data::ReadVcf
                          Sanger::CGP::Vaf::Process::AbstractVariant
                          Sanger::CGP::Vaf::Process::Variant
                         );

my $init_cwd = getcwd;

my $bail_count;

my @modules;
try {
  chdir($lib_path);
  find({ wanted => \&compile_modules, no_chdir => 1 }, './');
} finally {
  chdir $init_cwd;
  die "The try block died with: @_\n" if(@_);
};

if($bail_count) {
  BAIL_OUT("Modules failed to compile, see above errors");
}

for my $mod(@modules) {
  use_ok($mod);
  ok($mod->VERSION, "Check version inheritance exists ($mod)");
  if($mod->can('new')) { # only try new on things that have new defined
    new_ok($mod) unless( first {$mod eq $_} (@USE_SKIP, @NEW_SKIP) );
  }
}

done_testing();

sub compile_modules {
#return unless($_ =~ m/SequenceOntologyClassifier.pm/);
  if($_ =~ m/\.pm$/) {
    my ($stdout, $stderr, $exit) = capture { system("$^X -c $_"); };
    warn "\n$stderr" if($exit);
    is($exit, 0, "Module fails to compile cleanly: $_");

    $bail_count += $exit;
  }
}

sub build_module_set {
  if($_ =~ m/\.pm$/) {
    my ($dir_str,$file) = (File::Spec->splitpath( $_ ))[1,2];
    $file =~ s/\.pm$//;
    my @dirs = File::Spec->splitdir( $dir_str );
    shift @dirs;
    push @modules, (join '::', @dirs).$file;
  }
}
