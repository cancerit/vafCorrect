#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014-2016
#
# This file is part of cgpVaf
# cgpVaf is free software: you can redistribute it and/or modify it under
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


use ExtUtils::MakeMaker;

WriteMakefile(
  NAME          => 'Vaf',
  LICENSE       => 'agpl_3', # http://search.cpan.org/~dagolden/CPAN-Meta-2.142690/lib/CPAN/Meta/Spec.pm#license
  VERSION_FROM  => 'lib/Sanger/CGP/Vaf.pm',
  EXE_FILES     => [qw( bin/cgpVaf.pl
                        bin/createVafCmd.pl
												)],
  PREREQ_PM     => {
                    'Attribute::Abstract' => '1.02',
										'Bio::DB::HTS' => '1.10',
                    'Bio::Root::Version' => '1.006009001',
                    #'Bio::DB::Sam' => '1.42',
                    'Capture::Tiny' => '0.24',
                    'Config::IniFiles' => '2.83',
                    'Const::Fast' => '0.014',
										'Data::Dumper' => '2.145',
                    'File::Remove' => '0.31',
                    'File::Temp' => '0.2304',
                    'File::Type' => '0.22',
                    'File::ShareDir' => '1.102',
										'File::Path' => '2.09',
										'FindBin' => '1.50',
										'File::Basename' => '2.82',
										'File::Spec' => '3.40',
                    'IPC::System::Simple' => '1.25',
                    'List::Util' => '1.41',
                    'Log::Log4perl' => '1.44',
                    'LWP::Simple' => '6.00',
                    'Pod::Usage' => '1.64',
                    'TAP::Harness' => '3.33',
                    'Try::Tiny' => '0.22',
                   }
);

