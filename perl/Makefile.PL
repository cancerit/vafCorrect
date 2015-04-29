#!/usr/bin/perl

use ExtUtils::MakeMaker;
use File::ShareDir::Install;

WriteMakefile(
  NAME          => 'vcfCommon',
  VERSION_FROM  => 'lib/Sanger::CGP::VcfCompare',
  EXE_FILES     => [qw( 
                        bin/mergeAndPileup.pl
                        bin/createConfig.pl
                      )
                   ],
  PREREQ_PM     => {
                     'Module::Build' => 0.42,
                     'IPC::System::Simple' => 1.25,
                     'Capture::Tiny' => 0.24,
                     'Const::Fast' => 0.014,
                     'Bio::DB::Sam' => 1.39,
                     'Try::Tiny' => 0.19,
                     'List::Util' => 1.38,
                     'Test::Fatal' => 0.013,
		     'Math::Round' => 0.06,
                     'Config::IniFiles'=>2.83
);

package MY;
use File::ShareDir::Install qw( postamble );
