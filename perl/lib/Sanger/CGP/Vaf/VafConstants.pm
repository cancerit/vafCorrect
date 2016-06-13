package Sanger::CGP::Vaf::VafConstants;

use strict;

use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);

use Sanger::CGP::Vaf; # exports VERSION
use FindBin qw($Bin);

####
#   Base constants
####

const our $LIB_MEAN_INS_SIZE => 'mean_insert_size';
const our $LIB_SD => 'insert_size_sd';
const our $EXECUTE_EXTERNAL => 1;
const our $EXONERATE_SCORE_MULTIPLIER => 5;
const our $EXONERATE_SCORE_FACTOR => 70;
const our $READ_LENGTH_CUTOFF => 2;
const our $READ_LENGTH => 100;
const our $SEP => "/";
const our $test_mode => 0;

const our $INSERT_SIZE_FACTOR => 1;
const our $SPANNING_SEQ_DENOMINATOR => 2;
# alignment constants
const our $PROPER_PAIRED => 0x2;
const our $UNMAPPED => 0x4;
const our $REVERSE_STRAND => 0x10;
const our $MATE_REVERSE_STRAND => 0x20;
const our $NOT_PRIMARY_ALIGN => 0x100;
const our $MATE_UNMAPPED => 0x0008;
const our $READ_PAIRED => 0x0001;
const our $FIRST_IN_PAIR => 0x40;
const our $MAX_PILEUP_DEPTH => '1000000';
const our $SUPP_ALIGNMENT => 0x800;
const our $DUP_READ => 0x400;
const our $VENDER_FAIL => 0x200;
const our @FORMAT_TYPE => qw(MTR WTR AMB UNK VAF);

const our	@BED_HEADER_SNP=> qw(chr pos ref alt FAZ FCZ FGZ FTZ RAZ RCZ RGZ RTZ MTR_NORMAL WTR_NORMAL AMB_NORMAL VAF_NORMAL FAZ FCZ FGZ FTZ RAZ RCZ RGZ RTZ MTR_TUMOUR WTR_TUMOUR AMB_TUMOUR VAF_TUMOUR);
const our @BED_HEADER_INDEL=> qw(chr pos ref alt  MTR_NORMAL WTR_NORMAL AMB_NORMAL UNK_NORMAL VAF_NORMAL MTR_TUMOUR WTR_TUMOUR AMB_TUMOUR UNK_TUMOUR VAF_TUMOUR);

const our $SNP_TAGS => ['FAZ','FCZ','FGZ','FTZ','RAZ','RCZ','RGZ','RTZ','MTR','WTR','DEP','MDR','WDR','VAF','OFS'];
const our $INDEL_TAGS => ['MTR','WTR','DEP','AMB','UNK','MDR','WDR','VAF','OFS'];
const our $VERSION => Sanger::CGP::Vaf->VERSION;

const our $BASIC_COLUMN_TITLES => ['Normal', 'VariantID','Chrom','Pos','Ref','Alt','Qual','Filter','Gene','Transcript','RNA', 'CDS','Protein', 'Type', 'Effect'];
const our $BASIC_COLUMN_DESCS => ['Normal sample name', 'ID of the Variant','Chromosome','Position, the position of a sub or the position immediatly before an indel',
                        'Reference sequence of the Variant','Alternative sequence of the Variant','VCF Quality field','VCF Filter field',
                        'Gene affected, according to the VD tag','Transcript affected, according to the VD tag','Variant description at RNA level',
                        'Variant description at CDS level','Variant description at protein level','Variant type','Summary variants of overall effect'];


1;
