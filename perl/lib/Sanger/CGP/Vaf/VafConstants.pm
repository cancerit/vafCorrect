package Sanger::CGP::Vaf::VafConstants;

use strict;

use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);

use Sanger::CGP::Vaf; # exports VERSION
use FindBin qw($Bin);
use Bio::DB::HTS::Constants;

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
const our $DEFAULT_LIB_SIZE => 200;
const our $INSERT_SIZE_FACTOR => 1;
const our $SPANNING_SEQ_DENOMINATOR => 2;
const our $SPANNING_SEQ => 200;
# alignment constants RFLAGS is exported by Bio::DB::HTS::Constants
const our $PROPER_PAIRED => RFLAGS->{'MAP_PAIR'};
const our $UNMAPPED => RFLAGS->{'UNMAPPED'};
const our $REVERSE_STRAND => RFLAGS->{'REVERSED'};
const our $MATE_REVERSE_STRAND => RFLAGS->{'M_REVERSED'};
const our $NOT_PRIMARY_ALIGN => RFLAGS->{'NOT_PRIMARY'};
const our $MATE_UNMAPPED => RFLAGS->{'M_UNMAPPED'};
const our $READ_PAIRED => RFLAGS->{'PAIRED'};
const our $FIRST_IN_PAIR => RFLAGS->{'FIRST_MATE'};
const our $SECOND_IN_PAIR => RFLAGS->{'SECOND_MATE'};
const our $SUPP_ALIGNMENT => RFLAGS->{'SUPPLEMENTARY'};
const our $DUP_READ => RFLAGS->{'DUPLICATE'};
const our $VENDER_FAIL => RFLAGS->{'QC_FAILED'};

const our @FORMAT_TYPE => qw(MTR WTR AMB UNK VAF);
const our $NO_READS_READLEN => 50000;
const our $MAX_PILEUP_DEPTH => '1000000';


const our $DEFAULT_READS_EXCLUDE_FETCH_MATE => RFLAGS->{'NOT_PRIMARY'} +
                                        RFLAGS->{'SUPPLEMENTARY'} +
                                        RFLAGS->{'UNMAPPED'} +
                                        RFLAGS->{'QC_FAILED'};



const our $DEFAULT_FETCH_UNMAPPED =>    RFLAGS->{'NOT_PRIMARY'} +
                                        RFLAGS->{'SUPPLEMENTARY'} +
                                        RFLAGS->{'QC_FAILED'};

const our $DEFAULT_READS_EXCLUDE_PILEUP => RFLAGS->{'NOT_PRIMARY'} +
                                        RFLAGS->{'SUPPLEMENTARY'} +
                                        RFLAGS->{'DUPLICATE'} +
                                        RFLAGS->{'UNMAPPED'} +
                                        RFLAGS->{'QC_FAILED'};
                                        
const our $DEFAULT_READLEN_INCLUDE => RFLAGS->{'MAP_PAIR'} + 
                                        RFLAGS->{'FIRST_MATE'} + 
                                        RFLAGS->{'PAIRED'};

const our @BED_HEADER_SNP=> qw(chr pos ref alt FAZ FCZ FGZ FTZ RAZ RCZ RGZ RTZ MTR_NORMAL WTR_NORMAL AMB_NORMAL VAF_NORMAL FAZ FCZ FGZ FTZ RAZ RCZ RGZ RTZ MTR_TUMOUR WTR_TUMOUR AMB_TUMOUR VAF_TUMOUR);
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
