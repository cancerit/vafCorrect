#!/bin/bash

set -ex

if [[ -z "${TMPDIR}" ]]; then
  TMPDIR=/tmp
fi

set -u

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path such as /opt/ICGC"
  exit 1
fi

# get path to this script
SCRIPT_PATH=`dirname $0`;
SCRIPT_PATH=`(cd $SCRIPT_PATH && pwd)`

# get the location to install to
INST_PATH=$1
mkdir -p $1
INST_PATH=`(cd $1 && pwd)`
echo $INST_PATH

# get current directory
INIT_DIR=`pwd`

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

mkdir -p $INST_PATH/bin

# make sure tools installed can see the install loc of libraries
set +u
export LD_LIBRARY_PATH=`echo $INST_PATH/lib:$LD_LIBRARY_PATH | perl -pe 's/:\$//;'`
export PATH=`echo $INST_PATH/bin:$PATH | perl -pe 's/:\$//;'`
export MANPATH=`echo $INST_PATH/man:$INST_PATH/share/man:$MANPATH | perl -pe 's/:\$//;'`
set -u

cd perl

mkdir -p cache_tmp
curl -sSL https://raw.githubusercontent.com/samtools/samtools/1.9/misc/seq_cache_populate.pl > cache_tmp/seq_cache_populate.pl

mkdir -p CACHE
zcat t/testData/genome.fa.gz > cache_tmp/genome.fa
perl ./cache_tmp/seq_cache_populate.pl -root ./CACHE cache_tmp/genome.fa

export REF_CACHE=$PWD/CACHE/%2s/%2s/%s
export REF_PATH=$REF_CACHE

cpanm --no-wget --no-interactive --notest --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
cpanm -v --no-wget --no-interactive --mirror http://cpan.metacpan.org -l $INST_PATH .

rm -rf CACHE cache_tmp

mkdir -p $INST_PATH/config
cp config/log4perl.vaf.conf $INST_PATH/config
chmod +r $INST_PATH/config/log4perl.vaf.conf

cd $HOME
