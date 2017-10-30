#!/bin/bash 

##########LICENCE##########
# Copyright (c) 2014-2016 Genome Research Ltd.,
#
# This file is part of cgpVAF.
#
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

SOURCE_VCFTOOLS="http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12a.tar.gz/download"
SOURCE_BIOBDHTS="https://github.com/Ensembl/Bio-HTS/archive/2.3.tar.gz"
SOURCE_HTSLIB="https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2"
SOURCE_SAMTOOLS="https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2"
SOURCE_EXONERATE="http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0.tar.gz"

done_message () {
    if [ $? -eq 0 ]; then
        echo " done."
        if [ "x$1" != "x" ]; then
            echo $1
        fi
    else
        echo " failed.  See setup.log file for error messages." $2
        echo "    Please check INSTALL file for items that should be installed by a package manager"
        exit 1
    fi
}

get_distro () {
  EXT=""
  DECOMP="gunzip -f"
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
    DECOMP="bzip2 -fd"
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi
  rm -f $1.$EXT
  if hash curl 2>/dev/null; then
    curl --retry 10 -sS -o $1.$EXT -L $2
  else
    wget --tries=10 -nv -O $1.$EXT $2
  fi
}

get_file () {
# output, source
  if hash curl 2>/dev/null; then
    curl --insecure -sS -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

if [[ ($# -ne 1 && $# -ne 2) ]] ; then
  echo "Please provide an installation path and optionally perl lib paths to allow, e.g."
  echo "  ./setup.sh /opt/myBundle"
  echo "OR all elements versioned:"
  echo "  ./setup.sh /opt/cgpVafCorrect-X.X.X /opt/PCAP-X.X.X/lib/perl"
  exit 0
fi

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

INST_PATH=$1

# get current directory
INIT_DIR=`pwd`

# create inst_path
mkdir -p $INST_PATH/bin
mkdir -p $INST_PATH/lib
mkdir -p $INST_PATH/config
mkdir -p $INST_PATH/bin/hdr
cp $INIT_DIR/perl/config/log4perl.vaf.conf $INST_PATH/config/
cp -rp $INIT_DIR/README.md	$INST_PATH/README.md

cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
unset PERL5LIB
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"

# log information about this system

    echo '============== System information ===='
    set -x
    lsb_release -a
    uname -a
    sw_vers
    system_profiler
    grep MemTotal /proc/meminfo
    set +x
    echo


perlmods=( "File::ShareDir" "File::ShareDir::Install" "Bio::Root::Version@1.006924" "Module::Build~0.42" )

set -e
for i in "${perlmods[@]}" ; do
  echo -n "Installing build prerequisite $i..."

    set -x
    $INIT_DIR/perl/bin/cpanm -v --mirror http://cpan.metacpan.org -l $INST_PATH $i
    set +x
    echo; echo

  done_message "" "Failed during installation of $i."
done

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

cd $SETUP_DIR

done_message "" "Failed to build $CURR_TOOL."

if [ -e $SETUP_DIR/htslibGet.success ]; then
  echo " already staged ...";
else
  cd $SETUP_DIR
  get_distro "htslib" $SOURCE_HTSLIB
  touch $SETUP_DIR/htslibGet.success
fi

echo -n "Building Bio::DB::HTS ..."
if [ -e $SETUP_DIR/biohts.success ]; then
  echo " previously installed ...";
else
  cd $SETUP_DIR
  rm -rf bioDbHts
  get_distro "bioDbHts" $SOURCE_BIOBDHTS
  echo ls
  mkdir -p bioDbHts/htslib
  tar --strip-components 1 -C bioDbHts -zxf bioDbHts.tar.gz
  tar --strip-components 1 -C bioDbHts/htslib -jxf $SETUP_DIR/htslib.tar.bz2
  cd bioDbHts/htslib
  perl -pi -e 'if($_ =~ m/^CFLAGS/ && $_ !~ m/\-fPIC/i){chomp; s/#.+//; $_ .= " -fPIC -Wno-unused -Wno-unused-result\n"};' Makefile
  make -s -j$CPU
  rm -f libhts.so*
  cd ../
  env HTSLIB_DIR=$SETUP_DIR/bioDbHts/htslib perl Build.PL --install_base=$INST_PATH
  ./Build test
  ./Build install
  cd $SETUP_DIR
  rm -f bioDbHts.tar.gz
  touch $SETUP_DIR/biohts.success
fi

echo -n "Building htslib ..."
if [ -e $SETUP_DIR/htslib.success ]; then
  echo " previously installed ...";
else
  (
  mkdir -p htslib
  tar --strip-components 1 -C htslib -jxf htslib.tar.bz2
  cd htslib
  ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
  make -s -j$CPU
  make install
  cd $SETUP_DIR
  touch $SETUP_DIR/htslib.success
  ) >/dev/null
fi

export HTSLIB=$INST_PATH

if [[ ",$COMPILE," == *,samtools,* ]] ; then
  echo -n "Building samtools ..."
  if [ -e $SETUP_DIR/samtools.success ]; then
    echo " previously installed ...";
  else
    cd $SETUP_DIR
    rm -rf samtools
    get_distro "samtools" $SOURCE_SAMTOOLS
    mkdir -p samtools
    tar --strip-components 1 -C samtools -xjf samtools.tar.bz2
    cd samtools
    ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
    make -s -j$CPU all all-htslib
    make install all all-htslib
    cd $SETUP_DIR
    rm -f samtools.tar.bz2
    touch $SETUP_DIR/samtools.success
  fi
else
  echo "samtools - No change between vafCorrect versions"
fi

cd $SETUP_DIR

CURR_TOOL="vcftools"
CURR_SOURCE=$SOURCE_VCFTOOLS
echo -n "Building $CURR_TOOL ..."
if [ -e $SETUP_DIR/$CURR_TOOL.success ]; then
  echo -n " previously installed ..."
else
    set -ex
    get_distro $CURR_TOOL $CURR_SOURCE
    cd $SETUP_DIR
    mkdir vcftools
    tar --strip 1 -C vcftools -zxf $CURR_TOOL.tar.gz
    cd $CURR_TOOL
    patch Makefile < $INIT_DIR/patches/vcfToolsInstLocs.diff
    patch perl/Vcf.pm < $INIT_DIR/patches/vcfToolsProcessLog.diff
    make -j$CPU PREFIX=$INST_PATH
    touch $SETUP_DIR/$CURR_TOOL.success
fi

done_message "" "Failed to build $CURR_TOOL."

cd $SETUP_DIR

CURR_TOOL="exonerate"
CURR_SOURCE=$SOURCE_EXONERATE
echo -n "Building exonerate..."
  if [ -e $SETUP_DIR/$CURR_TOOL.success ]; then
    echo -n " previously installed ..."
  else
    set -ex
    get_distro $CURR_TOOL $CURR_SOURCE 
    tar zxf exonerate.tar.gz
    cd exonerate-2.2.0
    cp $INIT_DIR/patches/exonerate_pthread-asneeded.diff .
    patch -p1 < exonerate_pthread-asneeded.diff 
    ./configure --prefix=$INST_PATH 
    make -s
    make check 
    make install
    cd $INIT_DIR
    touch $SETUP_DIR/exonerate.success
  fi
  
done_message "" "Failed to build exonerate."


export PATH=$PATH:$INST_PATH/bin
export PERL5LIB=$PERL5LIB:$PERLROOT:$PERLARCH 

cd "$INIT_DIR/perl"

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi

  set -x
  perl $INIT_DIR/perl/bin/cpanm -v --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps . < /dev/null
  set +x

done_message "" "Failed during installation of core dependencies."

echo -n "Installing cgpVaf ..."

  set -e
  cd "$INIT_DIR/perl"
	echo -n `pwd`
  perl Makefile.PL INSTALL_BASE=$INST_PATH
  make
  make test
  make install

done_message "" " cgpVaf install failed."

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of PATH:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo "  $PERLARCH"
echo

exit 0
