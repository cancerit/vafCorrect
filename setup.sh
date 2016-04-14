

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

SOURCE_SAMTOOLS="https://github.com/samtools/samtools/archive/0.1.20.tar.gz"
SOURCE_TABIX="https://github.com/sb43/tabix/archive/0.2.6.tar.gz"
#SOURCE_VCFTOOLS="https://github.com/vcftools/vcftools/archive/v0.1.14.tar.gz"
SOURCE_VCFTOOLS="http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12a.tar.gz/download"

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
  if hash curl 2>/dev/null; then
    curl -sS -o $1.$EXT -L $2
  else
    wget -nv -O $1.$EXT $2
  fi
  mkdir -p $1
  `$DECOMP $1.$EXT`
  tar --strip-components 1 -C $1 -xf $1.tar
}

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /software/peakrescue"
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
cp -rp $INIT_DIR/perl/bin/hdr	$INST_PATH/bin/
cp -rp $INIT_DIR/README.md	$INST_PATH/README.md

cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
unset PERL5LIB
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"

# log information about this system
(
    echo '============== System information ===='
    set -x
    lsb_release -a
    uname -a
    sw_vers
    system_profiler
    grep MemTotal /proc/meminfo
    set +x
    echo
) >>$INIT_DIR/setup.log 2>&1

perlmods=( "File::ShareDir" "File::ShareDir::Install" )

set -e
for i in "${perlmods[@]}" ; do
  echo -n "Installing build prerequisite $i..."
  (
    set -x
    $INIT_DIR/perl/bin/cpanm -v --mirror http://cpan.metacpan.org -l $INST_PATH $i
    set +x
    echo; echo
  ) >>$INIT_DIR/setup.log 2>&1
  done_message "" "Failed during installation of $i."
done

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

cd $SETUP_DIR

CURR_TOOL="tabix-0.2.6"
CURR_SOURCE=$SOURCE_TABIX
echo -n "Building $CURR_TOOL ..."
if [ -e $SETUP_DIR/$CURR_TOOL.success ]; then
  echo -n " previously installed ..."
else
  (
    set -ex
    get_distro $CURR_TOOL $CURR_SOURCE
    cd $SETUP_DIR/$CURR_TOOL
    make -j$CPU
    cp tabix $INST_PATH/bin/.
    cp bgzip $INST_PATH/bin/.
    cd perl
    perl Makefile.PL INSTALL_BASE=$INST_PATH
    make
    make test
    make install
    touch $SETUP_DIR/$CURR_TOOL.success
  ) >>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build $CURR_TOOL."

# need to build samtools as has to be compiled in correct way for perl bindings
# does not deploy binary in this form
echo -n "Building samtools ..."
if [ -e $SETUP_DIR/samtools.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  (
  set -x
  if [ ! -e samtools ]; then
    get_distro "samtools" $SOURCE_SAMTOOLS
    perl -i -pe 's/^CFLAGS=\s*/CFLAGS=-fPIC / unless /\b-fPIC\b/' samtools/Makefile
  fi
  make -C samtools -j$CPU
	cp $SETUP_DIR/samtools/samtools $INST_PATH/bin/
  touch $SETUP_DIR/samtools.success
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build samtools."

export SAMTOOLS="$SETUP_DIR/samtools"

export PATH="$INST_PATH/bin:$PATH"

cd $SETUP_DIR

CURR_TOOL="vcftools"
CURR_SOURCE=$SOURCE_VCFTOOLS
echo -n "Building $CURR_TOOL ..."
if [ -e $SETUP_DIR/$CURR_TOOL.success ]; then
  echo -n " previously installed ..."
else
  (
    set -ex
    get_distro $CURR_TOOL $CURR_SOURCE
    cd $SETUP_DIR/$CURR_TOOL
    patch Makefile < $INIT_DIR/patches/vcfToolsInstLocs.diff
    patch perl/Vcf.pm < $INIT_DIR/patches/vcfToolsProcessLog.diff
    make -j$CPU PREFIX=$INST_PATH
    touch $SETUP_DIR/$CURR_TOOL.success
  ) >>$INIT_DIR/setup.log 2>&1
fi

done_message "" "Failed to build $CURR_TOOL."

cd "$INIT_DIR/perl"

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
(
  set -x
  perl $INIT_DIR/perl/bin/cpanm -v --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps . < /dev/null
  set +x
) >>$INIT_DIR/setup.log 2>&1
done_message "" "Failed during installation of core dependencies."

echo -n "Installing cgpVaf ..."
(
  set -e
  cd "$INIT_DIR/perl"
	echo -n `pwd`
  perl Makefile.PL INSTALL_BASE=$INST_PATH
  make
  make test
  make install
) >>$INIT_DIR/setup.log 2>&1
done_message "" " cgpVaf install failed."

# cleanup all junk
rm -rf $SETUP_DIR



echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo "  $PERLARCH"
echo

exit 0
