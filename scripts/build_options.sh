#!/bin/bash

# Copyright [2015-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# Script for testing/installing various build configurations of Bio::DB::HTS
# $1 - clone command for Bio::DB::HTS from GitHub
#    e.g. "git clone -b master https://github.com/Ensembl/Bio-DB-HTS.git"
#    Set this to be an alternative clone command if required
# $2 - the test to be run, to match one of the options below
#
export PERL5LIB_ORIG=$PERL5LIB

#
# Tests the Build.PL with various options
#
if [ "$2" = "BUILD_SYSTEM_INSTALLED_HTSLIB" ]; then
    echo Installs htslib, then runs Bio::DB::HTS Build process
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    sudo make install
    cd ..
    $1
    cd Bio-DB-HTS
    perl Build.PL
    ./Build
    export PERL5LIB=$PERL5LIB:$(pwd -P)/lib:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/Faidx
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    export PERL5LIB=$PERL5LIB_ORIG
    exit 0
fi

if [ "$2" = "BUILD_SYSTEM_INSTALL_ALL" ]; then
    echo Installs htslib, then builds and installs Bio::DB::HTS
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    sudo make install
    sudo ldconfig
    cd ..
    $1
    cd Bio-DB-HTS
    perl Build.PL
    sudo ./Build install
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    exit 0
fi


if [ "$2" = "BUILD_LOCAL_INSTALLED_HTSLIB" ]; then
    echo Installs htslib and Bio::DB::HTS to a local dir
    echo Specifies htslib dir using --prefix
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    make prefix=~/localsw install
    export HTSLIB_DIR=
    cd ..
    $1
    cd Bio-DB-HTS
    perl Build.PL --prefix=~/localsw
    ./Build
    export PERL5LIB=$PERL5LIB:$(pwd -P)/lib:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/Faidx
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    export PERL5LIB=$PERL5LIB_ORIG
    exit 0
fi



if [ "$2" = "BUILD_HTSLIB_DIR_ENV" ]; then
    echo Builds htslib, then runs Bio::DB::HTS Build process. Should run from this location.
    echo Specifies htslib location with HTSLIB_DIR environment variable
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    make
    export HTSLIB_DIR="$PWD"
    echo $HTSLIB_DIR
    cd ..
    $1
    cd Bio-DB-HTS
    perl Build.PL
    ./Build
    export PERL5LIB=$PERL5LIB:$(pwd -P)/lib:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/Faidx
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    export PERL5LIB=$PERL5LIB_ORIG
    exit 0
fi

if [ "$2" = "BUILD_HTSLIB_DIR_FLAG" ]; then
    echo Makes htslib, then runs Bio::DB::HTS Build process
    echo Specifies htslib location with --htslib flag
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    make
    export HTSLIB_DIR_FOR_FLAG="$PWD"
    echo $HTSLIB_DIR_FOR_FLAG
    cd ..
    $1
    cd Bio-DB-HTS
    perl Build.PL --htslib=$HTSLIB_DIR_FOR_FLAG
    ./Build
    export PERL5LIB=$PERL5LIB:$(pwd -P)/lib:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/Faidx
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    export PERL5LIB=$PERL5LIB_ORIG
    exit 0
fi


#TODO Alien::HTSlib dependency resolver
#TODO pkg-config test
#Test the static option using the INSTALL_STATIC_FLAG option

#
#test the INSTALL.pl script with various options
#

if [ "$2" = "INSTALL_WITH_SYSTEM_HTSLIB" ]; then
    echo INSTALL.pl with system install of htslib for running
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    sudo make install
    sudo ldconfig
    cd ..
    $1
    cd Bio-DB-HTS
    perl INSTALL.pl
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    exit 0
fi

if [ "$2" = "INSTALL_WITH_OTHER_HTSLIB" ]; then
    echo INSTALL.pl, with a htslib installed elsewhere for running
    export LD_LIBRARY_PATH_ORIG=$LD_LIBRARY_PATH
    git clone -b master --depth=1 https://github.com/samtools/htslib.git htslib_run_location
    cd htslib_run_location
    make
    cd ..
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd -P)/htslib_run_location
    $1
    cd Bio-DB-HTS
    perl INSTALL.pl
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ORIG
    echo "Completed $2"
    exit 0
fi

if [ "$2" = "INSTALL_PREFIX_PATH" ]; then
    echo INSTALL.pl with prefix at end of line
    $1
    cd Bio-DB-HTS
    perl INSTALL.pl ~/prefix_path_test
    export PERL5LIB=$PERL5LIB:~/prefix_path_test/lib/perl5/x86_64-linux-gnu-thread-multi/:~/prefix_path_test/lib/perl5/x86_64-linux-gnu-thread-multi/auto
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    export PERL5LIB=$PERL5LIB_ORIG
    exit 0
fi

if [ "$2" = "INSTALL_PREFIX_FLAG" ]; then
    echo INSTALL.pl with prefix at end of line
    $1
    cd Bio-DB-HTS
    perl INSTALL.pl --prefix=~/prefix_flag_test
    export PERL5LIB=$PERL5LIB:~/prefix_flag_test/lib/perl5/x86_64-linux-gnu-thread-multi/:~/prefix_flag_test/lib/perl5/x86_64-linux-gnu-thread-multi/auto
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    export PERL5LIB=$PERL5LIB_ORIG
    exit 0
fi


if [ "$2" = "INSTALL_STATIC_FLAG" ]; then
    echo INSTALL.pl with static option as flag
    $1
    cd Bio-DB-HTS
    perl INSTALL.pl --static
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    exit 0
fi

if [ "$2" = "INSTALL_HTSLIB_VERSION" ]; then
    echo INSTALL.pl built against a specific release of HTSlib
    $1
    cd Bio-DB-HTS
    perl INSTALL.pl --htslib_version 1.3
    export PERL5LIB=$PERL5LIB:$(pwd -P)/lib:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/Faidx
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    exit 0
fi


echo Build test option $2 not found
