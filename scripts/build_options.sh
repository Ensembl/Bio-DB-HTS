

#
# Script for testing various build configurations
# $1 - clone command -eg "git clone -b static-install https://github.com/drjsanger/Bio-HTS.git"
# $2 - the test to be run
#
export PERL5LIB_ORIG=$PERL5LIB

#
#test the Build.PL with various options
#
if [ "$2" = "BUILD_SYSTEM_INSTALLED_HTSLIB" ]; then
    echo Installs htslib, then runs Build process
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    sudo make install
    cd ..
    $1
    cd Bio-HTS
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
    echo Installs htslib, then runs Build process
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    sudo make install
    sudo ldconfig
    cd ..
    $1
    cd Bio-HTS
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
    echo Installs htslib to a local dir, then runs Build process
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    make prefix=~/localsw install
    export HTSLIB_DIR=
    cd ..
    $1
    cd Bio-HTS
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
    echo Makes htslib, then runs Build process. Should run from this location.
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    make
    export HTSLIB_DIR="$PWD"
    echo $HTSLIB_DIR
    cd ..
    $1
    cd Bio-HTS
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
    echo makes htslib, then runs Build process
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    make
    export HTSLIB_DIR_FOR_FLAG="$PWD"
    echo $HTSLIB_DIR_FOR_FLAG
    cd ..
    $1
    cd Bio-HTS
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

if [ "$2" = "BUILD_HTSLIB_DIR_WITH_STATIC_FLAG" ]; then
    echo Installs htslib, then runs Bio::DB::HTS Build process with static flag
    git clone -b master --depth=1 https://github.com/samtools/htslib.git
    cd htslib
    make
    rm -f libhts.so*
    export HTSLIB_DIR_FOR_FLAG="$PWD"
    echo $HTSLIB_DIR_FOR_FLAG
    cd ..
    $1
    cd Bio-HTS
    export HTSLIB_DIR=$HTSLIB_DIR_FOR_FLAG
    perl Build.PL --static=1 --install_base=~/localsw
    ./Build
    echo "Build of Bio::DB::HTS completed"
    rm -rf $HTSLIB_DIR_FOR_FLAG
    export PERL5LIB=$PERL5LIB:~/localsw
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
    cd Bio-HTS
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
    cd Bio-HTS
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
    cd Bio-HTS
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
    cd Bio-HTS
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
    cd Bio-HTS
    perl INSTALL.pl --static
    cd t
    for f in $(ls *.t) ;
    do
        perl $f
    done
    echo "Completed $2"
    exit 0
fi

echo Build test option $2 not found
