#!/bin/bash
export PERL5LIB=$PERL5LIB:$PWD/ensembl-test/modules:$PWD/lib:$PWD/blib/arch/auto/Bio/DB/HTS/:$PWD/blib/arch/auto/Bio/DB/HTS/Faidx

export TEST_AUTHOR=$USER

export WORK_DIR=$PWD

echo "Running test suite"
echo "Using PERL5LIB:$PERL5LIB"
echo "Current Directory list"
pwd
ls

echo "Test list"
ls -l t

echo "COVERALLS value=$COVERALLS"
echo "HTSLIB_VERSION value=$HTSLIB_VERSION"

perl $PWD/ensembl-test/scripts/runtests.pl t $SKIP_TESTS


rt=$?
if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit $rt
fi
