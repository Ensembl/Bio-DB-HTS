#!/bin/bash
export PERL5LIB=$PERL5LIB:$PWD/ensembl-test/modules:~/biodbhts/lib/perl5/x86_64-linux/

export TEST_AUTHOR=$USER

export WORK_DIR=$PWD

echo "Running test suite"
echo "Using $PERL5LIB"

echo "Test list"
pwd
ls -l t

echo "COVERALLS value=$COVERALLS"

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