# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
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

use strict;
use warnings;


use Cwd;
use File::Spec;
use File::Basename qw/dirname/;
use Test::More;
use Test::Warnings;
use Bio::EnsEMBL::Test::TestUtils;

if ( not $ENV{TEST_AUTHOR} ) {
  my $msg = 'Author test. Set $ENV{TEST_AUTHOR} to a true value to run.';
  plan( skip_all => $msg );
}


#chdir into the file's target & request cwd() which should be fully resolved now.
#then go back
my $file_dir = dirname(__FILE__);
my $original_dir = cwd();
chdir($file_dir);
my $cur_dir = cwd();
chdir($original_dir);
my $root = File::Spec->catdir($cur_dir, File::Spec->updir());

my $max_lines = 20;
my $count = 0;
my $found_embl_ebi_year = 0;
my $current_year = Time::Piece->new()->year();

my $notice_file = File::Spec->catfile($root, "NOTICE");
my $skip_copyright = undef;

# Annoingly I reproduce code from Bio::EnsEMBL::Test::TestUtils
# below - sub is_notice_file_good - because of the slightly different
# format of the NOTICE file here
# TO DO: re-write 'is_notice_file_good' to make it more flexible
if (-e $notice_file) {
    fail("$notice_file is not a file, cannot be read or is not a text file")
        unless (-f $notice_file && -r $notice_file && -T $notice_file);
    open my $fh, '<', $notice_file or die "Cannot open $notice_file: $!";
    while(my $line = <$fh>) {
        last if $count >= $max_lines;
        $found_embl_ebi_year = 1 
            if $line =~ /Copyright \[2015\-$current_year\] EMBL-European Bioinformatics Institute/;
    }
    if ($found_embl_ebi_year) {
        ok(1, "$notice_file has the correct Copyright year [2015-$current_year]");
        $skip_copyright = 1;
    } else {
        my $msg = "$notice_file is missing the correct Copyright year [2015-$current_year] in the first $max_lines lines";
        ok(0, $msg); 
    }
}

my @source_files = all_source_code(File::Spec->catfile($root));
#Find all files & run
foreach my $f (@source_files) {
    next if $f =~ /\/travisci\//;
    next if $f =~ /\/(kseq_)?data\//;
    next if $f =~ /\/blib\//;
    next if $f =~ /\/_build\//;
    next if $f =~ /\/ppport\.h$/;
    next if $f =~ /\/CLEAN\b/;
    next if $f =~ /\.(tmpl|hash|nw|ctl|txt|html|textile|md)$/;
    has_apache2_licence($f, $skip_copyright);
}

done_testing();
