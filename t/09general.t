#-*-Perl-*-
# Copyright [2015-2017] EMBL-European Bioinformatics Institute
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
use ExtUtils::MakeMaker;
use FindBin '$Bin';
use constant TEST_COUNT => 5;

use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch";

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if ($@) {
        use lib 't';
    }
    use Test;
    plan test => TEST_COUNT;
}

use Bio::DB::HTS;

{
    #no state required for these so call them as static methods
    #these get called in scalar context so the result is the number of
    #elements in the list rather than 'http' or whatever was matched
    ok(Bio::DB::HTS->is_remote( 'http://definitely.remote.com/' ), 1);
    ok(Bio::DB::HTS->is_remote( 'https://definitely.remote.com/' ), 1);
    ok(Bio::DB::HTS->is_remote( 'ftp://definitely.remote.com/' ), 1);
    ok(Bio::DB::HTS->is_remote( 's3://definitely-remote/' ), 1);
    ok(Bio::DB::HTS->is_remote( '/definitely/local/' ) == 0);
}

1;
