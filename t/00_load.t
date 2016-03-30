use Test::More tests => 3;

BEGIN { use_ok 'Bio::DB::HTS'; } 
BEGIN { use_ok 'Bio::DB::HTS::Tabix'; }
BEGIN { use_ok 'Bio::DB::HTS::Tabix::Iterator'; }

done_testing;

1;
