package Bio::DB::HTS::VCF::RowPtr;

use base Bio::DB::HTS::VCF::Row;

$Bio::DB::HTS::VCF::RowPtr::VERSION = '2.11';

sub DESTROY {
    # do nothing (overwrite subroutine in base class)
}

1;
