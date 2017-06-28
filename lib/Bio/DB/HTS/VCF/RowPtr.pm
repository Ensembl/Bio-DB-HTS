package Bio::DB::HTS::VCF::RowPtr;

use base Bio::DB::HTS::VCF::Row;

sub DESTROY {
    # do nothing (overwrite subroutine in base class)
}

1;
