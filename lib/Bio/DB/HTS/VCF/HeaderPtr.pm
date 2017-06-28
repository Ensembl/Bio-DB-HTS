package Bio::DB::HTS::VCF::HeaderPtr;

use base Bio::DB::HTS::VCF::Header;

sub DESTROY {
    # do nothing (overwrite subroutine in base class)
}

1;
