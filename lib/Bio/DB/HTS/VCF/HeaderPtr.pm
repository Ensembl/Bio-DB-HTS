package Bio::DB::HTS::VCF::HeaderPtr;

use base Bio::DB::HTS::VCF::Header;

$Bio::DB::HTS::VCF::HeaderPtr::VERSION = '3.01';

sub DESTROY {
    # do nothing (overwrite subroutine in base class)
}

1;
