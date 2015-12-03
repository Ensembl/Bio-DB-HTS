package Bio::DB::HTS::ReadIterator;

use strict;

sub new
{
    my $self = shift;
    my ($sam,$hts_file,$filter,$header) = @_;
    return bless {sam => $sam,
		  hts_file => $hts_file,
		  filter => $filter,
      header =>$header,
     },ref $self || $self;
}

sub next_seq
{
    my $self = shift;
    my $header = $self->{header} ;
    while (my $b = $self->{hts_file}->read1($header))
    {
      return Bio::DB::HTS::AlignWrapper->new($b,$self->{sam})
        if $self->{filter}->($b);
    }
    return;
}

1;
