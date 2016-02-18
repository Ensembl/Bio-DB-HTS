package Bio::DB::HTS::VCF;

use Bio::DB::HTS; #load XS


sub new {
  my $class         = shift;
  my (%args) = @_;
  my $filename = $args{filename};
  my $warnings = $args{warnings};

  my $reader = bcf_sr_open($filename);
  die "Error getting reader" unless $reader;

  my $self = bless {
                    bcf_reader => $reader,
                    filename => $filename,
                   }, ref $class || $class;
  return $self;
}


sub num_variants {
    my $self = shift;
    return bcf_num_variants($self->{bcf_reader});
}

sub DEMOLISH {
    my $self = shift;

    if ( $self->{bcf_reader} ) {
        bcf_sr_close($self->{bcf_reader});
    }
}

1;
