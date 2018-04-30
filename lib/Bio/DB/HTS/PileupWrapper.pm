
=head1 LICENSE

Copyright [2015-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 NAME

Bio::DB::HTS::PileupWrapper -- Add high-level methods to Bio::DB::HTS::Pileup

=head1 SYNOPSIS

See L<Bio::DB::HTS/The generic fetch() and pileup() methods> for usage of the pileup() method.

=head1 DESCRIPTION

See L<Bio::DB::HTS::Pileup> for documentation of this object's
methods. This class is used by the high-level API to return
Bio::DB::HTS::AlignWrapper objects from the call to alignment() rather
than Bio::DB::HTS::Alignment.

=head1 AUTHOR

Rishi Nag E<lt>rishi@ebi.ac.uk<gt>

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::HTS>, L<Bio::DB::HTS::Constants>

=cut

package Bio::DB::HTS::PileupWrapper;
$Bio::DB::HTS::PileupWrapper::VERSION = '2.11';

use strict;
use warnings;

use Bio::DB::HTS::AlignWrapper;

our $AUTOLOAD;
use Carp 'croak';

sub new {
    my $package = shift;
    my ( $align, $sam ) = @_;
    return bless { sam => $sam, pileup => $align }, ref $package || $package;

}

sub AUTOLOAD {
    my ( $pack, $func_name ) = $AUTOLOAD =~ /(.+)::([^:]+)$/;
    return if $func_name eq 'DESTROY';

    no strict 'refs';
    $_[0] or die "autoload called for non-object symbol $func_name";
    croak qq(Can't locate object method "$func_name" via package "$pack")
      unless $_[0]->{pileup}->can($func_name);

    *{"${pack}::${func_name}"} = sub { shift->{pileup}->$func_name(@_) };

    shift->$func_name(@_);
}

sub can {
    my $self = shift;
    return 1 if $self->SUPER::can(@_);
    return $self->{pileup}->can(@_);
}

sub alignment {
    my $self = shift;
    return Bio::DB::HTS::AlignWrapper->new( $self->{pileup}->b, $self->{sam} );
}

1;
