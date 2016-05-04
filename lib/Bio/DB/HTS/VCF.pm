=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 AUTHOR

Alex Hodgkins
Rishi Nag E<lt>rishi@ebi.ac.ukE<gt>

=head1 NAME

Bio::DB::HTS::VCF -- Read VCF/BCF data files

=head1 DESCRIPTION

This module provides a Perl interface to the HTSlib library for
reading variant calls stored in VCF and BCF file databases.

Individual rows can be read in and fields accessed.

A sweep set of methods allows running through rows one by one, either backwards or
forwards through the file.

=head1 SYNOPSIS

=head2 VCF Row Objects

  use Bio::DB::HTS::VCF ;

  # VCF Row objects


=head2 VCF Sweep Objects

Open the file and process using sweeps. Note that the two methods maintain pointers that are
independant of one another. Using the next_row() will start at the first row in the file
and go on to the next row in subsequent reads. This is independant of previous_row() calls.
Similarly previous_row() will start at the last row and read backwards.

  use Bio::DB::HTS::VCF ;

  my $sweep = Bio::DB::HTS::VCFSweep->new(filename => "data/test.vcf.gz");
  $sweep->header;
  my $row_forwards = $sweep->next_row(); #returns first row in file
  my $row_backwards = $sweep->previous_row(); #returns last row in file
  my $row_forwards = $sweep->next_row(); # returns second row in file
  my $row_backwards = $sweep->previous_row(); #returns penultimate row in file



=cut

package Bio::DB::HTS::VCF;
$Bio::DB::HTS::VCF::VERSION = '1.12';

use Bio::DB::HTS;
use strict;
use warnings;

sub new {
  my $class         = shift;
  my (%args) = @_;
  my $filename = $args{filename};
  my $warnings = $args{warnings};
  die "Filename needed for VCF access" unless $filename;
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

package Bio::DB::HTS::VCFSweep ;
$Bio::DB::HTS::VCFSweep::VERSION = '1.12';

use Bio::DB::HTS;
use strict;
use warnings;

sub new
{
  my $class         = shift;
  my (%args) = @_;
  my $filename = $args{filename};
  my $warnings = $args{warnings};

  die "Filename needed for VCF sweep access" unless $filename;

  my $reader = sweep_open($filename);
  die "Error opening file for VCF sweep" unless $reader;
  my $header = $reader->header_read() ;
  my $self = bless {
                    sweep => $reader,
                    filename => $filename,
                    header => $header,
                   }, ref $class || $class;
  return $self;
}

sub header
{
  shift->{header} ;
}

sub next_row
{
    my $self = shift;
    my $b    = $self->{sweep};
    return $b->sweep_next();
}

sub previous_row
{
    my $self = shift;
    my $b    = $self->{sweep};
    return $b->sweep_previous();
}

sub close {
    my $self = shift;

    if ( $self->{sweep} ) {
        sweep_close($self->{sweep});
    }
}

1;
