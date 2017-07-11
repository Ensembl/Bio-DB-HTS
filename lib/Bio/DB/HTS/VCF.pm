=head1 LICENSE

Copyright [2015-2017] EMBL-European Bioinformatics Institute

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


=head2 VCF use

The functions provided are for opening a VCF/BCF file, reading header values and then reading row values.

  use Bio::DB::HTS::VCF ;

  # File Open and header read functions
  my $v = Bio::DB::HTS::VCF->new( filename => "/path/to/myfile.vcf.gz" ) ;
  my $h = $v->header();

Once the file has been opened, various global values can be read from the header.

  $h->get_seqnames() ;
  $h->version() ;  # read the VCF file version

  $h->num_samples() ;
  $h->get_sample_names() ; #return an array of sample names

  $h->num_seqnames() ;
  $h->get_seqnames() ; # return an array of sequence names

To read a row use the next function. Various other fields can then be read from the row object. Some
of the functions to read these fields will need the header object supplied.

  my $row = $v->next() ;

  # row functions
  $row->chromosome($h) ;
  $row->position() ;
  $row->id() ;
  $row->num_filters() ;
  $row->quality() ;

The alleles can be retrieved with the following functions. get_alleles() will return the alleles in an array
reference.
  my $num_alleles = $row->num_alleles() ;
  my $alleles_ref = $row->get_alleles() ;
  my @alleles = @$alleles_ref ;
  for my $a (@alleles)
  {
    print($a.",") ;
  }

The variant type of an allele can be determined using the index of the allele. The index starts
from 1 for the first allele.
   $row->get_variant_type($allele_index)

This will return one of the values as defined in htslib. As of v1.3.1 these are as follows.
   #define VCF_REF   0
   #define VCF_SNP   1
   #define VCF_MNP   2
   #define VCF_INDEL 4
   #define VCF_OTHER 8

Each row object has filters that may or may not have been applied to it.
has_filter will return 0 if the filter is not present, 1 if it is present.
The PASS filter is represented by a dot.
  $row->has_filter($h,"DP50") ;
  $row->has_filter($h,".") ; #PASS filter check

Each row may also have additional info fields associated with each allele in the row.
Each info tag is specified as a string. If the row does not have an item of that info, or
it does not exist in the file, a string ID_NOT_FOUND will be returned by the get_info() function.

These are returned in an array reference. There
will be one value per allele in the row. The get_info_type() function can be used to clarify
the type of info. This will be the type of the info as specified in the VCF file header,
one of String, Integer, Float or Flag). The flag will be signified by a value of 1.

  $row->get_info_type($h,"INFO_ID") ;
  $info_result = $row->get_info($h,"INFO_ID") ;

Formats are dealt with similarly. A string ID_NOT_FOUND will be returned by the get_format() function
if the specified format ID is not found.

  $row->get_format_type($h,"FORMAT_ID") ;
  $row->get_format($h,"FORMAT_ID") ;

Genotype records are currently returned as a series of integers, across all the samples for the row.
  $row->get_genotypes($h)

Finally, free memory associated with the row.
  Bio::DB::HTS::VCF::Row->destroy($row) ;


=head2 VCF Sweep Objects

Open the file and process using sweeps. Note that the two methods maintain pointers that are
independant of one another. Using the next_row() will start at the first row in the file
and go on to the next row in subsequent reads. This is independant of previous_row() calls.
Similarly previous_row() will start at the last row and read backwards. However a call to next_row()
is needed beforehand as the read fails otherwise.

At the time of writing there are issues which seem to be due to the underlying HTSlib API calls,
so using the next() function is preferable to using sweeps.

  use Bio::DB::HTS::VCF ;

  my $sweep = Bio::DB::HTS::VCFSweep->new(filename => "data/test.vcf.gz");
  $sweep->header;
  my $row_forwards = $sweep->next_row(); #returns first row in file
  my $row_backwards = $sweep->previous_row(); #returns last row in file
  my $row_forwards = $sweep->next_row(); # returns second row in file
  my $row_backwards = $sweep->previous_row(); #returns penultimate row in file



=cut

package Bio::DB::HTS::VCF;
$Bio::DB::HTS::VCF::VERSION = '2.9';

use Bio::DB::HTS;
use Scalar::Util qw/reftype/;
use strict;
use warnings;
use Carp 'croak';

sub new
{
  my $class         = shift;
  my (%args) = @_;
  my $filename = $args{filename};
  my $warnings = $args{warnings};
  die "Filename needed for VCF access" unless $filename;
  my $reader = Bio::DB::HTS::VCFfile->open($filename) or
    croak "Error getting VCF file reader: $!" ;
  my $h = $reader->header_read() or
    croak "Error getting VCF file header: $!" ;

  my $self = bless {
                    vcf_file => $reader,
                    filename => $filename,
                    header => $h,
                   }, ref $class || $class;
  return $self;
}


sub header
{
    my $self = shift;
    return $self->{header};
}

sub next
{
    my $self = shift;
    return $self->{vcf_file}->read1($self->{header});
}

sub num_variants
{
    my $self = shift;
    return $self->{vcf_file}->num_variants($self->{filename});
}

sub close
{
    my $self = shift;
    if ( $self->{vcf_file} )
    {
        $self->{vcf_file}->vcf_close();
        delete $self->{vcf_file};
    }
}

sub DESTROY {
    my $self = shift;
    return if reftype($self) ne 'HASH';
    $self->close();
    return;
}

package Bio::DB::HTS::VCF::Sweep ;
$Bio::DB::HTS::VCF::Sweep::VERSION = '2.9';

use Bio::DB::HTS;
use Scalar::Util qw/reftype/;
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

sub close
{
    my $self = shift;

    if ( $self->{sweep} )
    {
        sweep_close($self->{sweep});
        delete $self->{sweep};
    }
}

sub DESTROY {
    my $self = shift;
    return if reftype($self) ne 'HASH';
    $self->close();
    return;
}

1;
