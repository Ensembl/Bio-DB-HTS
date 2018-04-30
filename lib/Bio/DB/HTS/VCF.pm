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

=head1 AUTHORS

Rishi Nag E<lt>rishi@ebi.ac.ukE<gt>, original author.

Alessandro Vullo C<< <avullo at cpan.org> >>, the current developer and maintainer.

=head1 NAME

Bio::DB::HTS::VCF -- Read VCF/BCF data files

=head1 DESCRIPTION

This module provides a Perl interface to the HTSlib library for reading variant calls 
stored in VCF and BCF file databases.

The functions provided are for opening a VCF/BCF file, reading header values, querying 
specific chromosome intervals and then reading row values.

A sweep set of methods allows running through rows one by one, either backwards or
forwards through the file.

=head1 SYNOPSIS

  use Bio::DB::HTS::VCF ;

  ### File Open ###
  my $v = Bio::DB::HTS::VCF->new( filename => "/path/to/myfile.vcf.gz" );

  # once the file has been opened, various global values can be read from the header
  my $h = $v->header();

  $h->get_seqnames() ;
  $h->version() ;  # read the VCF file version

  $h->num_samples() ;
  $h->get_sample_names() ; #return an array of sample names

  $h->num_seqnames() ;
  $h->get_seqnames() ; # return an array of sequence names

  ### Individual rows can be read in and fields accessed ###
  my $row = $v->next() ;

  # row functions
  $row->chromosome($h) ;
  $row->position() ;
  $row->id() ;
  $row->num_filters() ;
  $row->quality() ;

  # retrieve alleles
  my $num_alleles = $row->num_alleles();
  my $alleles = $row->get_alleles();
  my $allele_index = 1;
  for my $a (@$alleles) {
    printf( "(%s, %s)\n", $a, $row->get_variant_type($allele_index++) ) ;
  }

  # query filters
  $row->has_filter($h,"DP50");
  $row->has_filter($h,"."); # PASS filter check

  $row->get_info_type($h, "AF"); # one of "String", "Integer", "Float" or "Flag".
  $info_result = $row->get_info($h, "NS"); # [3]

  $row->get_format_type($h, "GT") ; # "String"
  $row->get_format($h, "DP") ; # [ 1, 8, 5 ]

  ### free memory associated with the row
  Bio::DB::HTS::VCF::Row->destroy($row);

  ### query specific locations
  my $iter = $v->query("20:1000000-1231000");
  while (my $result = $iter->next) {
    print $result->chromosome($h), $result->position(), $result->id(), $result->num_filters(), $result->quality(), "\n";
  }

=head1 METHODS

=head2 C<new>

Opens a VCF/BCF file for reading. If the file is indexed (i.e. tabix for VCF, csi for BCF) 
the index is opened and used for querying arbitrary locations on the chromosomes.

=over 1

=item $vcf = Bio::DB::HTS::VCF->new($filepath)

Returns an instance of Bio::DB::HTS::VCF.

=back

=head2 C<header>

Returns instance of Bio::DB::HTS::VCF::Header, representing the header of the file.

=over 1

=item $header = $vcf->header()

=back

=head2 C<num_variants>

Returns the number of variants (i.e. rows) of the file.

=over 1

=item $nv = $vcf->num_variants();

=back

=head2 C<close>

Close the VCF/BCF file, allocated memory will be released, included the index, if present.

=over 1

=item $vcf->close()

=back

=head2 C<next>

Returns the next row (starting from the first one) read from the file.

=over

=item $row = $vcf->next()

Returns an instance of Bio::DB::HTS::VCF::Row or undef if end of file is reached.

=back

=head1 Querying an indexed VCF/BCF file

If the file is indexed, the file can be queried for variants on a specified region.
Regions can be specified using either the "chr", "chr:start" or "chr:start-end" format,
with start <= end.

Once an iterator is obtained, individual rows belonging to the result set can be
sequentially accessed by iteratively invoking the iterator next method until it
returns nothing.

=head2 C<query>

=over 1

=item $iterator = $vcf->query($region);

Returns an instance of Bio::DB::HTS::VCF::Iterator or undef if the chromosome is not 
found in the index or raises an exception in case the underlying HTSlib library cannot 
parse the region.

=back

=head1 HEADER METHODS

Once the file has been opened, various global values can be read from the header.

=head2 C<version>

Returns the VCF file version, as a string

=over 1

=item $h->version()

=back

=head2 C<num_samples>

Returns the number of samples

=over 1

=item $h->num_samples()

=back

=head2 C<get_sample_names>

Returns the list of sample names

=over 1

=item $sample_names = $h->get_sample_names()

Returns an array ref of strings representing the sample names

=back

=head2 C<get_seqnames>

Returns the number of sequence names

=over 1

=item $h->num_seqnames()

=back

=head2 C<get_seqnames>

Returns the list of sequence names

=over 1

=item $h->get_seqnames()

Returns an array ref of strings representing the sequence names

=back

=head2 C<fmt_text>

Get header formatted text, as a string

=over 1

=item $h->fmt_text()

Returns the text string representing the content of the header

=back

=head1 ROW METHODS

Individual rows can be read in and fields accessed. To read a row use the next function,
which returns a Bio::DB::HTS::VCF::Row instance.

Various fields can then be read from the row object. Some of the functions to read these 
fields will need the header object supplied.

=over 6

=item $row->print($header)

Returns a formatted textual representation of the row.

=item $row->chromosome($header)

=item $row->position()

=item $row->id()

=item $row->quality()

=item $row->reference()

=back

=head2 Accessing alleles information

=over 2

=item $row->num_alleles()

Returns the number of alleles

=item $row->get_alleles()

Returns the alleles as strings in an array ref

=back

The variant type of an allele can be determined using the index of the allele. The index starts
from 1 for the first allele:

=over 2

=item $row->is_snp()

Returns a true value if the row refers to a SNP.

=item $row->get_variant_type($allele_index)

This will return one of the values as defined in htslib. As of v1.3.1 these are as follows.

=over 5

=item VCF_REF   0

=item VCF_SNP   1

=item VCF_MNP   2

=item VCF_INDEL 4

=item VCF_OTHER 8

=back

=back

=head2 Row filters

Each row object has filters that may or may not have been applied to it.

=over 2

=item $row->num_filters()

Returns the number of filters of the row.

=item $row->has_filter($header, $filter)

Returns 0 if the filter is not present, 1 if it is present. The PASS filter is represented by a dot.

=back

=head2 Accessing info fields

Each row may have additional info fields associated with each allele in the row. 

=over 2

=item $row->get_info_type($header, $info_id)

Returns the type of the info ID as specified in the VCF file header, one of "String", "Integer", "Float" or "Flag".

=item $row->get_info($header, $info_id) or $row->get_info($header)

If an info_id string is passed, returns an array ref of values for that particular info field, 
one for each allele in the row. If the row does not have an item of that info, or it does not 
exist in the file, a string "ID_NOT_FOUND" will be returned.

Alternatively, the get_info() method can be invoked by just passing the header. In this case,
the whole info field is returned organised as a hash ref where keys are the info IDs and values
are the info fields for the corresponding ID.

=back

=head2 Accessing format fields

Formats are dealt with similarly to info fields.

=over 2

=item $row->get_format_type($header, $format_id)

Returns the type of the format ID as specified in the VCF file header, one of "String", "Integer", "Float" or "Flag".

=item $row->get_format($header, $format_id) or $row->get_format($header)

If a format_id string is passed, returns an array ref of values for that particular format ID. 
If the row does not have an item of that format, or it does not exist in the file, a string 
"ID_NOT_FOUND" will be returned.

Alternatively, the get_format() method can be invoked by just passing the header. In this case,
it returns the complete format specification as a hash ref of FORMAT_ID => [ FORMAT_ID_VALUE, ... ].

=back

=head2 Accessing genotypes

Genotype records are currently returned as a series of integers, across all the samples for the row.

=over 1

=item $row->get_genotypes($header)

Returns an array reference of integers representing genotype records.

=back

=head1 VCF SWEEP OBJECTS

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
$Bio::DB::HTS::VCF::VERSION = '2.11';

use Bio::DB::HTS;
use Scalar::Util qw/reftype/;
use strict;
use warnings;
use Carp 'croak';

use Bio::DB::HTS::VCF::Iterator;

sub new {
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
		    # can either have a tabix indexed VCF or a (CSI) indexed BCF file
		    index => Bio::DB::HTS::VCFfile->tbx_index_load($filename) || Bio::DB::HTS::VCFfile->bcf_index_load($filename)
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

sub query {
  my ($self, $region) = @_;

  die "Please provide a region" unless defined $region;

  my ( $chr, $start, $end ) = $region =~ /^([^:]+)(?::(\d+))?(?:-(\d+))?$/;
  die "You must specify a region in the format chr, chr:start or chr:start-end"
    unless defined $chr;
  die "Couldn't get start of the region" unless defined $start;
  die "End in $region is less than the start" if $end < $start;

  my $iter = Bio::DB::HTS::VCFfile->query($region, $self->{index}, $self->{header});
  unless ( $iter ) {
    # this likely means the chromosome wasn't found in the index, or it couldn't parse the provided region.
    my $seqnames_hash = { map { $_ => 1 } @{ $self->{header}->get_seqnames } };

    if ( not exists $seqnames_hash->{ $chr } ) {
      # warn("Specified chromosome '$chr' does not exist in file " . $self->{_filename});
      return undef;
    } else {
      die "Unable to get iterator for region '$region' in file ". $self->{filename} . " -- htslib couldn't parse your region string";
    }
  }

  return Bio::DB::HTS::VCF::Iterator->new( iter => $iter, file => $self->{vcf_file}, header => $self->{header}, index => $self->{index} );

}

sub close
{
    my $self = shift;
    if ( $self->{vcf_file} ) {
        $self->{vcf_file}->vcf_close();
        delete $self->{vcf_file};
      }
    if ($self->{tbx_index}) {
      Bio::DB::HTS::Tabix::tbx_close($self->{tbx_index});
      delete $self->{tbx_index};
    }
    if ($self->{bcf_index}) {
      Bio::DB::HTS::VCFfile->bcf_index_close($self->{bcf_index});
      delete $self->{bcf_index};
    }
    #
}

sub DESTROY {
    my $self = shift;
    return if reftype($self) ne 'HASH';
    $self->close();
    return;
}

package Bio::DB::HTS::VCF::Sweep ;
$Bio::DB::HTS::VCF::Sweep::VERSION = '2.11';

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
