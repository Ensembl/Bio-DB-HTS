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

Faidx -- Perl extension for accessing bgzip compressed and indexed FASTA using htslib

=head1 SYNOPSIS

 #include the module
 use Bio::DB::HTS::Faidx;

 #create the index object
 my $fasta = "$Bin/data/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa.gz" ;
 my $location = "I:1-100" ;
 my $seq_id = "I" ;
 my $index = Bio::DB::HTS::Faidx->new($fasta);

 #get sequence using the location
 my $seq = "" ;
 my $length = 0 ;
 ($seq, $length) = $index->get_sequence($location);

 $seq = $index->get_sequence_no_length($location);
 $length = $index->length($seq_id);

 #get sequence using the separate parameters.
 #Note here that the sequence start and end points are zero indexed, so code accordingly.
 ($seq, $length) = $index->get_sequence2("I",1,99);
 $seq = $index->get_sequence2_no_length("I",1,99);

 my @seq_ids = $index->get_all_sequence_ids();

 #returns 1 if sequence ID is present, 0 if not
 my $has_seq = $index->has_sequence('I');

=head1 AUTHOR

Rishi Nag

=cut


package Bio::DB::HTS::Faidx;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use Faidx ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(

) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(

);

our $VERSION = '2.11';

require XSLoader;
XSLoader::load('Bio::DB::HTS::Faidx', $VERSION);

# Preloaded methods go here.
1;
__END__
