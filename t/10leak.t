#-*-Perl-*-
# Copyright [2015-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;
use File::Temp qw(tempfile);
use FindBin '$Bin';

use constant HAS_LEAKTRACE => eval{ require Test::LeakTrace };
use Test::More HAS_LEAKTRACE ? (tests => 6) : (skip_all => 'require Test::LeakTrace');
use Test::LeakTrace;

use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch";

use_ok('Bio::DB::HTS');
use_ok('Bio::DB::HTS::AlignWrapper');
use_ok('Bio::DB::HTS::VCF');

my $cramfile = "$Bin/data/yeast.sorted.cram";
my $fastafile = "$Bin/data/yeast.fasta";

# low level tests (defined in lib/Bio/DB/HTS.pm)
no_leaks_ok {
    my $hts_file = Bio::DB::HTSfile->open($cramfile);
    my $header  = $hts_file->header_read();
    my $targets = $header->n_targets;

    my $target_names = $header->target_name;
    my $target_lens = $header->target_len;
    my $text = $header->text;

    my $c = "\@CO\tThis is a comment\n";
    $header->text($c);

    my $region = 'X:51-1000' ;
    my $fai = Bio::DB::HTS::Fai->open("$Bin/data/yeast.fasta");
    my $seq = $fai->fetch($region);

    my $count;
    while ( my $b = $hts_file->read1($header) )
    {
        $count++;
    }

    my @result = $header->parse_region($region);
    @result = $header->parse_region('seq_invalid:51-1000');

    Bio::DB::HTSfile->index_build($cramfile);
    my $index = Bio::DB::HTSfile->index_load($hts_file);

    my @a;
    my $print_region = sub {
        my ( $alignment, $data ) = @_;
        push @a, $alignment;
        return;
    };

    $index->fetch( $hts_file, $header->parse_region('X'),
                   $print_region, "foobar" );

    my %matches;
    my $fetch_back = sub {
        my ( $tid, $pos, $p ) = @_;
        my $p1 = $pos + 1;
        my $r  = $fai->fetch( $header->target_name->[$tid] . ":$p1-$p1" );
        for my $pileup (@$p) {
            my $b    = $pileup->b;
            my $qpos = $pileup->qpos;
            my $base =
              $pileup->indel == 0 ? substr( $b->qseq, $qpos, 1 ) :
              $pileup->indel > 0 ? '*' :
              '-';
            $matches{matched}++ if $r eq $base;
            $matches{total}++;
        }
    };

    $index->pileup( $hts_file, $header->parse_region($region), $fetch_back );
} 'Low level';

# high level tests (defined in lib/Bio/DB/HTS.pm)

#
# Note that the high level API does not reset the CRAM file pointer to the start
# of the file as the method to do so is (at time or writing) not easily accessible.
# Therefore some of these tests create a new HTS object.
#

no_leaks_ok {
  for my $use_fasta (0,1) {
    my $hts = Bio::DB::HTS->new(
				-fasta => $fastafile,
				-bam          => $cramfile,
				-expand_flags => 1,
				-autoindex    => 1,
				-force_refseq => $use_fasta, );
    my $seg = $hts->segment('I');
    my $seq = $seg->seq;

    my @alignments =
      $hts->get_features_by_location( -seq_id => 'I',
                                      -start  => 500,
                                      -end    => 20000 );

    my @keys = $alignments[0]->get_all_tags;

    my %att = $alignments[0]->attributes;

    my $query = $alignments[0]->query;
    my $target = $alignments[0]->target;

    my @pads = $alignments[0]->padded_alignment;

    # There is an issue to be resolved still with fetching features by name for CRAM files
    # so opening a new object each time to ensure the tests take place.
    $hts = Bio::DB::HTS->new(
			     -fasta => $fastafile,
			     -bam          => $cramfile,
			     -expand_flags => 1,
			     -autoindex    => 1,
			     -force_refseq => $use_fasta, );

    my @f = $hts->features( -name => 'SRR507778.24982' );

    $hts = Bio::DB::HTS->new(
			     -fasta => $fastafile,
			     -bam          => $cramfile,
			     -expand_flags => 1,
			     -autoindex    => 1,
			     -force_refseq => $use_fasta, );
    @f = $hts->features(
			-filter => sub {
			  my $a = shift;
			  return 1 if $a->display_name eq 'SRR507778.24982';
			} );

    $hts = Bio::DB::HTS->new(
			     -fasta => $fastafile,
			     -bam          => $cramfile,
			     -expand_flags => 1,
			     -autoindex    => 1,
			     -force_refseq => $use_fasta, );
    @f = $hts->features(
			-seq_id => 'XIII',
			-filter => sub {
			  my $a = shift;
			  return 1 if $a->display_name =~ /^SRR507778/;
			} );

    $hts = Bio::DB::HTS->new(
			     -fasta => $fastafile,
			     -bam          => $cramfile,
			     -expand_flags => 1,
			     -autoindex    => 1,
			     -force_refseq => $use_fasta, );
    @f = $hts->features(
			-filter => sub {
			  my $a = shift;
			  return 1 if $a->display_name =~ /^SRR507778/;
			} );

    $hts = Bio::DB::HTS->new(
			     -fasta => $fastafile,
			     -bam          => $cramfile,
			     -expand_flags => 1,
			     -autoindex    => 1,
			     -force_refseq => $use_fasta, );
    @f = $hts->features( -name => 'SRR507778.24982',
                         -tags => { FIRST_MATE => 1 } );

    # try iteration
    my $i =
      $hts->get_seq_stream( -seq_id => 'XIII', -start => 200, -end => 10000 );
    my $count = 0;
    while ( $i->next_seq ) {
      $count++;
    }

    # try tam fh retrieval
    my $fh = $hts->get_seq_fh( -seq_id => 'XIII', -start => 200, -end => 10000, );
    $count = 0;
    $count++ while <$fh>;
    $fh->close;

    $hts = Bio::DB::HTS->new(
			     -fasta => $fastafile,
			     -bam          => $cramfile,
			     -expand_flags => 1,
			     -autoindex    => 1,
			     -force_refseq => $use_fasta, );
    $i = $hts->get_seq_stream(); # all features!
    $count = 0;
    while ( $i->next_seq ) {
      $count++;
    }

    $hts = Bio::DB::HTS->new(
			     -fasta => $fastafile,
			     -bam          => $cramfile,
			     -expand_flags => 1,
			     -autoindex    => 1,
			     -force_refseq => $use_fasta, );
    $i = $hts->get_seq_stream( -max_features => 200, -seq_id => 'XIII' );
    $count = 0;
    while ( $i->next_seq ) {
      $count++;
    }

    # try the read_pair aggregation
    my @pairs = $hts->features( -type => 'read_pair', -seq_id => 'XIII' );

    # test high level API version of pileup
    $hts = Bio::DB::HTS->new(
			     -fasta => $fastafile,
			     -bam          => $cramfile,
			     -expand_flags => 1,
			     -autoindex    => 1,
			     -force_refseq => $use_fasta, );
    my @coverage = $hts->features( -type => 'coverage', -seq_id => 'XIII' );
    my ($c) = $coverage[0]->get_tag_values('coverage');

    my %matches;
    my $fetch_back = sub {
      my ( $seqid, $pos, $p ) = @_;
      my $r = $hts->segment( $seqid, $pos, $pos )->dna;
      for my $pileup (@$p) {
	my $a    = $pileup->alignment;
	my $qpos = $pileup->qpos;
	my $dna  = $a->query->dna;
	my $base =
	  $pileup->indel == 0 ? substr( $dna, $qpos, 1 ) :
	  $pileup->indel > 0 ? '*' :
	  '-';
	$matches{matched}++ if $r eq $base;
	$matches{total}++;
      }
    };
    $hts->pileup( 'XIII:1-1000', $fetch_back );

    # test filtered coverage (no filtering)
    $hts = Bio::DB::HTS->new(
			     -fasta => $fastafile,
			     -bam          => $cramfile,
			     -expand_flags => 1,
			     -autoindex    => 1,
			     -force_refseq => $use_fasta, );
    @coverage = $hts->features( -type => 'coverage', -seq_id => 'XIII',
                                -filter => sub { return 1; } );
    ($c) = $coverage[0]->get_tag_values('coverage');

    # test filtered coverage (really filtering)
    $hts = Bio::DB::HTS->new(
			     -fasta => $fastafile,
			     -bam          => $cramfile,
			     -expand_flags => 1,
			     -autoindex    => 1,
			     -force_refseq => $use_fasta, );
    @coverage = $hts->features( -type => 'coverage', -seq_id => 'XIII',
                                -filter => sub { my $a = shift;
						 return defined $a->start && $a->start < 100000; } );
    ($c) = $coverage[0]->get_tag_values('coverage');

  } # end for my $use_fasta ( 0, ...)
} 'High level';

# get_info/format have been modified to allow not to specify ID,
# lots of allow/dealloc going underneath, test SVs are not leaking
no_leaks_ok {
  my $v = Bio::DB::HTS::VCF->new( filename => $Bin . "/data/test.bcf" );
  my $h = $v->header();
  my $row = $v->next();
  
  my $info_result = $row->get_info($h);
  $info_result = $row->get_info($h, 'DUMMY');
  $info_result = $row->get_info($h, 'NS');

  my $fmt_result = $row->get_format($h, 'DUMMY');
  $fmt_result = $row->get_format($h, 'HQ');
  $fmt_result = $row->get_format($h);

  # Query tests
  $v = Bio::DB::HTS::VCF->new( filename => $Bin . "/data/test.vcf.gz" );
  my $iter = $v->query("20:1000000-1231000");
  my $i = 0;
  $i++ while $row = $iter->next;

} 'VCF/BCF reading/querying';

done_testing();
