#-*-Perl-*-

use strict;
use ExtUtils::MakeMaker;
use File::Temp qw(tempfile);
use FindBin '$Bin';
use constant TEST_COUNT => 235;

use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch";

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if ($@) {
        use lib 't';
    }
    use Test;
    plan test => TEST_COUNT;
}

use Bio::DB::HTS;
use Bio::DB::HTS::AlignWrapper;


{
    my $sam = Bio::DB::Tam->open("$Bin/data/ex1.sam.gz");
    ok($sam);
    my $align = Bio::DB::Bam::Alignment->new();
    ok($align);

    # quench annoying stderr message from library here
    open my $saverr,">&STDERR";
    open STDERR,">/dev/null";
    my $head  = Bio::DB::Tam->header_read2("$Bin/data/ex1.fa.fai");
    open STDERR,">&",$saverr;
    ok($head);

    my $result = $sam->read1($head,$align);
    ok($result>0);
    ok($align->qseq,'CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCG');
    ok($align->start,1);
    ok($sam->read1($head,$align)>0);
    ok($align->start,3);
    ok($header->target_name->[$align->tid],'seq1');

    # test ability to write a BAM file
    my (undef,$filename) = tempfile('BAM_TEST_XXXXXXX',UNLINK=>1);
    $sam = Bio::DB::Tam->open("$Bin/data/ex1.sam.gz");
    $bam = Bio::DB::Bam->open($filename,'w');
    ok($bam);
    ok($bam->header_write($head),0);
    $count = 0;
    while ($sam->read1($head,$align) > 0) {
	$count++;
	$bam->write1($align);
    }
    ok($count,3307);

}
