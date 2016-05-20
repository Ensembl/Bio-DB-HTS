package include::TestHelper;
use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(test_output_file test_output_dir test_input_file);

=head2 test_output_file

 Title   : test_output_file
 Usage   : my $output_file = test_output_file();
 Function: Get the full path of a file suitable for writing to.
           When your test script ends, the file will be automatically deleted.
 Returns : string (file path)
 Args    : none

=cut

sub test_output_file {
    die "test_output_file takes no args\n" if @_;
    my $tmp = File::Temp->new();
    close($tmp); # Windows needs this
    return $tmp->filename;
}

=head2 test_output_dir

 Title   : test_output_dir
 Usage   : my $output_dir = test_output_dir();
 Function: Get the full path of a directory suitable for storing temporary files
           in.
           When your test script ends, the directory and its contents will be
           automatically deleted.
 Returns : string (path)
 Args    : none

=cut

sub test_output_dir {
    die "test_output_dir takes no args\n" if @_;

    return tempdir(CLEANUP => 1);
}

=head2 test_input_file

 Title   : test_input_file
 Usage   : my $input_file = test_input_file();
 Function: Get the path of a desired input file stored in the standard location
           (currently t/data), but correct for all platforms.
 Returns : string (file path)
 Args    : list of strings (ie. at least the input filename, preceded by the
           names of any subdirectories within t/data)
           eg. for the file t/data/in.file pass 'in.file', for the file
           t/data/subdir/in.file, pass ('subdir', 'in.file')

=cut

sub test_input_file {
  return File::Spec->catfile('t', 'kseq_data', @_);
}

1;