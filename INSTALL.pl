#!/usr/bin/perl

use strict;
use File::Temp 'tempdir';

my $use_prefix = shift @ARGV if(@ARGV > 0);

prompt_yn("This will install Bio-HTSTools and its dependencies. Continue?") or exit 0;
my $prefix_path = prompt_prefix($use_prefix);

# STEP 0: various dependencies
my $git = `which git`;
$git or die <<END;
'git' command not in path. Please install git and try again.
On Debian/Ubuntu systems you can do this with the command:

  apt-get install git
END


`which cc` or die <<END;
'cc' command not in path. Please install it and try again.
On Debian/Ubuntu systems you can do this with the command:

  apt-get install build-essential
END

`which make` or die <<END;
'make' command not in path. Please install it and try again.
On Debian/Ubuntu systems you can do this with the command:

  apt-get install build-essential
END

-e '/usr/include/zlib.h' or die <<END;
zlib.h library header not found in /usr/include. Please install it and try again.
On Debian/Ubuntu systems you can do this with the command:

  apt-get install zlib1g-dev
END
    ;

eval "require Bio::SeqFeature::Lite" or die <<END;
BioPerl does not seem to be installed. Please install it and try again.
On Debian/Ubuntu systems you can do this with the command:

    apt-get install bioperl

On other systems use the CPAN shell:

    perl -MCPAN -e 'install Bio::Perl'
END
    ;

# STEP 1: Create a clean directory for building
my $install_dir = tempdir(CLEANUP => 1);
info("Performing build in $install_dir");


# STEP 2: Check out HTSlib
info("Checking out HTSlib");
chdir $install_dir;
system "git clone https://github.com/samtools/htslib.git";
-d './htslib' or die "git clone seems to have failed. Could not find $install_dir/htslib directory";
chdir './htslib';
system "git checkout master";

# STEP 3: Check out Bio-HTS
info("Checking out Bio-HTS");
chdir $install_dir;
system "git clone https://github.com/Ensembl/Bio-HTS.git";
-d './Bio-HTS' or die "git clone seems to have failed. Could not find $install_dir/Bio-HTS directory";
chdir "./Bio-HTS";
system "git checkout master";

# Step 4: Build libhts.a
info("Building HTSlib");
chdir "$install_dir/htslib";
# patch makefile
open my $in, '<','Makefile'     or die "Couldn't open Makefile for reading: $!";
open my $out,'>','Makefile.new' or die "Couldn't open Makefile.new for writing: $!";
while (<$in>) {
    chomp;
    if (/^CFLAGS/ && !/-fPIC/) {
	s/#.+//;  # get rid of comments
	$_ .= " -fPIC -Wno-unused -Wno-unused-result";
    }
} continue {
    print $out $_,"\n";
}

close $in;
close $out;
rename 'Makefile.new','Makefile' or die "Couldn't rename Makefile.new to Makefile: $!";
if($use_prefix) {
  system "make";
  system "make install prefix=$prefix_path";
}
else {
  system "make";
}
-e 'libhts.a' or die "Compile didn't complete. No libhts.a library file found";

# Step 5: Build Bio::DB::HTSlib
info("Building Bio::DB::HTSlib");
chdir "$install_dir/Bio-HTS";
if($use_prefix) {
  system "env HTSLIB_DIR=$prefix_path/lib perl Build.PL --install_base=$prefix_path";
}
else {
  system "env HTSLIB_DIR=$install_dir/htslib perl Build.PL";
}
-e "./Build" or die "Build.PL didn't execute properly: no Build file found";
system "./Build";
`./Build test` =~ /Result: PASS/ or die "Build test failed. Not continuing";

# Step 6: Install
if($use_prefix) {
  info("Installing Bio-HTSTools to $prefix_path.");
  system "./Build install";
}
else {
  info("Installing Bio-HTSTools using sudo. You will be asked for your password.");
  info("If this step fails because sudo isn't installed, go back and run this script again as superuser.");
  system "sudo ./Build install";
}

# Step 7: Yay!
info("Bio-HTSTools is now installed.");
chdir '/';

exit 0;

sub prompt_prefix {
  my $run = shift;
  return 0 unless($run);
  print STDERR "Install to local/custom path requested [~/local/]: ";
  my $input = <>;
  chomp $input;
  $input= "$ENV{HOME}/local/" unless($input);
  return $input;
}

sub prompt_yn {
    my $msg = shift;
    print STDERR "$msg [Y/n]: ";
    my $input = <>;
    chomp $input;
    return 1 unless $input;
    return $input =~ /^[yY]/;
}

sub info {
    my $msg = shift;
    chomp $msg;
    print STDERR "\n*** $msg ***\n";
}
