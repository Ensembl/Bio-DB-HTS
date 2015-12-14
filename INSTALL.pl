#!/usr/bin/perl

use strict;
use File::Temp 'tempdir';

prompt_yn("This will install Bio::DB::HTS and its dependencies. Continue?") or exit 0;

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
info( 'Stage 1 Completed' ) ;

# STEP 2: Check out htslib
info("Checking out htslib");
chdir $install_dir;
system "git clone https://github.com/samtools/htslib.git";
-d './htslib' or die "git clone seems to have failed. Could not find $install_dir/htslib directory";
chdir './htslib';
system "git checkout master";
info( 'Stage 2 Completed' ) ;


# Step 3: Build libhts.a
info("Building htslib");
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
system "make";
-e 'libhts.a' or die "Compile didn't complete. No libhts.a library file found";
info( 'Stage 3 Completed' ) ;


# Step 4: Build Bio::DB::HTS
# currently this fails
info("Building Bio::DB::HTS");
system "git clone https://github.com/rishidev/Bio-HTS.git";
chdir "$install_dir/Bio-HTS" ;
system "env HTSLIB_DIR=$install_dir/htslib perl Build.PL";
info( 'Stage 4 Not Completed' ) ;
info( 'Please build htslib and then call perl Build.PL manually' ) ;

# Step 5: Install
#info("Installing Bio::DB::Sam using sudo. You will be asked for your password.");
#info("If this step fails because sudo isn't installed, go back and run this script again as superuser.");
#system "sudo ./Build install";

# Step 6: Yay!
#info("Bio::DB::HTS is now installed.");
chdir '/';

exit 0;

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
