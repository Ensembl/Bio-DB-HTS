#!/usr/bin/perl

use strict;
use File::Temp 'tempdir';

prompt_yn("This will install Bio::DB::Sam and its dependencies. Continue?") or exit 0;

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


# STEP 2: Check out samtools
info("Checking out samtools 0.1.19");
chdir $install_dir;
system "git clone https://github.com/samtools/samtools.git";
-d './samtools' or die "git clone seems to have failed. Could not find $install_dir/samtools directory";
chdir './samtools';
system "git checkout 0.1.19";

# STEP 3: Check out Bio-SamTools
info("Checking out Bio-SamTools");
chdir $install_dir;
system "git clone https://github.com/GMOD/GBrowse-Adaptors.git";
-d './GBrowse-Adaptors' or die "git clone seems to have failed. Could not find $install_dir/GBrowse-Adaptors directory";
chdir "./GBrowse-Adaptors/Bio-SamTools";
system "git checkout release-1_39";

# Step 4: Build libbam.a
info("Building samtools");
chdir "$install_dir/samtools";
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
-e 'libbam.a' or die "Compile didn't complete. No libbam.a library file found";

# Step 5: Build Bio::DB::Sam
info("Building Bio::DB::Sam");
chdir "$install_dir/GBrowse-Adaptors/Bio-SamTools";
system "env SAMTOOLS=$install_dir/samtools perl Build.PL";
-e "./Build" or die "Build.PL didn't execute properly: no Build file found";
system "./Build";
`./Build test` =~ /Result: PASS/ or die "Build test failed. Not continuing";

# Step 6: Install
info("Installing Bio::DB::Sam using sudo. You will be asked for your password.");
info("If this step fails because sudo isn't installed, go back and run this script again as superuser.");
system "sudo ./Build install";

# Step 7: Yay!
info("Bio::DB::Sam is now installed.");
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
