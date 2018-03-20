#!/usr/bin/perl

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
use File::Temp 'tempdir';
use Cwd;
use File::Path qw(make_path);
use Getopt::Long;

my $htslib_version = "1.5";
my $version = get_version();

my $help = "INSTALL.pl [-h|--help] [--prefix=filepath] [--static] [~/prefix/path]\n";
$help .= "--help (-h)  - this help message\n";
$help .= "--prefix     - Path to install Bio::DB::HTS.\n";
$help .= "                  Alternative to providing a path at the end of comandline\n";
$help .= "                  If neither are provided defaults to root install (requires root access)\n";
$help .= "--static     - Build Bio::DB::HTS using the static libhts.a library and do not install htslib\n";
$help .= "--htslib_version     - Build Bio::DB::HTS using the specified HTSlib version (Default ".$htslib_version.")\n";

my $cwd = system 'pwd';


my $opts = parse_options();
my $prefix_path;
$prefix_path = $opts->{'prefix'} if(exists($opts->{'prefix'}) && defined($opts->{'prefix'}));
$htslib_version = $opts->{'htslib_version'} if($opts->{'htslib_version'}) ;

my ($version_major,$version_minor,$version_revision) = split /\./, $htslib_version ;


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

# The following libraries are version number dependant for HTSlib
if( $version_major >= 2 || ($version_major==1 && $version_minor>=5) )
{
-e '/usr/include/lzma.h' or die <<END;
lzma.h library header not found in /usr/include. Please install it and try again.
On Debian/Ubuntu systems you can do this with the command:

  apt-get install liblzma-dev
END

-e '/usr/include/bzlib.h' or die <<END;
zlib.h library header not found in /usr/include. Please install it and try again.
On Debian/Ubuntu systems you can do this with the command:

  apt-get install libbz2-dev
END
}
else
{
-e '/usr/include/zlib.h' or die <<END;
zlib.h library header not found in /usr/include. Please install it and try again.
On Debian/Ubuntu systems you can do this with the command:

  apt-get install zlib1g-dev
END

    ;
}


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


# STEP 2: Download HTSlib
info("Checking out HTSlib v".$htslib_version);
chdir $install_dir;
my $htslib_archive = $htslib_version.".zip" ;
my $htslib_archive_url = "https://github.com/samtools/htslib/archive/".$htslib_version.".zip" ;
system "wget -O " . $htslib_archive . " " .$htslib_archive_url  ;
if ( $? == -1 )
{
  die "HTSlib fetch command failed: $!\n";
}
-f './'.$htslib_archive or die "Could not fetch HTSlib archive ".$htslib_archive_url ;
system "unzip ".$htslib_archive ;
system "mv htslib-$htslib_version htslib" ;
-d './htslib' or die "Unzip seems to have failed. Could not find $install_dir/htslib directory";


# STEP 3: Download Bio-DB-HTS
info("Fetching version $version of Bio-DB-HTS from GitHub");
chdir $install_dir;
my $biodbhts_archive = "$version.zip" ;
my $biodbhts_archive_url = "https://github.com/Ensembl/Bio-DB-HTS/archive/$version.zip" ;
system "wget -O " . $biodbhts_archive . " " .$biodbhts_archive_url  ;
if ( $? == -1 )
{
  die "Bio-DB-HTS fetch command failed: $!\n";
}
-f './'.$biodbhts_archive or die "Could not fetch Bio::DB::HTS archive ".$biodbhts_archive_url ;
system "unzip ".$biodbhts_archive ;
system "mv Bio-DB-HTS-$version Bio-DB-HTS" ;
-d './Bio-DB-HTS' or die "Unzip seems to have failed. Could not find $install_dir/Bio-DB-HTS directory";


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
warn "***MAKE***\n";
system "make";
if(defined $prefix_path && !$opts->{'static'}) {
  warn "***MAKE INSTALL***\n";
  system "make install prefix=$prefix_path";
}
if($opts->{'static'}){
  system "rm -f $install_dir/htslib/libhts.so*";
}
-e 'libhts.a' or die "Compile didn't complete. No libhts.a library file found";

# Step 5: Build Bio::DB::HTS
info("Building Bio::DB::HTS");
chdir "$install_dir/Bio-DB-HTS";
my $cmd;
if(defined $prefix_path) {
  if(!$opts->{'static'}){
    $cmd = "env HTSLIB_DIR=$prefix_path/lib perl Build.PL --install_base=$prefix_path";
  }else{
    $cmd = "env HTSLIB_DIR=$install_dir/htslib perl Build.PL --install_base=$prefix_path";
  }
}else{
  $cmd = "env HTSLIB_DIR=$install_dir/htslib perl Build.PL";
}
if($opts->{'static'}){
  $cmd .= " --static=1";
}
warn "***CMD*** : $cmd\n";
system $cmd;
-e "./Build" or die "Build.PL didn't execute properly: no Build file found";
system "./Build";
`./Build test` =~ /Result: PASS/ or die "Build test failed. Not continuing";

# Step 6: Install
if(defined $prefix_path) {
  info("Installing Bio::DB::HTS to $prefix_path.");
  system "./Build install";
}
else {
  info("Installing Bio::DB::HTS using sudo. You will be asked for your password.");
  info("If this step fails because sudo isn't installed, go back and run this script again as superuser.");
  system "sudo ./Build install";
}
if($opts->{'static'}){
  system "rm -f $install_dir/htslib/libhts.a";
}

# Step 7: Yay!
info("Bio::DB::HTS v$version is now installed.");
chdir '/';

exit 0;

sub parse_options {
  my ($factory) = @_;
	my %opts = ();

  my $result = &GetOptions (
    'static' => \$opts{'static'},
    'prefix=s' => \$opts{'prefix'},
    'htslib_version=s' => \$opts{'htslib_version'},
    'h|help' => \$opts{'h'},
  );

  if($opts{'h'}){
    warn $help,"\n";
    exit(0);
  }

  if(defined($opts{'prefix'})){
    $opts{'prefix'} = prefix_install($opts{'prefix'});
  }elsif(!defined($opts{'prefix'}) && @ARGV > 0){
    $opts{'prefix'} = prefix_install(shift(@ARGV));
  }elsif(exists($ENV{PERL_LOCAL_LIB_ROOT}) && defined($ENV{PERL_LOCAL_LIB_ROOT})){
    $opts{'prefix'} = prefix_install($ENV{PERL_LOCAL_LIB_ROOT});
  }

  $opts{'static'} = 0 if(!defined($opts{'static'}));

  $opts{'static'} = $ENV{STATIC_HTS} if((!$opts{'static'}) &&
                                      exists($ENV{STATIC_HTS}) && defined($ENV{STATIC_HTS}));


  return \%opts;
}

sub info {
    my $msg = shift;
    chomp $msg;
    print STDERR "\n*** $msg ***\n";
}

sub prefix_install {
  my $prefix_path = shift;
  if($prefix_path =~ s/^\~//) {
    $prefix_path = $ENV{HOME}.$prefix_path
  }
  elsif($prefix_path !~ m|^/|) {
    $prefix_path = getcwd().'/'.$prefix_path;
  }
  my $err;
  unless(-e $prefix_path) {
    make_path($prefix_path, {verbose => 1, error => \$err});
    if (@$err) {
      for my $diag (@$err) {
        my ($dir, $message) = %$diag;
        if ($dir eq '') { print "general error: $message\n"; }
        else { print "problem creating dir $dir: $message\n"; }
      }
    }
  }
  return $prefix_path;
}

sub get_version {
  use FindBin '$Bin';
  my $version = `grep -m 1 "HTS::VERSION" $Bin/lib/Bio/DB/HTS.pm`;
  $version =~ /(\d+?\.\d+?)\D/;
  $version = $1;

  use Carp;
  croak "Couldn't get a meaningful library version: $version"
    unless $version =~ /\d+?\.\d+?/;

  return $version;
}
