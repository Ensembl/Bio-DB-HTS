#!/usr/bin/perl

use strict;
use Module::Build;
use Module::Load::Conditional qw(can_load);

my $HeaderFile = "hts.h";
my $LibFile    = "libhts.a";
my $ReadLine;

my ($sam_include,$sam_lib) = find_sam(); # may exit with error here

my $class = Module::Build->subclass(code=><<EOF);


EOF


my $build = $class->new(
    module_name        => 'Bio::HTSTools',
    dist_version_from  => 'lib/Bio/DB/HTS.pm',
    dist_author        => 'Rishi Nag',
    dist_abstract      => 'Perl interface to HTS library for DNA sequencing',
    license            => 'Apache Licence 2.0',
    include_dirs       => [$sam_include],
    extra_linker_flags => ["-L$sam_lib",'-lhts','-lpthread','-lz'],

    extra_compiler_flags=>[

        # must match DFLAGS in Samtools Makefile
	'-D_IOLIB=2','-D_FILE_OFFSET_BITS=64',

	# turn off warnings originating in Perl's Newx* calls
	'-Wformat=0',
    ],


    build_requires => {
      'ExtUtils::CBuilder' => 0,
    },
    requires     => {
	'perl'                => '5.008',
	'Bio::Root::Version'  => '1.006009001',
    },
#    create_makefile_pl => 'passthrough',
    );

$build->create_build_script;

exit 0;

sub find_sam
{
    my ($sam_include,$sam_lib);

    if (my $samtools = _samtools())
    {
      $sam_include = $samtools if -e "$samtools/$HeaderFile";
      $sam_include = "$samtools/include" if -e "$samtools/include/$HeaderFile";
      $sam_include = "$samtools/htslib" if -e "$samtools/htslib/$HeaderFile";
      $sam_lib     = $samtools if -e "$samtools/$LibFile";
      $sam_lib     = "$samtools/lib" if -e "$samtools/lib/$LibFile";
    }

    my @search_path = qw(/ /usr /usr/share /usr/local);

    unless ($sam_include)
    {
      for my $p (@search_path)
       {
         $sam_include ||= "$p/include" if -e "$p/include/$HeaderFile";
       }
    }

    unless ($sam_lib)
    {
      for my $p (@search_path)
      {
        $sam_lib ||= "$p/lib" if -e "$p/lib/$LibFile";
      }
    }


    unless ($sam_include && $sam_lib)
    {
      print STDOUT "This module requires htslib (http://htslib.org).\n";
      my $prompt = "Please enter the location of the $HeaderFile and compiled $LibFile files: ";
      my $found;
      while (!$found)
      {
        my $path = prompt($prompt);
        print STDOUT "\n";
        last unless $path;
        $sam_include = $path if -e "$path/$HeaderFile";
        $sam_include = "$path/include" if -e "$path/include/$HeaderFile";
        $sam_include = "$path/htslib" if -e "$path/htslib/$HeaderFile";
        $sam_lib = $path if -e "$path/$LibFile";
        $sam_lib = "$path/lib" if -e "$path/lib/$LibFile";
        $found = $sam_include && $sam_lib;
        unless ($found)
        {
          print STDOUT "That didn't seem to be right.\n";
          $prompt = "Try again, or hit <enter> to cancel: ";
        }
      }
    }

    unless ($sam_include && $sam_lib)
    {
      die <<END;
Can\'t find $LibFile and/or $HeaderFile!

If you haven\'t done so already, please compile htslib from
http://htslib.org and set the HTSLIB_DIR
environment variable to point to a  distribution directory
containing the compiled $LibFile and $HeaderFile files.
END
    }

    print STDOUT "Found $sam_include/$HeaderFile and $sam_lib/$LibFile\n";
    return ($sam_include,$sam_lib);
}


sub prompt {
    my $msg = shift;

    unless (defined $ReadLine) {
	eval "require Term::ReadLine";
	$ReadLine   = Term::ReadLine->can('new') || 0;
	$ReadLine &&= Term::ReadLine->new(\*STDOUT);
	eval {readline::rl_set('TcshCompleteMode','On')};
    }


    unless ($ReadLine) {
	print STDOUT $msg;
	chomp (my $in = <>);
	return $in;
    }

    my $in = $ReadLine->readline($msg);
    chomp $in;
    $in=~ s/\s+//;
    $ReadLine->addhistory($in) if $in =~ /\S/;
    return $in;
}

sub _samtools
{
    $ENV{HTSLIB_DIR} ||
    ( can_load(modules => {'Alien::SamTools' => undef, 'File::ShareDir' => undef}) &&
      File::ShareDir::dist_dir('Alien-SamTools'));
}
