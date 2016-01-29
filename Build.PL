#!/usr/bin/perl

use strict;
use Module::Build;
use Module::Load::Conditional qw(can_load);

my $HeaderFile = "hts.h";
my $LibFile    = "libhts.a";
my $ReadLine;

my ( $hts_include, $hts_lib ) = find_hts();    # may exit with error here

my $class = Module::Build->subclass( code => <<EOF);


EOF

my $build = $class->new(
    module_name        => 'Bio::DB::HTS',
    dist_version_from  => 'lib/Bio/DB/HTS.pm',
    dist_author        => 'Rishi Nag',
    dist_abstract      => 'Perl interface to HTS library for DNA sequencing',
    license            => 'Apache_2_0',
    include_dirs       => [$hts_include],
    extra_linker_flags => [ "-L$hts_lib", '-lhts', '-lpthread', '-lz' ],

    extra_compiler_flags => [

        # must match DFLAGS in HTSlib Makefile
        '-D_IOLIB=2', '-D_FILE_OFFSET_BITS=64',

        # turn off warnings originatiNg in Perl's Newx* calls
        '-Wformat=0', ],

    build_requires => { 'ExtUtils::CBuilder' => 0, },
    configure_requires => { 'Module::Build' => 0.38, },
    requires => { 'perl' => '5.008', 'Bio::Root::Version' => '1.006009001', },
    #    create_makefile_pl => 'passthrough',
);

$build->create_build_script;

exit 0;

sub find_hts {
    my ( $hts_include, $hts_lib );

    if ( my $htslib_dir = _htstools() ) {
        $hts_include = $htslib_dir if -e "$htslib_dir/$HeaderFile";
        $hts_include = "$htslib_dir/include"
          if -e "$htslib_dir/include/$HeaderFile";
        $hts_include = "$htslib_dir/htslib" if -e "$htslib_dir/htslib/$HeaderFile";
        $hts_lib     = $htslib_dir          if -e "$htslib_dir/$LibFile";
        $hts_lib     = "$htslib_dir/lib"    if -e "$htslib_dir/lib/$LibFile";
    }

    my @search_path = qw(/ /usr /usr/share /usr/local);

    unless ($hts_include) {
        for my $p (@search_path) {
            $hts_include ||= "$p/include" if -e "$p/include/$HeaderFile";
        }
    }

    unless ($hts_lib) {
        for my $p (@search_path) {
            $hts_lib ||= "$p/lib" if -e "$p/lib/$LibFile";
        }
    }

    unless ( $hts_include && $hts_lib ) {
        print STDOUT "This module requires htslib (http://htslib.org).\n";
        my $prompt =
"Please enter the location of the $HeaderFile and compiled $LibFile files: ";
        my $found;
        while ( !$found ) {
            my $path = prompt($prompt);
            print STDOUT "\n";
            last unless $path;
            $hts_include = $path           if -e "$path/$HeaderFile";
            $hts_include = "$path/include" if -e "$path/include/$HeaderFile";
            $hts_include = "$path/htslib"  if -e "$path/htslib/$HeaderFile";
            $hts_lib     = $path           if -e "$path/$LibFile";
            $hts_lib     = "$path/lib"     if -e "$path/lib/$LibFile";
            $found = $hts_include && $hts_lib;

            unless ($found) {
                print STDOUT "That didn't seem to be right.\n";
                $prompt = "Try again, or hit <enter> to cancel: ";
            }
        }
    }

    unless ( $hts_include && $hts_lib ) {
        die <<END;
Can\'t find $LibFile and/or $HeaderFile!

If you haven\'t done so already, please compile htslib from
http://htslib.org and set the HTSLIB_DIR
environment variable to point to a  distribution directory
containing the compiled $LibFile and $HeaderFile files.
END
    }

    print STDOUT "Found $hts_include/$HeaderFile and $hts_lib/$LibFile\n";
    return ( $hts_include, $hts_lib );
} ## end sub find_hts

sub prompt {
    my $msg = shift;

    unless ( defined $ReadLine ) {
        eval "require Term::ReadLine";
        $ReadLine = Term::ReadLine->can('new') || 0;
        $ReadLine &&= Term::ReadLine->new( \*STDOUT );
        eval { readline::rl_set( 'TcshCompleteMode', 'On' ) };
    }

    unless ($ReadLine) {
        print STDOUT $msg;
        chomp( my $in = <> );
        return $in;
    }

    my $in = $ReadLine->readline($msg);
    chomp $in;
    $in =~ s/\s+//;
    $ReadLine->addhistory($in) if $in =~ /\S/;
    return $in;
}

sub _htstools {
    $ENV{HTSLIB_DIR} || (
        can_load(
            modules => { 'Alien::HTSlib' => undef, 'File::ShareDir' => undef }
        ) &&
        File::ShareDir::dist_dir('Alien-HTSlib') );
}
