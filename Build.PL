#!/usr/bin/perl

# Copyright [2015-2023] EMBL-European Bioinformatics Institute
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

use Module::Build;

my $class = Module::Build->subclass(
    class => 'Module::Build::HTS',
    );

my $build = $class->new(

    module_name        => 'Bio::DB::HTS',
    dist_version_from  => 'lib/Bio/DB/HTS.pm',
    dist_author        => 'Alessandro Vullo',
    dist_abstract      => 'Perl interface to HTS library for DNA sequencing',
    license            => 'Apache_2_0',

    extra_compiler_flags => [

        # must match DFLAGS in HTSlib Makefile
        '-D_IOLIB=2', '-D_FILE_OFFSET_BITS=64',
        # warnings not treated as errors
        '-Wno-error',
        # Don't care about unused results from function calls
        '-Wno-unused-result',  ],

    build_requires => { 'ExtUtils::CBuilder' => 0, },
    configure_requires => { 'Module::Build' => 0.38, },
    requires => { 'perl' => '5.008', 'Bio::SeqFeature::Lite' => 0 },
    meta_merge => {
        'resources' => {
            'repository' => 'https://github.com/Ensembl/Bio-DB-HTS',
            },
        },

    add_to_cleanup => [
        't/data/*ai',
        't/data/*.gzi',
        't/data/write_test_07.cram',
        't/data/write_test_08.bam',
        ],
);

$build->find_hts;
$build->set_include_and_compiler_flags;
$build->create_build_script;

exit 0;


package Module::Build::HTS;

use Module::Load::Conditional qw(can_load);
use base 'Module::Build';

sub find_hts {
    my ($self) = @_;

    # If either of these are set, we expect to find the htslib files there:
    # (They're explicitly set by the user, so we shouldn't fall back to
    # finding another copy somewhere else.)
    my $incdir = $self->args('htslib-includedir');
    my $libdir = $self->args('htslib-libdir');
    if ($incdir && $libdir) {
        return 1 if $self->find_hts_in_split_install_dirs($incdir, $libdir);
        $self->die_hts_not_found(
            "--htslib-includedir '$incdir' or --htslib-libdir '$libdir' command line parameters do not contain expected files\n"
        );
    }
    elsif (my $dir = $self->args('htslib')) {
        return 1 if $self->find_hts_in_build_dir($dir);
        return 1 if $self->find_hts_in_install_dir($dir);
        $self->die_hts_not_found(
            "--htslib '$dir' command line parameter does not contain expected files\n"
        );
    }
    elsif ($dir = $ENV{'HTSLIB_DIR'}) {
        return 1 if $self->find_hts_in_build_dir($dir);
        return 1 if $self->find_hts_in_install_dir($dir);
        $self->die_hts_not_found(
            "HTSLIB_DIR=$ENV{HTSLIB_DIR} environment variable does not contain expected files\n"
        );
    }

    # Search through remaining possible (but not fatal) locations:
    my $found = 0;
    foreach my $dir (
        $self->prefix,
        from_Alien(),
        scalar `pkg-config --variable=libdir htslib 2>/dev/null`,
        qw{ /usr /usr/local /usr/share /opt/local },
    ) {
        if ($dir and $self->find_hts_in_install_dir($dir)) {
            $found = 1;
            last;
        }
    }
    return 1 if $found;

    # Try pkgconfig again but this time trust whatever it returns, without extra verification
    return 1 if $self->find_hts_with_pkgconfig();

    $self->die_hts_not_found();
}

sub set_include_and_compiler_flags {
    my ($self) = @_;

    my $hts_include = $self->config_data('hts_include');
    my $hts_lib     = $self->config_data('hts_lib');
    my $static      = $self->args('static');

    $self->include_dirs([$hts_include]) if $hts_include;

    my @linker_dirflags;
    if ($hts_lib) {
      push @linker_dirflags, "-L$hts_lib";
      push @linker_dirflags, "-Wl,-rpath,$hts_lib" unless $static;
    }
    $self->extra_linker_flags(@linker_dirflags, '-lhts', '-lpthread', '-lz');
}

sub hts_dev_files_exist {
    my ($lib, $include) = @_;
    return (-f "$lib/libhts.a" || -f "$lib/libhts.so") && -f "$include/htslib/hts.h";
}

# Look for the library and header in the location where htslib was compiled
sub find_hts_in_build_dir {
    my ($self, $root) = @_;

    chomp($root);
    $root =~ s{/$}{};
    $root =~ s{/(lib|include|include/htslib)$}{};

    if (hts_dev_files_exist($root, $root)) {
        $self->config_data('hts_lib'     => $root);
        $self->config_data('hts_include' => $root);
        return 1;
    }
    else {
        return 0;
    }
}


sub find_hts_in_install_dir {
    my ($self, $root) = @_;

    chomp($root);
    $root =~ s{/$}{};
    $root =~ s{/(lib|include|include/htslib)$}{};

    my $hts_lib     = "$root/lib";
    my $hts_include = "$root/include";
    if (hts_dev_files_exist($hts_lib, $hts_include)) {
        $self->config_data('hts_lib'     => $hts_lib);
        $self->config_data('hts_include' => $hts_include);
        return 1;
    }
    else {
        return 0;
    }
}

sub find_hts_in_split_install_dirs {
    my ($self, $hts_include, $hts_lib) = @_;

    chomp($hts_lib);
    chomp($hts_include);
    $hts_include =~ s{include/htslib$}{include};

    if (hts_dev_files_exist($hts_lib, $hts_include)) {
        $self->config_data('hts_lib'     => $hts_lib);
        $self->config_data('hts_include' => $hts_include);
        return 1;
    }
    else {
        return 0;
    }
}

sub find_hts_with_pkgconfig {
    my ($self) = @_;

    return 0 unless can_load(
        modules => { 'ExtUtils::PkgConfig' => undef }
    );
    my $pkg_name = 'htslib';
    return 0 unless ExtUtils::PkgConfig->exists($pkg_name);

    if (my $libs_only_L = ExtUtils::PkgConfig->libs_only_L($pkg_name)) {
        # For compatibility with other htslib search methods. Note that this
        # assumes there will be at most one each of -I/-L directories, which
        # is true as of htslib-1.5 but might change in the future.
        $libs_only_L =~ s{^-L}{};
        $self->config_data('hts_lib' => $libs_only_L);
    }
    if (my $cflags_only_I = ExtUtils::PkgConfig->cflags_only_I($pkg_name)) {
        # See above
        $cflags_only_I =~ s{^-I}{};
        $self->config_data('hts_include' => $cflags_only_I);
    }

    return 1;
}

sub die_hts_not_found {
    my ($self, $msg) = @_;

    $msg ||= '';
    die $msg, <<END;

This module requires HTSlib (http://htslib.org/)
Install it if you have not done so already.

This script will attempt to locate HTSlib by looking for htslib/hts.h
and libhts.a / libhts.so in:

  1. --htslib-includedir and --htslib-libdir command line arguments
  2. --htslib command line argument
  3. HTSLIB_DIR environment variable
  4. --prefix command line argument (which also sets installation location)
  5. Alien::HTSlib dependency resolver
  6. pkg-config (extra directories can be set in PKG_CONFIG_PATH environment variable)
  7. common library locations: /usr /usr/local, /usr/share, /opt/local

If none of the above succeeds but htslib is registered with pkg-config, the script
will try using pkg-config paths (via ExtUtils::PkgConfig) without checking if header
and library files exist.

END

}

sub from_Alien {
    can_load(
        modules => { 'Alien::HTSlib' => undef, 'File::ShareDir' => undef }
    ) && File::ShareDir::dist_dir('Alien-HTSlib');
}
