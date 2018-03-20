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
use warnings;

use Module::Build;

my $class = Module::Build->subclass(
    class => 'Module::Build::HTS',
    );

my $build = $class->new(

    module_name        => 'Bio::DB::HTS',
    dist_version_from  => 'lib/Bio/DB/HTS.pm',
    dist_author        => 'Rishi Nag',
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
    configure_requires => { 'Module::Build' => 0.42, },
    requires => { 'perl' => '5.008', 'Bio::Root::Version' => '1.006001', },
    meta_merge => {
        'resources' => {
            'repository' => 'https://github.com/Ensembl/Bio-DB-HTS',
            },
        },
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
    if (my $dir = $self->args('htslib')) {
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

    $self->die_hts_not_found();
}

sub set_include_and_compiler_flags {
    my ($self) = @_;

    my $hts_include = $self->config_data('hts_include');
    my $hts_lib     = $self->config_data('hts_lib');
    my $static      = $self->args('static');
    $self->include_dirs([$hts_include]);
    if($static){
      $self->extra_linker_flags("-L$hts_lib", '-lhts', '-lpthread', '-lz');
    }else{
      $self->extra_linker_flags("-L$hts_lib", "-Wl,-rpath,$hts_lib", '-lhts', '-lpthread', '-lz');
    }
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

sub die_hts_not_found {
    my ($self, $msg) = @_;

    $msg ||= '';
    die $msg, <<END;

This module requires HTSlib (http://htslib.org/)
Install it if you have not done so already.

This script will attempt to locate HTSlib by looking for htslib/hts.h
and libhts.a / libhts.so in:

  1. --htslib command line argument
  2. HTSLIB_DIR environment variable
  3. --prefix command line argument (which also sets installation location)
  4. Alien::HTSlib dependency resolver
  5. pkg-config (extra directories can be set in PKG_CONFIG_PATH environment variable)
  6. common library locations: /usr /usr/local, /usr/share, /opt/local

END

}

sub from_Alien {
    can_load(
        modules => { 'Alien::HTSlib' => undef, 'File::ShareDir' => undef }
    ) && File::ShareDir::dist_dir('Alien-HTSlib');
}
