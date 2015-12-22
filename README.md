This is a Perl interface to the HTS Library. See http://htslib.org

## ONE-STEP INSTALLATION

In the root directory of this distribution you will find the script
INSTALL.pl. Running this will download the latest versions of this
module and HTSlib into a temporary directory, compile them, test and
install Bio-HTSTools. Note that HTSlib will need to be installed
on your system, or the LD_LIBRARY_PATH (DYLD_LIBRARY_PATH on Mac OSX)
includes the location of libhts.a

Simply run:
  perl INSTALL.pl

## MULTI-STEP INSTALLATION

The more traditional install requires you to separately download,
unpack and compile HTSlib.

Then set the environment variable HTSLIB_DIR to point to this directory.

You will also need to install Bio::Perl from CPAN.

Now run:
```
  perl Build.PL
  ./Build
  ./Build test
  (sudo) ./Build install
```

## TROUBLESHOOTING:

If you encounter problems during compiling, you may need to edit
Build.PL so that extra_compiler_flags matches the CFLAGS and DFLAGS
settings in the HTSlib Makefile.  Here are some common problems:

1. When building this module, you get an error like the following:
relocation R_X86_64_32 against `a local symbol' can not be used when
making a shared object; recompile with -fPIC

To fix this, edit the Makefile in the Samtools distribution by adding
"-fPIC" to the CFLAGS line. While you're at it, you may also wish to
get rid of a bunch of unused variable warnings that appears under
recent versions of gcc. The modified CFLAGS will look like this

  CFLAGS= -g -Wall -Wno-unused -Wno-unused-result -O2 -fPIC #-m64 #-arch ppc

Then do "make clean; make" in the Samtools directory to recompile the
library. After this you should be able to build this module without
errors.

2. When building this module, you get an error about a missing math
library.

To fix this, follow the recipe in (1) except add -m64 to CFLAGS so it
looks like this:

  CFLAGS=	-g -Wall -O2 -fPIC #-m64 #-arch ppc

## TESTING AND CONTRIBUTING:

You can obtain the most recent development version of this module via
its GitHub site at https://github.com/Ensembl/Bio-HTS. Please
feel free to submit bug reports, patches, etc. via GitHub.

## AUTHOR:

Rishi Nag <rishi@ebi.ac.uk>

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
