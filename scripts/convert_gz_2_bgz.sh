#!/bin/bash

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

script_name="convert_gz_2_bgz.sh"


if [ -z "$1" ]; then
    echo "Usage: $script_name gzipped_file|directory_with_gz_files [bgzip_program]"
    exit 1
fi


if [ -z "$2" ]; then
    bgzip="bgzip"
else
    bgzip=$(cd "$(dirname "$2")" && pwd)/$(basename "$2") 
    echo "Using user specified bgzip:$bgzip"
fi

#check for presence of programs as required
if !(hash gzip 2>/dev/null;) then
    echo "$script_name: GZIP not found, exiting"
    exit 1
fi

if !(hash $bgzip 2>/dev/null;) then
    echo "$script_name: BGZIP not found, exiting"
    exit 1
fi

#process 
if [ -f "$1" ]; then
    dir="${1%/*}"
    file="${1##*/}"
    filename="${1%.*}"
    echo "Making bgzip file $filename from $file"
    cd $dir
    gzip -d $file
    $bgzip $filename
    cd -
elif [ -d "$1" ]; then
    cd $1
    find . -name "*.gz" -print0 | while read -d $'\0' file
    do
        echo $file
        filename="${file%.*}"
        gzip -d $file
        $bgzip $filename
    done
    cd -
else
    echo "$script_name: argument $1 not found"
    exit 1
fi

exit $?




