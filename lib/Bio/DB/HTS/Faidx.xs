/*
 * Copyright [2015-2018] EMBL-European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <stdio.h>
#include <string.h>

#include "htslib/faidx.h"
#ifndef Newx
#  define Newx(v,n,t) New(0,v,n,t)
#endif


// Code is written to use a blessed int pointer to this strut as an object
// You cannot use Data::Dumper to inspect the Faidx object. Sorry
typedef struct
{
  char* path;
  faidx_t* index;
} Faidx ;


SV* new(const char * classname, const char * path)
{
  Faidx   * faidx;
  SV      * obj;
  SV      * obj_ref;
  faidx_t * fai;

  Newx(faidx, 1, Faidx);

  fai = fai_load(path);
  faidx->path = savepv(path);
  faidx->index = fai;
  obj = newSViv((IV)faidx);
  obj_ref = newRV_noinc(obj);
  sv_bless(obj_ref, gv_stashpv(classname, GV_ADD));
  SvREADONLY_on(obj);

  return obj_ref;
}


void get_sequence(SV* obj, SV* location, SV** seq, int* seq_len)
{
  faidx_t *fai;
  char* char_seq;

  *seq = newSVpvn("",0);
  *seq_len = 0;

  fai = ((Faidx*)SvIV(SvRV(obj)))->index;
  //Fetch sequence
  char_seq = fai_fetch(fai, SvPV(location, PL_na), seq_len);

  //Push into a SV
  sv_catpv(*seq, char_seq);
  sv_2mortal(*seq);
  //Free the buffer created by faidx
  free(char_seq);
}

void get_sequence2(SV* obj, SV* seq_id, int start, int end, SV** seq, int* seq_len)
{
  faidx_t *fai;
  char* char_seq;

  *seq = newSVpvn("",0);
  *seq_len = 0;

  fai = ((Faidx*)SvIV(SvRV(obj)))->index;
  //Fetch sequence
  char_seq = faidx_fetch_seq(fai, SvPV(seq_id, PL_na), start, end, seq_len);

  //Push into a SV
  sv_catpv(*seq, char_seq);
  sv_2mortal(*seq);
  //Free the buffer created by faidx
  free(char_seq);
}


int has_sequence(SV* obj, SV* seq_id)
{
  int has_seq=-1 ;
  has_seq = faidx_has_seq(((Faidx*)SvIV(SvRV(obj)))->index, SvPV(seq_id, PL_na));
  return has_seq;
}


int length(SV* obj, char* seq_id)
{
    int length = 0 ;
    faidx_t *fai = ((Faidx*)SvIV(SvRV(obj)))->index ;
    length = faidx_seq_len(fai, seq_id) ;
    return length ;
}


void DESTROY(SV* obj)
{
  Faidx* faidx = (Faidx*)SvIV(SvRV(obj));
  Safefree(faidx->path);
  fai_destroy(faidx->index);
  Safefree(faidx);
}


MODULE = Bio::DB::HTS::Faidx    PACKAGE = Bio::DB::HTS::Faidx
PROTOTYPES: ENABLE


SV*
new(classname, path)
  char* classname
  char* path


void
get_sequence(obj, location, OUTLIST seq, OUTLIST length)
  SV* obj
  SV* location
  SV* seq
  int length = NO_INIT
  CODE:
     get_sequence(obj, location, &seq, &length) ;


void
get_sequence_no_length(obj, location, OUTLIST seq)
  SV* obj
  SV* location
  SV* seq
CODE:
  int seq_len=0 ;
  get_sequence(obj, location, &seq, &seq_len) ;


void
get_sequence2(obj, seq_id, start, end, OUTLIST seq, OUTLIST length)
  SV* obj
  SV* seq_id
  int start
  int end
  SV* seq
  int length = NO_INIT


void
get_sequence2_no_length(obj, seq_id, start, end, OUTLIST seq)
  SV* obj
  SV* seq_id
  int start
  int end
  SV* seq
CODE:
  int seq_len=0 ;
  get_sequence2(obj, seq_id, start, end, &seq, &seq_len) ;


int
has_sequence(obj, seq_id)
  SV* obj
  SV* seq_id


int
length(obj, seq_id)
  SV* obj
  char* seq_id


void
get_all_sequence_ids(obj)
   SV* obj
INIT:
   int num_seqs ;
   int i ;
    const char* faidx_name ;
PPCODE:
    num_seqs = 0 ;
    faidx_t *fai = ((Faidx*)SvIV(SvRV(obj)))->index ;
    num_seqs = faidx_nseq(fai) ;
    EXTEND(SP,num_seqs);
    for( i=0 ; i<num_seqs ; i++ )
    {
      faidx_name = faidx_iseq(fai,i) ;
      PUSHs(sv_2mortal(newSVpv(faidx_name,0))) ;
    }


void
DESTROY(obj)
  SV* obj
