/*
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
*/

#ifdef PERL_CAPI
#define WIN32IO_IS_STDIO
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef FCGI
 #include <fcgi_stdio.h>
#else
 #ifdef USE_SFIO
  #include <config.h>
 #else
  #include <stdio.h>
 #endif
 #include <perlio.h>
#endif

#ifndef Newx
#  define Newx(v,n,t) New(0,v,n,t)
#endif

#ifndef Newxz
#  define Newxz(v,n,t) Newz(0,v,n,t)
#endif

#include <unistd.h>
#include <math.h>
#include <string.h>
#include "kseq.h"
#include "hts.h"
#include "sam.h"
#include "khash.h"
#include "faidx.h"
#include "tbx.h"
#include "bgzf.h"
#include "vcf.h"
#include "synced_bcf_reader.h"

/* stolen from bam_aux.c */
#define BAM_MAX_REGION 1<<29

typedef htsFile*        Bio__DB__HTSfile;
typedef bam_hdr_t*      Bio__DB__HTS__Header;
typedef bam1_t*         Bio__DB__HTS__Alignment;
typedef hts_idx_t*      Bio__DB__HTS__Index;
typedef faidx_t*        Bio__DB__HTS__Fai;
typedef bam_pileup1_t*  Bio__DB__HTS__Pileup;
typedef tbx_t*          Bio__DB__HTS__Tabix;
typedef hts_itr_t*      Bio__DB__HTS__Tabix__Iterator;
typedef bcf_srs_t*      Bio__DB__HTS__VCF;
typedef bcf_hdr_t*      Bio__DB__HTS__VCF__Header;
typedef bcf1_t*         Bio__DB__HTS__VCF__Row;


typedef struct {
  SV* callback;
  SV* data;
} fetch_callback_data;
typedef fetch_callback_data *fetch_callback_dataptr;
typedef struct {
  int    start;
  int    end;
  double width;
  int    reads;
  int*   bin;
} coverage_graph;
typedef coverage_graph *coverage_graph_ptr;

static int MaxPileupCnt=8000;

void XS_pack_charPtrPtr( SV * arg, char ** array, int count) {
  int i;
  AV * avref;
  avref = (AV*)sv_2mortal((SV*)newAV());
  for (i=0; i<count; i++) {
    av_push(avref, newSVpv(array[i], strlen(array[i])));
  }
  SvSetSV( arg, newRV((SV*)avref));
}

int hts_fetch_fun (void *data, bam1_t *b)
{
  dSP;
  int count;

  fetch_callback_dataptr fcp;
  SV* callback;
  SV* callbackdata;
  SV* alignment_obj;
  bam1_t *b2;

  fcp          = (fetch_callback_dataptr) data;
  callback     = fcp->callback;
  callbackdata = fcp->data;

  /* turn the bam1_t into an appropriate object */
  /* need to dup it here so that the C layer doesn't reuse the address under Perl */
  b2 = bam_dup1(b);

  alignment_obj = sv_setref_pv(newSV(sizeof(bam1_t)),"Bio::DB::HTS::Alignment",(void*) b2);

  /* set up subroutine stack for the call */
  ENTER;
  SAVETMPS;

  PUSHMARK(SP);
  XPUSHs(sv_2mortal(alignment_obj));
  XPUSHs(callbackdata);
  PUTBACK;

  /* execute the call */
  count = call_sv(callback,G_SCALAR|G_DISCARD);

  FREETMPS;
  LEAVE;

  return 1;
}

int invoke_pileup_callback_fun(uint32_t tid,
			       uint32_t pos,
			       int n,
			       const bam_pileup1_t *pl,
			       void *data) {
  dSP;
  int count,i;
  fetch_callback_dataptr fcp;
  SV*  callback;
  SV*  callbackdata;
  SV*  pileup_obj;
  SV* p;
  SV** pileups;
  AV*  pileup;

  fcp          = (fetch_callback_dataptr) data;
  callback     = fcp->callback;
  callbackdata = fcp->data;

  /* turn the bam_pileup1_t into the appropriate object */
  /* this causes a compiler warning -- ignore it */
  pileup = newAV();
  av_extend(pileup,n);
  for (i=0;i<n;i++) {
    p = newSV(sizeof(bam_pileup1_t));
    sv_setref_pv(p,"Bio::DB::HTS::Pileup",(void*) &pl[i]);
    av_push(pileup,p);
  }

  /* set up subroutine stack for the call */
  ENTER;
  SAVETMPS;

  PUSHMARK(SP);
  XPUSHs(sv_2mortal(newSViv(tid)));
  XPUSHs(sv_2mortal(newSViv(pos)));
  XPUSHs(sv_2mortal(newRV_noinc((SV*)pileup)));
  XPUSHs(callbackdata);
  PUTBACK;

  /* execute the call */
  count = call_sv(callback,G_SCALAR|G_DISCARD);

  FREETMPS;
  LEAVE;
}

/*
   Declarations to allow add_pileup_line to work
   Ported from samtoosl v1 setup.
*/

/* start pileup support copy from bam.h in samtools */
/* but pileup functions are offered as bam_plp_auto_f in htslib */

typedef int (*bam_pileup_f)(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);

typedef struct
{
  bam_plp_t iter;
  bam_pileup_f func;
  void *data;
} hts_plbuf_t;


hts_plbuf_t *hts_plbuf_init(bam_pileup_f func, void *data)
{
    hts_plbuf_t *buf;
    buf = calloc(1, sizeof(hts_plbuf_t));
    buf->iter = bam_plp_init(0, 0);
    buf->func = func;
    buf->data = data;
    return buf;
}

void hts_plbuf_destroy(hts_plbuf_t *buf)
{
    bam_plp_destroy(buf->iter);
    free(buf);
}

int hts_plbuf_push(const bam1_t *b, hts_plbuf_t *buf)
{
    int ret, n_plp, tid, pos;
    const bam_pileup1_t *plp;
    ret = bam_plp_push(buf->iter, b);
    if (ret < 0) return ret;
    while ((plp = bam_plp_next(buf->iter, &tid, &pos, &n_plp)) != 0)
        buf->func(tid, pos, n_plp, plp, buf->data);
    return 0;
}


/* end pileup support copy from bam.h in samtools */

/**
   pileup support functions
*/
int add_pileup_line (void *data, bam1_t *b)
{
  hts_plbuf_t *pileup = (hts_plbuf_t*) data;
  hts_plbuf_push(b,pileup);
  return 0;
}



int coverage_from_pileup_fun (uint32_t tid,
			      uint32_t pos,
			      int n,
			      const bam_pileup1_t *pl,
			      void *data) {
  coverage_graph_ptr  cgp;
  int                 bin;
  int                 i;
  int                 valid;

  cgp = (coverage_graph_ptr) data;
  cgp->reads += n;

  valid = 0;
  for (i=0;i<n;i++) {
    if (!pl[i].is_del && !pl[i].is_refskip)
        valid++;
  }

  if (pos >= cgp->start && pos <= cgp->end) {
    bin = (pos-cgp->start)/cgp->width;
    cgp->bin[bin] += valid;
  }

  return 0;
}


/**
   From bam_aux.c in samtools. Needed to allow pileup function to work.
*/
int bam_parse_region(bam_hdr_t *header, const char *str, int *ref_id, int *beg, int *end)
{
    const char *name_lim = hts_parse_reg(str, beg, end);
    if (name_lim) {
        char *name = malloc(name_lim - str + 1);
        memcpy(name, str, name_lim - str);
        name[name_lim - str] = '\0';
        *ref_id = bam_name2id(header, name);
        free(name);
    }
    else {
        // not parsable as a region, but possibly a sequence named "foo:a"
        *ref_id = bam_name2id(header, str);
        *beg = 0; *end = INT_MAX;
    }
    if (*ref_id == -1) return -1;
    return *beg <= *end? 0 : -1;
}

/**
   From bam.c in samtools - these are wrappers that can be used OK here.
*/
char *bam_format1(const bam_hdr_t *header, const bam1_t *b)
{
    kstring_t str;
    str.l = str.m = 0; str.s = NULL;
    sam_format1(header, b, &str);
    return str.s;
}


void bam_view1(const bam_hdr_t *header, const bam1_t *b)
{
        char *s = bam_format1(header, b);
        puts(s);
        free(s);
}


/**
   Get the file extension for a filename
*/
int get_index_fmt_from_extension(const char * filename)
{
  char * ext = strrchr( filename, '.' ) ;
  if( strcmp(ext, ".cram")==0 )
  {
    return HTS_FMT_CRAI ;
  }
  if( strcmp(ext, ".bam")==0 )
  {
    return HTS_FMT_BAI ; //could also be HTS_FMT_CSI
  }
  return -1 ;
}

/**
   fetch function
*/
int hts_fetch(htsFile *fp, const hts_idx_t *idx, int tid, int beg, int end, void *data, bam_plp_auto_f func)
{
    int ret;
    hts_itr_t *iter ;
    bam1_t *b ;

    iter = sam_itr_queryi(idx, tid, beg, end);
    b = bam_init1();

    while((ret = sam_itr_next(fp, iter, b)) >= 0)
    {
        func(data,b);
    }
    hts_itr_destroy(iter);
    bam_destroy1(b);
    return (ret == -1)? 0 : ret;
}




MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::Fai PREFIX=fai_

Bio::DB::HTS::Fai
fai_load(packname="Bio::DB::HTS::Fai", filename)
  char * packname
  char * filename
 PROTOTYPE: $$
 CODE:
    RETVAL = fai_load(filename);
 OUTPUT:
    RETVAL

void
fai_destroy(fai)
  Bio::DB::HTS::Fai fai
  PROTOTYPE: $
  CODE:
    fai_destroy(fai);

SV*
fai_fetch(fai,reg)
  Bio::DB::HTS::Fai    fai
    const char *reg
  PROTOTYPE: $$$
  PREINIT:
    char     *seq;
    int       len;
  CODE:
    seq = fai_fetch(fai,reg,&len);
    if (seq == NULL)
       XSRETURN_EMPTY;
    RETVAL = newSVpv(seq,len);
    free((void*)seq);
  OUTPUT:
    RETVAL


MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTSfile PREFIX=hts_

int
max_pileup_cnt(packname,...)
CODE:
	if (items > 1)
	   MaxPileupCnt = SvIV(ST(1));
	RETVAL = MaxPileupCnt;
OUTPUT:
        RETVAL


Bio::DB::HTSfile
hts_open(packname, filename, mode="r")
      char * packname
      char * filename
      char * mode
      PROTOTYPE: $$$
      CODE:
        RETVAL = hts_open(filename,mode);
      OUTPUT:
      RETVAL


void
hts_close(htsfile)
   Bio::DB::HTSfile   htsfile
PROTOTYPE: $
CODE:
   hts_close(htsfile);


int
hts_index_build(packname, filename)
   char *      packname
   const char * filename
  CODE:
     RETVAL = sam_index_build(filename,0); //generate BAI for BAM files
  OUTPUT:
     RETVAL



Bio::DB::HTS::Index
hts_index_load(packname, htsfile)
    char *      packname
    Bio::DB::HTSfile htsfile
    PROTOTYPE: $$
    CODE:
      RETVAL = sam_index_load(htsfile, htsfile->fn) ;
    OUTPUT:
      RETVAL

void
hts_index_close(indexfile)
           Bio::DB::HTS::Index indexfile
    PROTOTYPE: $$
    CODE:
      hts_idx_destroy(indexfile) ;



Bio::DB::HTS::Header
hts_header_read(htsfile)
    Bio::DB::HTSfile htsfile
    PROTOTYPE: $$
    PREINIT:
      bam_hdr_t *bh;
      int64_t result ;
    CODE:
      if( htsfile->format.format == bam ) //enum value from htsExactFormat from hts.h
      {
        result = bgzf_seek(htsfile->fp.bgzf,0,0) ;
      }
      bh = sam_hdr_read(htsfile);
      RETVAL = bh ;
    OUTPUT:
      RETVAL


int
hts_header_write(hts,header)
    Bio::DB::HTSfile     hts
    Bio::DB::HTS::Header header
    PROTOTYPE: $$
    CODE:
      RETVAL= sam_hdr_write(hts,header);
    OUTPUT:
      RETVAL


Bio::DB::HTS::Alignment
hts_read1(htsfile,header)
    Bio::DB::HTSfile        htsfile
    Bio::DB::HTS::Header    header
  PROTOTYPE: $$
  PREINIT:
    bam1_t *alignment;
    CODE:
       alignment = bam_init1();
       if (sam_read1(htsfile,header,alignment) >= 0) {
         RETVAL = alignment ;
       }
       else
         XSRETURN_EMPTY;
    OUTPUT:
       RETVAL


MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::Alignment PREFIX=bama_

Bio::DB::HTS::Alignment
bama_new(package="Bio::DB::HTS::Alignment")
   char * package
   PROTOTYPE: $
   CODE:
      RETVAL = bam_init1();
   OUTPUT:
      RETVAL

void
bama_DESTROY(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
    bam_destroy1(b);

int
bama_tid(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.tid = SvIV(ST(1));
    RETVAL=b->core.tid;
OUTPUT:
    RETVAL

int
bama_pos(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.pos = SvIV(ST(1));
    RETVAL=b->core.pos;
OUTPUT:
    RETVAL

int
bama_calend(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
   RETVAL=bam_endpos(b);
OUTPUT:
   RETVAL

int
bama_cigar2qlen(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
   RETVAL=bam_cigar2qlen(b->core.n_cigar,bam_get_cigar(b));
OUTPUT:
   RETVAL

int
bama_qual(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.qual = SvIV(ST(1));
    RETVAL=b->core.qual;
OUTPUT:
    RETVAL

int
bama_flag(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.flag = SvIV(ST(1));
    RETVAL=b->core.flag;
OUTPUT:
    RETVAL

int
bama_n_cigar(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
  if (items > 1)
    b->core.n_cigar = SvIV(ST(1));
    RETVAL=b->core.n_cigar;
OUTPUT:
    RETVAL

int
bama_l_qseq(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.l_qseq = SvIV(ST(1));
    RETVAL=b->core.l_qseq;
OUTPUT:
    RETVAL

SV*
bama_qseq(b)
Bio::DB::HTS::Alignment b
PROTOTYPE: $
PREINIT:
    char* seq;
    int   i;
CODE:
    seq = Newxz(seq,b->core.l_qseq+1,char);
    for (i=0;i<b->core.l_qseq;i++) {
      seq[i]=seq_nt16_str[bam_seqi(bam_get_seq(b),i)];
    }
    RETVAL = newSVpv(seq,b->core.l_qseq);
    Safefree(seq);
OUTPUT:
    RETVAL

SV*
bama__qscore(b)
Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
    RETVAL = newSVpv(bam_get_qual(b),b->core.l_qseq);
OUTPUT:
    RETVAL

int
bama_mtid(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.mtid = SvIV(ST(1));
    RETVAL=b->core.mtid;
OUTPUT:
    RETVAL

int
bama_mpos(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.mpos = SvIV(ST(1));
    RETVAL=b->core.mpos;
OUTPUT:
    RETVAL

int
bama_isize(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.isize = SvIV(ST(1));
    RETVAL=b->core.isize;
OUTPUT:
    RETVAL

int
bama_l_aux(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    RETVAL=SvIV(newSViv(bam_get_l_aux(b)));
OUTPUT:
    RETVAL

char*
bama_aux(b)
   Bio::DB::HTS::Alignment b
PREINIT:
   uint8_t *s;
   uint8_t type, key[2];
   char    str[8192];
CODE:
   s = bam_get_aux(b);
   str[0] = '\0';

   int left  = sizeof(str) - strlen(str);
   while (left > 0 && (s < b->data + b->l_data)) {
        char* d   = str+strlen(str);

	key[0] = s[0];
	key[1] = s[1];
 	left -= snprintf(d, left, "%c%c:", key[0], key[1]);

	d    += 3;
	s    += 2;
	type = *s++;

	if (left <= 0) continue;

	if (type == 'A')      { left -= snprintf(d, left, "A:%c", *s);           s++; }
	else if (type == 'C') { left -= snprintf(d, left, "i:%u", *s);           s++; }
	else if (type == 'c') { left -= snprintf(d, left, "i:%d", *s);           s++; }
	else if (type == 'S') { left -= snprintf(d, left, "i:%u", *(uint16_t*)s);s += 2; }
	else if (type == 's') { left -= snprintf(d, left, "i:%d", *(int16_t*)s); s += 2; }
	else if (type == 'I') { left -= snprintf(d, left, "i:%u", *(uint32_t*)s);s += 4; }
	else if (type == 'i') { left -= snprintf(d, left, "i:%d", *(int32_t*)s); s += 4; }
	else if (type == 'f') { left -= snprintf(d, left, "f:%g", *(float*)s);   s += 4; }
	else if (type == 'd') { left -= snprintf(d, left, "d:%lg", *(double*)s); s += 8; }
	else if (type == 'Z' || type == 'H') { left -= snprintf(d, left, "%c:", type);
	                                       strncat(d,s,left);
					       while (*s++) {}
					       left = sizeof(str) - strlen(str);
	                                     }
	if (left <= 0) continue;
	strncat(d,"\t",left);
	left--;
   }
   str[strlen(str)-1] = '\0';
   RETVAL = str;
OUTPUT:
   RETVAL

SV*
bama_aux_get(b,tag)
   Bio::DB::HTS::Alignment b
   char*               tag
PROTOTYPE: $$
PREINIT:
   int           type;
   uint8_t       *s;
CODE:
   s    = bam_aux_get(b,tag);
   if (s==0)
      XSRETURN_EMPTY;
   type = *s++;
   switch (type) {
   case 'c':
     RETVAL = newSViv((int32_t)*(int8_t*)s);
     break;
   case 'C':
     RETVAL = newSViv((int32_t)*(uint8_t*)s);
     break;
   case 's':
     RETVAL = newSViv((int32_t)*(int16_t*)s);
     break;
   case 'S':
     RETVAL = newSViv((int32_t)*(uint16_t*)s);
     break;
   case 'i':
     RETVAL = newSViv(*(int32_t*)s);
     break;
   case 'I':
     RETVAL = newSViv((int32_t)*(uint32_t*)s);
     break;
   case 'f':
     RETVAL = newSVnv(*(float*)s);
     break;
   case 'Z':
   case 'H':
     RETVAL = newSVpv((char*)s,0);
     break;
   case 'A':
     RETVAL = newSVpv((char*)s,1);
     break;
   default:
     XSRETURN_EMPTY;
   }
OUTPUT:
   RETVAL

void
bama_aux_keys(b)
Bio::DB::HTS::Alignment b
PROTOTYPE: $
PREINIT:
   uint8_t *s;
   uint8_t type;
PPCODE:
   {
     s = bam_get_aux(b);  /* s is a khash macro */
     while (s < b->data + b->l_data) {
       XPUSHs(sv_2mortal(newSVpv(s,2)));
       s   += 2;
       type = *s++;
       if      (type == 'A') { ++s; }
       else if (type == 'C') { ++s; }
       else if (type == 'c') { ++s; }
       else if (type == 'S') { s += 2; }
       else if (type == 's') { s += 2; }
       else if (type == 'I') { s += 4; }
       else if (type == 'i') { s += 4; }
       else if (type == 'f') { s += 4; }
       else if (type == 'Z' || type == 'H') { while (*s) ++(s); ++(s); }
     }
   }

SV*
bama_data(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
PREINIT:
    STRLEN  len;
CODE:
    if (items > 1) {
      b->data     = SvPV(ST(1),len);
      b->l_data = len;
    }
    RETVAL=newSVpv(b->data,b->l_data);
OUTPUT:
    RETVAL

int
bama_data_len(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->l_data = SvIV(ST(1));
    RETVAL=b->l_data;
OUTPUT:
    RETVAL

int
bama_m_data(b,...)
    Bio::DB::HTS::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1) {
      b->m_data = SvIV(ST(1));
    }
    RETVAL=b->m_data;
OUTPUT:
    RETVAL

SV*
bama_qname(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=newSVpv(bam_get_qname(b),0);
OUTPUT:
    RETVAL

int
bama_paired(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FPAIRED) != 0;
OUTPUT:
  RETVAL

int
bama_proper_pair(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FPROPER_PAIR) != 0;
OUTPUT:
  RETVAL

int
bama_unmapped(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FUNMAP) != 0;
OUTPUT:
  RETVAL

int
bama_munmapped(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FMUNMAP) != 0;
OUTPUT:
  RETVAL

int
bama_reversed(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
  RETVAL=bam_is_rev(b);
OUTPUT:
  RETVAL

int
bama_mreversed(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
CODE:
  RETVAL=bam_is_mrev(b);
OUTPUT:
  RETVAL

SV*
bama_cigar(b)
  Bio::DB::HTS::Alignment b
PROTOTYPE: $
PREINIT:
    int        i;
    uint32_t  *c;
    AV        *avref;
CODE:
    avref = (AV*) sv_2mortal((SV*)newAV());
    c     = bam_get_cigar(b);
    for (i=0;i<b->core.n_cigar;i++)
      av_push(avref, newSViv(c[i]));
    RETVAL = (SV*) newRV((SV*)avref);
OUTPUT:
  RETVAL

MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::Header PREFIX=bam_

Bio::DB::HTS::Header
bam_new(packname=Bio::DB::HTS::Header)
PROTOTYPE: $
CODE:
    RETVAL = bam_hdr_init();
OUTPUT:
    RETVAL

int
bam_n_targets(bamh)
  Bio::DB::HTS::Header bamh
  PROTOTYPE: $
  CODE:
    RETVAL = bamh->n_targets;
  OUTPUT:
    RETVAL

SV*
bam_target_name(bamh)
  Bio::DB::HTS::Header bamh
  PROTOTYPE: $
  PREINIT:
    int i;
    AV * avref;
  CODE:
    avref = (AV*) sv_2mortal((SV*)newAV());
    for (i=0;i<bamh->n_targets;i++)
      av_push(avref, newSVpv(bamh->target_name[i],0));
    RETVAL = (SV*) newRV((SV*)avref);
  OUTPUT:
    RETVAL

SV*
bam_target_len(bamh)
    Bio::DB::HTS::Header bamh
  PROTOTYPE: $
  PREINIT:
    int i;
    AV * avref;
  CODE:
    avref = (AV*) sv_2mortal((SV*)newAV());
    for (i=0;i<bamh->n_targets;i++)
       av_push(avref, newSViv(bamh->target_len[i]));
    RETVAL = (SV*) newRV((SV*)avref);
  OUTPUT:
    RETVAL

SV*
bam_text(bamh, ...)
  Bio::DB::HTS::Header bamh
  PREINIT:
    char   *newtext;
    STRLEN n;
  CODE:
    /* in case text is not null terminated, we copy it */
    RETVAL = newSVpv(bamh->text,bamh->l_text);
    if (items > 1) {
      newtext = (char*) SvPV(ST(1),n);
      strcpy(bamh->text,newtext);
      bamh->l_text = n;
    }
  OUTPUT:
    RETVAL


void
bam_parse_region(bamh,region)
    Bio::DB::HTS::Header bamh
    char*            region
    PROTOTYPE: $
    PREINIT:
       int seqid,start,end;
    PPCODE:
    {
      bam_parse_region(bamh,
		       region,
		       &seqid,
		       &start,
		       &end);
      if (seqid < 0)
	XSRETURN_EMPTY;
      else {
	EXTEND(sp,3);
	PUSHs(sv_2mortal(newSViv(seqid)));
	PUSHs(sv_2mortal(newSViv(start)));
	PUSHs(sv_2mortal(newSViv(end)));
      }
    }

void
bam_view1(bamh,alignment)
     Bio::DB::HTS::Header     bamh
     Bio::DB::HTS::Alignment  alignment
     PROTOTYPE: $$
     CODE:
       bam_view1(bamh,alignment);


void
bam_DESTROY(bamh)
  Bio::DB::HTS::Header bamh
  PROTOTYPE: $
  CODE:
    bam_hdr_destroy(bamh);



MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::Index PREFIX=bami_

int
bami_fetch(bai,hfp,ref,start,end,callback,callbackdata=&PL_sv_undef)
  Bio::DB::HTS::Index bai
  Bio::DB::HTSfile    hfp
  int   ref
  int   start
  int   end
  CV*   callback
  SV*   callbackdata
PREINIT:
  fetch_callback_data fcd;
CODE:
  {
    fcd.callback = (SV*) callback;
    fcd.data     = callbackdata;
    RETVAL = hts_fetch(hfp,bai,ref,start,end,&fcd,hts_fetch_fun);
  }
OUTPUT:
    RETVAL


void
bami_pileup(bai,hfp,ref,start,end,callback,callbackdata=&PL_sv_undef)
  Bio::DB::HTS::Index bai
  Bio::DB::HTSfile    hfp
  int   ref
  int   start
  int   end
  CV*   callback
  SV*   callbackdata
PREINIT:
  fetch_callback_data fcd;
  hts_plbuf_t        *pileup;
CODE:
  fcd.callback = (SV*) callback;
  fcd.data     = callbackdata;
  pileup       = hts_plbuf_init(invoke_pileup_callback_fun,(void*)&fcd);
  bam_plp_set_maxcnt(pileup->iter,MaxPileupCnt);
  hts_fetch(hfp,bai,ref,start,end,(void*)pileup,add_pileup_line);
  hts_plbuf_push(NULL,pileup);
  hts_plbuf_destroy(pileup);

AV*
bami_coverage(bai,hfp,ref,start,end,bins=0,maxcnt=8000)
    Bio::DB::HTS::Index bai
    Bio::DB::HTSfile    hfp
    int             ref
    int             start
    int             end
    int             bins
    int             maxcnt
PREINIT:
    coverage_graph  cg;
    hts_plbuf_t    *pileup;
    AV*             array;
    SV*             cov;
    int             i;
    bam_hdr_t      *bh;
CODE:
  {
      /* TODO:can we do away with this check by a move to CSI as the standard for BAM indices */
      if (end >= BAM_MAX_REGION)
      {
        if( hfp->format.format == bam ) //enum value from htsExactFormat from hts.h
        {
          bgzf_seek(hfp->fp.bgzf,0,0);
          bh = sam_hdr_read(hfp);
          end = bh->target_len[ref];
          bam_hdr_destroy(bh);
        }
      }
      if ((bins==0) || (bins > (end-start)))
         bins = end-start;

      /* coverage graph used to communicate to our callback
	  the region we are sampling */
      cg.start = start;
      cg.end   = end;
      cg.reads = 0;
      cg.width = ((double)(end-start))/bins;
      Newxz(cg.bin,bins+1,int);

      /* accumulate coverage into the coverage graph */
      pileup   = hts_plbuf_init(coverage_from_pileup_fun,(void*)&cg);
      if (items >= 7)
            bam_plp_set_maxcnt(pileup->iter,maxcnt);
      else
            bam_plp_set_maxcnt(pileup->iter,MaxPileupCnt);
      hts_fetch(hfp,bai,ref,start,end,(void*)pileup,add_pileup_line);
      hts_plbuf_push(NULL,pileup);
      hts_plbuf_destroy(pileup);

      /* now normalize to coverage/bp and convert into an array */
      array = newAV();
      av_extend(array,bins);
      for  (i=0;i<bins;i++)
           av_store(array,i,newSVnv(((float)cg.bin[i])/cg.width));
      Safefree(cg.bin);
      RETVAL = array;
      sv_2mortal((SV*)RETVAL);  /* this fixes a documented bug in perl typemap */
  }
OUTPUT:
    RETVAL


void
bami_close(hts_idx)
  Bio::DB::HTS::Index hts_idx
  CODE:
    hts_idx_destroy(hts_idx) ;



MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::Pileup PREFIX=pl_

int
pl_qpos(pl)
  Bio::DB::HTS::Pileup pl
  CODE:
    RETVAL = pl->qpos;
  OUTPUT:
    RETVAL

int
pl_pos(pl)
  Bio::DB::HTS::Pileup pl
  CODE:
    RETVAL = pl->qpos+1;
  OUTPUT:
    RETVAL

int
pl_indel(pl)
  Bio::DB::HTS::Pileup pl
  CODE:
    RETVAL = pl->indel;
  OUTPUT:
    RETVAL

int
pl_level(pl)
  Bio::DB::HTS::Pileup pl
  CODE:
    RETVAL = pl->level;
  OUTPUT:
    RETVAL

int
pl_is_del(pl)
  Bio::DB::HTS::Pileup pl
  CODE:
    RETVAL = pl->is_del;
  OUTPUT:
    RETVAL

int
pl_is_refskip(pl)
  Bio::DB::HTS::Pileup pl
  CODE:
    RETVAL = pl->is_refskip;
  OUTPUT:
    RETVAL

int
pl_is_head(pl)
  Bio::DB::HTS::Pileup pl
  CODE:
    RETVAL = pl->is_head;
  OUTPUT:
    RETVAL

int
pl_is_tail(pl)
  Bio::DB::HTS::Pileup pl
  CODE:
    RETVAL = pl->is_tail;
  OUTPUT:
    RETVAL

Bio::DB::HTS::Alignment
pl_b(pl)
  Bio::DB::HTS::Pileup pl
  CODE:
    RETVAL = bam_dup1(pl->b);
  OUTPUT:
     RETVAL

Bio::DB::HTS::Alignment
pl_alignment(pl)
  Bio::DB::HTS::Pileup pl
  CODE:
    RETVAL = bam_dup1(pl->b);
  OUTPUT:
     RETVAL



MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::Tabix PREFIX = tabix_

Bio::DB::HTS::Tabix
tabix_tbx_open(fname)
    char *fname
  CODE:
    RETVAL = tbx_index_load(fname);
  OUTPUT:
    RETVAL


void
tabix_tbx_close(t)
    Bio::DB::HTS::Tabix t
  CODE:
    tbx_destroy(t);
  OUTPUT:

Bio::DB::HTS::Tabix::Iterator
tabix_tbx_query(t, region)
    Bio::DB::HTS::Tabix t
    char *region
  CODE:
    RETVAL = tbx_itr_querys(t, region);
  OUTPUT:
    RETVAL



SV*
tabix_tbx_header(fp, tabix)
    Bio::DB::HTSfile fp
    Bio::DB::HTS::Tabix tabix
  PREINIT:
    int num_header_lines = 0;
    AV *av_ref;
    kstring_t str = {0,0,0};
  CODE:
    av_ref = newAV();
    while ( hts_getline(fp, KS_SEP_LINE, &str) >= 0 ) {
        if ( ! str.l ) break; //no lines left so we are done
        if ( str.s[0] != tabix->conf.meta_char ) break;

        //the line begins with a # so add it to the array
        ++num_header_lines;
        av_push(av_ref, newSVpv(str.s, str.l));
    }

    if ( ! num_header_lines )
        XSRETURN_EMPTY;

    RETVAL = newRV_noinc((SV*) av_ref);
  OUTPUT:
    RETVAL

SV*
tabix_tbx_seqnames(t)
    Bio::DB::HTS::Tabix t
  PREINIT:
    const char **names;
    int i, num_seqs;
    AV *av_ref;
  CODE:
    names = tbx_seqnames(t, &num_seqs); //call actual tabix method

    //blast all the values onto a perl array
    av_ref = newAV();
    for (i = 0; i < num_seqs; ++i) {
        SV *sv_ref = newSVpv(names[i], 0);
        av_push(av_ref, sv_ref);
    }

    free(names);

    //return a reference to our array
    RETVAL = newRV_noinc((SV*)av_ref);
  OUTPUT:
    RETVAL

MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::Tabix::Iterator PREFIX = tabix_

SV*
tabix_tbx_iter_next(iter, fp, t)
    Bio::DB::HTS::Tabix::Iterator iter
    Bio::DB::HTSfile fp
    Bio::DB::HTS::Tabix t
  PREINIT:
    kstring_t str = {0,0,0};
  CODE:
    if (tbx_itr_next(fp, t, iter, &str) < 0)
        XSRETURN_EMPTY;

    RETVAL = newSVpv(str.s, str.l);
  OUTPUT:
    RETVAL

void
tabix_tbx_iter_free(iter)
	Bio::DB::HTS::Tabix::Iterator iter
  CODE:
	tbx_itr_destroy(iter);


MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::VCF PREFIX = vcf_

Bio::DB::HTS::VCF
vcf_bcf_sr_open(filename)
    char* filename
    PREINIT:
        bcf_srs_t* sr = bcf_sr_init();
    CODE:
        bcf_sr_add_reader(sr, filename);
        RETVAL = sr;
    OUTPUT:
        RETVAL


Bio::DB::HTS::VCF::Header
vcf_bcf_header(vcf)
    Bio::DB::HTS::VCF vcf
    PREINIT:
        bcf_hdr_t* h;
    CODE:
        h = vcf->readers[0].header;
        RETVAL = h;
    OUTPUT:
        RETVAL


Bio::DB::HTS::VCF::Row
vcf_bcf_next(vcf)
    Bio::DB::HTS::VCF vcf
    PREINIT:
        bcf1_t* line;
    CODE:
        if ( bcf_sr_next_line(vcf) ) {
            line = bcf_sr_get_line(vcf, 0); //0 being the first and only reader
            RETVAL = line;
        }
        else {
            XSRETURN_EMPTY;
        }
    OUTPUT:
        RETVAL


SV*
vcf_bcf_num_variants(vcf)
    Bio::DB::HTS::VCF vcf
    PREINIT:
        int n_records = 0;
    CODE:
        //loop through all the lines but don't do anything with them
        while ( bcf_sr_next_line(vcf) ) {
            ++n_records;
        }

        RETVAL = newSViv(n_records);
    OUTPUT:
        RETVAL


void
vcf_bcf_sr_close(vcf)
    Bio::DB::HTS::VCF vcf
    CODE:
        bcf_sr_destroy(vcf);
  OUTPUT:
