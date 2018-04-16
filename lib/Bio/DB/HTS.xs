/*
Copyright [2015-2018] EMBL-European Bioinformatics Institute

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

#define TRACEME(x) do {						\
    if (SvTRUE(perl_get_sv("Bio::DB::HTS::ENABLE_DEBUG", TRUE)))	\
      { PerlIO_stdoutf (x); PerlIO_stdoutf ("\n"); }		\
  } while (0)

#ifndef Newx
#  define Newx(v,n,t) New(0,v,n,t)
#endif

#ifndef Newxz
#  define Newxz(v,n,t) Newz(0,v,n,t)
#endif

#include <unistd.h>
#include <math.h>
#include <string.h>
#include <zlib.h>

#include "htslib/kseq.h"
#include "htslib/hts.h"
#include "htslib/hfile.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/vcf_sweep.h"
#include "htslib/synced_bcf_reader.h"

#include "htslib/khash.h"
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

/* stolen from bam_aux.c */
#define BAM_MAX_REGION 1<<29

typedef htsFile*        Bio__DB__HTSfile;
typedef htsFile*        Bio__DB__HTS__VCFfile;
typedef bam_hdr_t*      Bio__DB__HTS__Header;
typedef bam1_t*         Bio__DB__HTS__Alignment;
typedef hts_idx_t*      Bio__DB__HTS__Index;
typedef faidx_t*        Bio__DB__HTS__Fai;
typedef bam_pileup1_t*  Bio__DB__HTS__Pileup;
typedef tbx_t*          Bio__DB__HTS__Tabix;
typedef hts_itr_t*      Bio__DB__HTS__Tabix__Iterator;
typedef vcfFile*        Bio__DB__HTS__VCFfile;
typedef hts_itr_t*      Bio__DB__HTS__VCF__Iterator;
typedef bcf_hdr_t*      Bio__DB__HTS__VCF__Header;
typedef bcf_hdr_t*      Bio__DB__HTS__VCF__HeaderPtr;
typedef bcf1_t*         Bio__DB__HTS__VCF__Row;
typedef bcf1_t*         Bio__DB__HTS__VCF__RowPtr;
KSEQ_INIT(gzFile, gzread)
typedef gzFile          Bio__DB__HTS__Kseq;
typedef kseq_t*         Bio__DB__HTS__Kseq__Iterator;
typedef kstream_t*      Bio__DB__HTS__Kseq__Kstream;
typedef kstring_t*      Bio__DB__HTS__Kseq__Kstring;
typedef bcf_sweep_t*    Bio__DB__HTS__VCF__Sweep;

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

static int invoke_sv_to_int_fun(SV *func, SV *arg)
{
  dSP;
  int count, ret;

  ENTER;
  SAVETMPS;

  PUSHMARK(SP);
  XPUSHs(arg);
  PUTBACK;

  count = call_sv(func, G_SCALAR);
  SPAGAIN;

  if (count != 1) return -1;

  ret = POPi;
  PUTBACK;
  FREETMPS;
  LEAVE;

  return ret;
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

  /* The underlying bam1_t will be bam_destroy1()ed by alignment_obj's
   * destructor, so we need to duplicate it here. We could create the Perl SV
   * alongside the C bam1_t (cf bami_coverage), but note that a new b & b_sv
   * would be needed for each iteration, as some callback functions will
   * expect distinct references to distinct alignment objects each time.
   */
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
  SV* p;
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

  return 0;
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


MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS

SV*
htslib_version(packname="Bio::DB::HTS")
    char * packname
  PROTOTYPE: $
  CODE:
    RETVAL = newSVpv(hts_version(), 0);
  OUTPUT:
    RETVAL

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
fai_DESTROY(fai)
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


int
hts_is_remote(packname, filename)
      char * packname
      char * filename
      PROTOTYPE: $$
      CODE:
        RETVAL = hisremote(filename) ;
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
hts_DESTROY(htsfile)
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
    PREINIT:
      SV *htsfile_sv = SvRV(ST(1));
      HV *assocfile = get_hv("Bio::DB::HTS::_associated_file", GV_ADD);
    CODE:
      RETVAL = sam_index_load(htsfile, htsfile->fn) ;

      /* For CRAM, it is important that this hts_idx_t is destroyed *before*
       * the associated htsFile. We hold on to a reference to the htsFile
       * (which we'll release in our destructor) to ensure it outlives us.
       */
      SvREFCNT_inc(htsfile_sv);
      hv_store(assocfile, (char *) &RETVAL, sizeof RETVAL, htsfile_sv, 0);
    OUTPUT:
      RETVAL

Bio::DB::HTS::Header
hts_header_read(htsfile)
    Bio::DB::HTSfile htsfile
    PROTOTYPE: $$
    PREINIT:
      bam_hdr_t *bh;
      int64_t result ;
      const htsFormat *format ;
    CODE:
      format = hts_get_format( htsfile ) ;
      if( format->format == bam  ) //enum value from htsExactFormat from hts.h
        result = bgzf_seek(htsfile->fp.bgzf,0,0) ;
      /*
       * https://github.com/Ensembl/Bio-DB-HTS/issues/54
       * must seek at beginning of file for sam as well
       */
      else if ( format->format == sam ) {
	/*
	 * Using hseek with htslib < 1.5 triggers segfault
	 * and couldn't find a valid alternative.
	 *
	 * WARNING: we're tied to buggy behaviour for htslib <= 1.3.1
	 */
	if ( strcmp(hts_version(), "1.5") >= 0 )
	  result = hseek(htsfile->fp.hfile, 0, SEEK_SET);
      }

      bh = sam_hdr_read(htsfile);
      RETVAL = bh ;
    OUTPUT:
      RETVAL


int
hts_header_write(htsfile,header, ...)
    Bio::DB::HTSfile     htsfile
    Bio::DB::HTS::Header header
    PREINIT:
      char *reference = "";
      const htsFormat *format ;
    PROTOTYPE: $$
    CODE:
      format = hts_get_format( htsfile ) ;
      if( format->format == cram )
      {
        if(items > 2)
        {
          reference = (char *)SvPV_nolen(ST(2));
          hts_set_fai_filename(htsfile, reference);
        }
        else
        {
          croak("Error: need reference sequence file for writing CRAM file '%s'", htsfile->fn);
        }
      }
      RETVAL= sam_hdr_write(htsfile,header);
      if (RETVAL != 0)
        croak("Error %d while creating file '%s'", RETVAL, htsfile->fn);
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
       else {
         bam_destroy1(alignment);
         XSRETURN_EMPTY;
       }
    OUTPUT:
       RETVAL

int
hts_write1(htsfile,header,align)
    Bio::DB::HTSfile            htsfile
    Bio::DB::HTS::Header        header
    Bio::DB::HTS::Alignment     align
    PROTOTYPE: $$
    CODE:
      RETVAL = sam_write1(htsfile,header,align);
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
    Newxz(seq,b->core.l_qseq+1,char);
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
    RETVAL = newSVpv((char *) bam_get_qual(b),b->core.l_qseq);
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
	                                       strncat(d, (char *) s, left);
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
       XPUSHs(sv_2mortal(newSVpv((char *) s, 2)));
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
      b->data     = (uint8_t *) SvPV(ST(1),len);
      b->l_data = len;
    }
    RETVAL=newSVpv((char *) b->data, b->l_data);
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
bami_coverage(bai,hfp,ref,start,end,bins=0,maxcnt=8000,filter=NULL)
    Bio::DB::HTS::Index bai
    Bio::DB::HTSfile    hfp
    int             ref
    int             start
    int             end
    int             bins
    int             maxcnt
    SV*             filter
PREINIT:
    coverage_graph  cg;
    hts_plbuf_t    *pileup;
    AV*             array;
    int             i, ret;
    bam_hdr_t      *bh;
    hts_itr_t      *iter;
    const htsFormat *format ;
CODE:
  {
      /* TODO:can we do away with this check by a move to CSI as the standard for BAM indices */
      if (end >= BAM_MAX_REGION)
      {
        format = hts_get_format( hfp ) ;
        if( format->format == bam ) //enum value from htsExactFormat from hts.h
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

      iter = sam_itr_queryi(bai, ref, start, end);

      if (items >= 8 && SvROK(filter) && SvTYPE(SvRV(filter)) == SVt_PVCV)
      {
        bam1_t *b = bam_init1();
        SV *b_sv = sv_setref_pv(newSV(sizeof b), "Bio::DB::HTS::Alignment", b);

        while ((ret = sam_itr_next(hfp, iter, b)) >= 0)
          if (invoke_sv_to_int_fun(filter, b_sv) != 0)
            hts_plbuf_push(b, pileup);

        SvREFCNT_dec(b_sv); /* b_sv's destructor will call bam_destroy1(b) */
      }
      else
      {
        bam1_t *b = bam_init1();
        while ((ret = sam_itr_next(hfp, iter, b)) >= 0)
          hts_plbuf_push(b, pileup);
        bam_destroy1(b);
      }

      hts_itr_destroy(iter);

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
bami_DESTROY(hts_idx)
  Bio::DB::HTS::Index hts_idx
  PREINIT:
    HV *assocfile = get_hv("Bio::DB::HTS::_associated_file", GV_ADD);
  CODE:
    hts_idx_destroy(hts_idx) ;
    // Now release our reference to the associated Bio::DB::HTSfile
    hv_delete(assocfile, (char *) &hts_idx, sizeof hts_idx, 0);



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
    free(str.s);
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
    if (tbx_itr_next(fp, t, iter, &str) < 0) {
        free(str.s);
        XSRETURN_EMPTY;
    }

    RETVAL = newSVpv(str.s, str.l);
    free(str.s);

  OUTPUT:
    RETVAL

void
tabix_tbx_iter_free(iter)
	Bio::DB::HTS::Tabix::Iterator iter
  CODE:
	tbx_itr_destroy(iter);


MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::VCFfile PREFIX = vcf_file_

Bio::DB::HTS::VCFfile
vcf_file_open(packname, filename, mode="r")
    char* packname
    char* filename
    char* mode
    PROTOTYPE: $$$
    CODE:
      RETVAL = bcf_open(filename, mode);
    OUTPUT:
      RETVAL

Bio::DB::HTS::Tabix
vcf_file_tbx_index_load(packname, fname)
  char* packname
  char *fname
  PROTOTYPE: $$
  CODE:
      htsFile *fp = hts_open(fname,"r");
      if ( !fp ) croak("Could not read %s\n", fname);
      enum htsExactFormat format = hts_get_format(fp)->format;
      if ( hts_close(fp) ) croak("hts_close returned non-zero status: %s\n", fname);

      if ( format != vcf ) XSRETURN_UNDEF;

      RETVAL = tbx_index_load(fname);
  OUTPUT:
      RETVAL

Bio::DB::HTS::Index
vcf_file_bcf_index_load(packname, filename)
     char* packname
     char* filename
     PROTOTYPE: $$
     CODE:
         htsFile *fp = hts_open(filename,"r");
         if ( !fp ) croak("Could not read %s\n", filename);
         enum htsExactFormat format = hts_get_format(fp)->format;
	 if ( hts_close(fp) ) croak("hts_close returned non-zero status: %s\n", filename);

	 if ( format != bcf ) XSRETURN_UNDEF;

         RETVAL = bcf_index_load(filename);

     OUTPUT:
         RETVAL

void
vcf_file_bcf_index_close(packname, bcf_idx)
     char* packname
     Bio::DB::HTS::Index bcf_idx
     PROTOTYPE: $$
     CODE:
         hts_idx_destroy(bcf_idx);


Bio::DB::HTS::VCF::Header
vcf_file_header_read(vfile)
    Bio::DB::HTS::VCFfile vfile
    PREINIT:
        bcf_hdr_t* h;
    CODE:
        h = bcf_hdr_read(vfile);
        RETVAL = h;
    OUTPUT:
        RETVAL


Bio::DB::HTS::VCF::Row
vcf_file_read1(vfile,header)
    Bio::DB::HTS::VCFfile vfile
    Bio::DB::HTS::VCF::Header header
    PREINIT:
        bcf1_t *rec;
    CODE:
        rec = bcf_init();
        if ( bcf_read(vfile, header, rec) == 0 )
        {
            bcf_unpack(rec, BCF_UN_ALL) ;
            RETVAL = rec ;
        }
        else
        {
            XSRETURN_EMPTY;
        }
    OUTPUT:
        RETVAL


SV*
vcf_file_num_variants(packname,filename)
    char* packname
    char* filename
    PROTOTYPE: $$$
    PREINIT:
        int n_records = 0;
        vcfFile* vfile;
        bcf_hdr_t* h;
        bcf1_t *rec;
    CODE:
        vfile = bcf_open(filename, "r");
        h = bcf_hdr_read(vfile);
        rec = bcf_init();

        //loop through all the lines but don't do anything with them
        while(bcf_read(vfile, h, rec) == 0)
        {
            ++n_records;
        }
        bcf_destroy(rec);
        bcf_hdr_destroy(h);
        bcf_close(vfile) ;
        RETVAL = newSViv(n_records);
    OUTPUT:
        RETVAL

Bio::DB::HTS::VCF::Iterator
vcf_file_query(packname, region, ...)
     char* packname
     char* region
     INIT:
         if( items < 4 )
           croak("Missing arguments");

         if( !(SvOK(ST(2)) && sv_isobject(ST(2))))
	   croak("Invalid index argument");

         if( !(SvOK(ST(3)) && sv_isobject(ST(3))) )
	   croak("Invalid header argument");

     CODE:
	 if ( sv_isa( ST(2), "Bio::DB::HTS::Tabix" ) ) {
	   RETVAL = tbx_itr_querys ( INT2PTR(tbx_t*, SvIV((SV *)SvRV(ST(2)))), region );
         } else if ( sv_isa( ST(2), "Bio::DB::HTS::Index" ) ) {
	   assert( sv_isa( ST(3), "Bio::DB::HTS::VCF::Header") );
	   RETVAL = bcf_itr_querys ( INT2PTR(hts_idx_t*, SvIV((SV *)SvRV(ST(2)))), INT2PTR(bcf_hdr_t*, SvIV((SV *)SvRV(ST(3)))), region );
         } else
           croak ( "Argument is not a valid index" );

         if ( RETVAL == NULL ) XSRETURN_UNDEF;

     OUTPUT:
         RETVAL


void
vcf_file_vcf_close(vfile)
    Bio::DB::HTS::VCFfile vfile
    CODE:
        bcf_close(vfile);

MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::VCF::Iterator PREFIX = vcf_

Bio::DB::HTS::VCF::Row
vcf_iter_next(iter, fp, hdr, ...)
    Bio::DB::HTS::VCF::Iterator iter
    Bio::DB::HTS::VCFfile fp
    Bio::DB::HTS::VCF::Header hdr
  PREINIT:
    kstring_t str = { 0, 0, 0 };
    bcf1_t *rec = bcf_init();
    int ret;

  INIT:
    if ( items < 4 )
      croak("Missing arguments");

    if( !(SvOK(ST(3)) && sv_isobject(ST(3))) )
      croak("Invalid index argument");

  CODE:
    if ( sv_isa( ST(3), "Bio::DB::HTS::Tabix" ) ) {
      if (tbx_itr_next(fp, INT2PTR(tbx_t*, SvIV((SV *)SvRV(ST(3)))), iter, &str) < 0 || vcf_parse1(&str, hdr, rec) < 0) {
        free(str.s);
	bcf_destroy(rec);
        XSRETURN_EMPTY;
      }

      free(str.s);

    } else if ( sv_isa( ST(3), "Bio::DB::HTS::Index" ) ) {
      if (bcf_itr_next(fp, iter, rec) < 0) {
        bcf_destroy(rec);
        XSRETURN_EMPTY;
      }
    } else
      croak ( "VCF/BCF file does not have a valid index" );

    bcf_unpack(rec, BCF_UN_ALL) ;
    RETVAL = rec;

  OUTPUT:
    RETVAL

void
vcf_iter_free(iter)
  Bio::DB::HTS::VCF::Iterator iter
  CODE:
	/* can call it also on a non-tabix index, 
	   see HTSlib synced_bcf_reader.c: bcf_sr_destroy1, _reader_fill_buffer */
	tbx_itr_destroy(iter);

MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::VCF::Header PREFIX = vcfh_

void
vcfh_DESTROY(h)
    Bio::DB::HTS::VCF::Header h
    CODE:
        bcf_hdr_destroy(h);

SV*
vcfh_version(header)
  Bio::DB::HTS::VCF::Header header
  PREINIT:
  CODE:
     RETVAL = newSVpv(bcf_hdr_get_version(header),0) ;
  OUTPUT:
     RETVAL


int
vcfh_num_samples(header)
  Bio::DB::HTS::VCF::Header header
  PREINIT:
  CODE:
     RETVAL = bcf_hdr_nsamples(header) ;
  OUTPUT:
     RETVAL


SV*
vcfh_get_sample_names(header)
    Bio::DB::HTS::VCF::Header header
    PREINIT:
        int nsamples = 0 ;
        int i ;
        AV *av_ref;
    CODE:
        av_ref = newAV();
        nsamples = bcf_hdr_nsamples(header) ;
        for (i=0 ; i<nsamples ; i++)
        {
            SV *sv_ref = newSVpv(header->samples[i], 0);
            av_push(av_ref, sv_ref);
        }
        RETVAL = newRV_noinc((SV*)av_ref);
   OUTPUT:
        RETVAL

int
vcfh_num_seqnames(header)
  Bio::DB::HTS::VCF::Header header
  PREINIT:
        int nseq = 0 ;
  CODE:
     bcf_hdr_seqnames(header, &nseq);
     RETVAL = nseq;
  OUTPUT:
     RETVAL


SV*
vcfh_get_seqnames(header)
    Bio::DB::HTS::VCF::Header header
    PREINIT:
        int nseq = 0 ;
        const char **seqnames ;
        int i = 0 ;
        AV *av_ref = newAV() ;
    CODE:
        seqnames = bcf_hdr_seqnames(header, &nseq);
        for (i = 0; i < nseq; i++)
        {
            SV *sv_ref = newSVpv(seqnames[i], 0);
            av_push(av_ref, sv_ref);
        }
        free(seqnames) ;
        RETVAL = newRV_noinc((SV*)av_ref);
   OUTPUT:
        RETVAL


SV*
vcfh_fmt_text(header)
  Bio::DB::HTS::VCF::Header header
  PREINIT:
  int len, is_bcf = 0; /* discard IDX fields */
  CODE:
    /*
     * get header formatted text
     * use deprecated (since 1.4) bcf_hdr_fmt_text to support
     * older htslib versions
     * NOTE:
     *   using bcf_hdr_format is claimed to be better (optimised for huge headers),
     *   but for htslib >= 1.4 bcf_hdr_fmt_text calls bcf_hdr_format underneath
     *   and for htslib 1.3.1 bcf_hdr_fmt_text execute the same statements
     *   as bcf_hdr_format
     */
    RETVAL = newSVpv(bcf_hdr_fmt_text(header, is_bcf, &len), 0);
  OUTPUT:
    RETVAL


MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::VCF::Row PREFIX = vcfrow_

void
vcfrow_print(row,header)
  Bio::DB::HTS::VCF::Row row
  Bio::DB::HTS::VCF::Header header
  PREINIT:
     int i ;
  CODE:
     printf("\nVCF data line:\n");
     printf("chromosome:%s\t", bcf_hdr_id2name(header,row->rid));
     printf("position:%d\t", (row->pos+1));
     printf("QUAL:%f\t", row->qual);
     printf("ID:%s\t", row->d.id );
     printf("REF:%s\n", row->d.als);
     printf("Num Alleles:%d\n", row->n_allele-1);
     for( i=1 ; i<row->n_allele ; i++ )
     {
       printf("ALT[%d]=%s\t", i, row->d.allele[i]);
     }
     printf("\nNum Filters:%d\n", row->d.n_flt);

#     printf("\nfilter:%s\t", row->d.id );
#     printf("\info:%s\n", row->d.als);
  OUTPUT:


SV*
vcfrow_chromosome(row,header)
  Bio::DB::HTS::VCF::Row row
  Bio::DB::HTS::VCF::Header header
  PREINIT:
  CODE:
     RETVAL = newSVpv(bcf_hdr_id2name(header,row->rid),0) ;
  OUTPUT:
     RETVAL


int
vcfrow_position(row)
  Bio::DB::HTS::VCF::Row row
  PREINIT:
  CODE:
     RETVAL = row->pos+1;
  OUTPUT:
     RETVAL

float
vcfrow_quality(row)
  Bio::DB::HTS::VCF::Row row
  PREINIT:
  CODE:
     RETVAL = row->qual;
  OUTPUT:
     RETVAL


SV*
vcfrow_id(row)
  Bio::DB::HTS::VCF::Row row
  PREINIT:
  CODE:
     RETVAL = newSVpv(row->d.id, 0) ;
  OUTPUT:
     RETVAL

SV*
vcfrow_reference(row)
  Bio::DB::HTS::VCF::Row row
  PREINIT:
  CODE:
     RETVAL = newSVpv(row->d.als, 0) ;
  OUTPUT:
     RETVAL


int
vcfrow_num_alleles(row)
  Bio::DB::HTS::VCF::Row row
  PREINIT:
  CODE:
     RETVAL = row->n_allele-1 ;
  OUTPUT:
     RETVAL


SV*
vcfrow_get_alleles(row)
  Bio::DB::HTS::VCF::Row row
  PREINIT:
     int i;
     AV *av_ref;
  CODE:
     av_ref = newAV();
     for (i = 1; i < row->n_allele; ++i) {
        SV *sv_ref = newSVpv(row->d.allele[i], 0);
        av_push(av_ref, sv_ref);
     }
     RETVAL = newRV_noinc((SV*)av_ref);
  OUTPUT:
     RETVAL

int
vcfrow_num_filters(row)
  Bio::DB::HTS::VCF::Row row
  PREINIT:
  CODE:
     RETVAL = row->d.n_flt ;
  OUTPUT:
     RETVAL

int
vcfrow_has_filter(row,header,filter)
  Bio::DB::HTS::VCF::Row row
  Bio::DB::HTS::VCF::Header header
  char* filter
  PREINIT:
  CODE:
     RETVAL = bcf_has_filter(header,row,filter) ;
  OUTPUT:
     RETVAL


int
vcfrow_is_snp(row)
  Bio::DB::HTS::VCF::Row row
  PREINIT:
  CODE:
     RETVAL = bcf_is_snp(row) ;
  OUTPUT:
     RETVAL


int
vcfrow_get_variant_type(row, allele_index)
  Bio::DB::HTS::VCF::Row row
  int allele_index
  CODE:
     RETVAL = bcf_get_variant_type(row, allele_index);
  OUTPUT:
     RETVAL


SV*
vcfrow_get_info_type(row,header,id)
  Bio::DB::HTS::VCF::Row row
  Bio::DB::HTS::VCF::Header header
  char* id
  PREINIT:
      bcf_info_t* info ;
  CODE:
      info = bcf_get_info(header, row, id);
      if( info == NULL )
      {
        RETVAL = newSVpv("",0);
      }
      else
      {
        switch( info->type )
        {
          case BCF_BT_FLOAT:
               RETVAL = newSVpv("Float",0);
               break ;
          case BCF_BT_NULL:
               RETVAL = newSVpv("Flag",0);
               break ;
          case BCF_BT_CHAR:
               RETVAL = newSVpv("String",0);
               break ;
          default:
               RETVAL = newSVpv("Integer",0);
        }
      }
  OUTPUT:
      RETVAL

SV*
vcfrow_get_info(row, header, ...)
  Bio::DB::HTS::VCF::Row row
  Bio::DB::HTS::VCF::Header header

  PREINIT:
      bcf_info_t* info ;
      int i = 0, avi;
      int strlength = 0;
      int* buf_i;
      float* buf_f;
      char* buf_c;
      int result;

      vdict_t *d;
      khint_t k;
      AV* row_ids;
      HV* info_data;

  INIT:
      if ( items < 2 )
	croak ( "Missing arguments" );

      row_ids = newAV();
      if ( items > 2 ) {
        if ( SvOK(ST(2)) && SvTYPE(ST(2)) == SVt_PV ) {
	  av_push ( row_ids, newSVpv ( SvPVX(ST(2)) , 0) );
        } else
	  croak ( "ID argument must be a valid string" );
      } else {
	d = (vdict_t*)header->dict[BCF_DT_ID];
	if ( d == 0 ) croak ( "Couldn't get ID dict" );

	for ( k = kh_begin(d); k != kh_end(d); ++k )
	  if ( kh_exist(d, k) && bcf_get_info(header, row, kh_key(d, k)) != NULL )
	    av_push ( row_ids, newSVpv ( kh_key(d, k), 0) );
      }

  CODE:

      info_data = newHV();

      for ( avi = 0; avi <= AvFILL(row_ids); ++avi ) {
	char* id = savepv( SvPV_nolen( *av_fetch ( row_ids, avi, 0 ) ) ) ;

        info = bcf_get_info(header, row, id);

        if( info == NULL ) { /* info null, nothing to return */
          hv_store( info_data, id, strlen(id), newSVpv("ID_NOT_FOUND", 0), 0 );
	} else {
          AV* av_ref = newAV();

          if( info->type == BCF_BT_NULL ) {
            buf_i = calloc(1, sizeof(int)) ;
            result = bcf_get_info_flag(header, row, id, &buf_i, &(info->len));

            if( result == 1 )
              av_push(av_ref, newSViv(1));
	    else
              av_push(av_ref, newSViv(0));

            free(buf_i);
          } else if( info->type == BCF_BT_FLOAT ) {
            buf_f = calloc(info->len, sizeof(float));
            result = bcf_get_info_float(header, row, id, &buf_f, &(info->len)) ;

            for( i=0 ; i<result ; i++ )
              av_push(av_ref, newSVnv(buf_f[i])) ;

            free(buf_f);
          } else if( info->type == BCF_BT_CHAR ) {
            strlength = info->len+1 ;
            buf_c = calloc(strlength, sizeof(char));
            result = bcf_get_info_string(header, row, id, &buf_c, &strlength) ;
            buf_c[info->len] = '\0' ;

            av_push(av_ref, newSVpv(buf_c,0));

            free(buf_c);
          } else if( info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32 ) {
            buf_i = calloc(info->len, sizeof(int));
            result = bcf_get_info_int32(header, row, id, &buf_i, &(info->len)) ;

            for( i=0 ; i<result ; i++ )
              av_push(av_ref, newSViv(buf_i[i])) ;

            free(buf_i);
          }

	  hv_store ( info_data, id, strlen(id), newRV_noinc((SV*)av_ref), 0 );
	}
      }

      if ( AvFILL(row_ids) == 0 ) {
	STRLEN len;
	char* key = SvPV(*av_fetch(row_ids, 0, 0), len);

	SV** svp = hv_fetch( info_data, key, len, 0);
	if (svp != NULL) {
	  RETVAL = newSVsv(*svp);

	  SvREFCNT_dec((SV*)info_data);
	} else
	  croak ("Couldn't find key");

      } else
	RETVAL = newRV_noinc((SV*)info_data);

      SvREFCNT_dec((SV*)row_ids);

  OUTPUT:
      RETVAL


SV*
vcfrow_get_format_type(row,header,id)
  Bio::DB::HTS::VCF::Row row
  Bio::DB::HTS::VCF::Header header
  char* id
  PREINIT:
      bcf_fmt_t* fmt ;
  CODE:
      fmt = bcf_get_fmt(header, row, id);
      if( fmt == NULL )
      {
        RETVAL = newSVpv("",0);
      }
      else
      {
        switch( fmt->type )
        {
          case BCF_BT_FLOAT:
               RETVAL = newSVpv("Float",0);
               break ;
          case BCF_BT_CHAR:
               RETVAL = newSVpv("String",0);
               break ;
          default:
               RETVAL = newSVpv("Integer",0);
        }
      }
  OUTPUT:
      RETVAL

SV*
vcfrow_get_format(row, header, ...)
  Bio::DB::HTS::VCF::Row row
  Bio::DB::HTS::VCF::Header header

  PREINIT:
      bcf_fmt_t* fmt;
      int i, avi;
      int* buf_i = NULL;
      float* buf_f = NULL;
      char* buf_c = NULL;
      int result;

      vdict_t *d;
      khint_t k;
      AV* row_ids;
      HV* fmt_data;

  INIT:
      if ( items < 2 )
	croak ( "Missing arguments" );

  CODE:

      if ( items > 2 ) {
        if ( SvOK(ST(2)) && SvTYPE(ST(2)) == SVt_PV ) {
	  // av_push ( row_ids, newSVpv ( SvPVX(ST(2)) , 0) );

	  fmt = bcf_get_fmt(header, row, SvPVX(ST(2)));

	  if( fmt == NULL ) { /* fmt null, nothing to return */
	    RETVAL = newSVpv("ID_NOT_FOUND", 0);
	  } else {
	    AV* av_ref = (AV*) sv_2mortal((SV*) newAV() );

	    int ndst = 0;
	    if( fmt->type == BCF_BT_FLOAT ) {
	      result = bcf_get_format_float(header, row, SvPVX(ST(2)), &buf_f, &ndst) ;
	      if ( result < 0 )
		croak ("Couldn't read float format");

	      for( i=0 ; i<ndst ; i++ )
		av_push(av_ref, newSVnv(buf_f[i]));

	      free(buf_f);

	    } else if( fmt->type == BCF_BT_CHAR ) {
	      result = bcf_get_format_char(header, row, SvPVX(ST(2)), &buf_c, &ndst) ;
	      if ( result < 0 )
		croak ("Couldn't read string format");

	      av_push(av_ref, newSVpv(buf_c, ndst+1));
	      free(buf_c);

	    } else if( fmt->type == BCF_BT_INT8 || fmt->type == BCF_BT_INT16 || fmt->type == BCF_BT_INT32 ) {
	      result = bcf_get_format_int32(header, row, SvPVX(ST(2)), &buf_i, &ndst) ;
	      if ( result < 0 )
		croak ("Couldn't read int format");

	      for( i=0 ; i<ndst ; i++ )
		av_push(av_ref, newSViv(buf_i[i]));

	      free(buf_i);
	    }

	    RETVAL = newRV((SV*)av_ref);
	  }
        } else
	  croak ( "ID argument must be a valid string" );

      } else {
	fmt_data = (HV*) sv_2mortal( (SV*) newHV() );

	d = (vdict_t*)header->dict[BCF_DT_ID];
	if ( d == 0 ) croak ( "Couldn't get ID dict" );

	for ( k = kh_begin(d); k != kh_end(d); ++k )
	  if ( kh_exist(d, k) && (fmt = bcf_get_fmt(header, row, kh_key(d, k) ) ) != NULL ) {

	    AV* av_ref = (AV*) sv_2mortal((SV*) newAV() );

	    int ndst = 0;
	    if( fmt->type == BCF_BT_FLOAT ) {
	      result = bcf_get_format_float(header, row, kh_key(d, k), &buf_f, &ndst) ;
	      if ( result < 0 )
		croak ("Couldn't read float format");

	      for( i=0 ; i<ndst ; i++ )
		av_push(av_ref, newSVnv(buf_f[i]));

	      free(buf_f);

	    } else if( fmt->type == BCF_BT_CHAR ) {
	      result = bcf_get_format_char(header, row, kh_key(d, k), &buf_c, &ndst) ;
	      if ( result < 0 )
		croak ("Couldn't read string format");

	      av_push(av_ref, newSVpv(buf_c, ndst+1));
	      free(buf_c);

	    } else if( fmt->type == BCF_BT_INT8 || fmt->type == BCF_BT_INT16 || fmt->type == BCF_BT_INT32 ) {
	      result = bcf_get_format_int32(header, row, kh_key(d, k), &buf_i, &ndst) ;
	      if ( result < 0 )
		croak ("Couldn't read int format");

	      for( i=0 ; i<ndst ; i++ )
		av_push(av_ref, newSViv(buf_i[i]));

	      free(buf_i);
	    }

	    char* key = savepv((const char*) kh_key(d, k));
	    hv_store ( fmt_data, key, (I32)strlen(key), newRV((SV*) av_ref), 0 );
	  }

	RETVAL = newRV((SV*) fmt_data);
      }

  OUTPUT:
      RETVAL

    
SV*
vcfrow_get_genotypes(row,header)
  Bio::DB::HTS::VCF::Row row
  Bio::DB::HTS::VCF::Header header
  PREINIT:
      int ngt ;
      int* gt_arr = NULL ;
      int ngt_arr = 0;
      AV* av_ref;
      int i=0 ;
  CODE:
      av_ref = newAV();
      /* Note the VCF header type treats this as a String but BCF treats as an int */
      ngt = bcf_get_genotypes(header, row, &gt_arr, &ngt_arr);
      for( i=0 ; i<ngt_arr ; i++ )
      {
        av_push(av_ref, newSViv(gt_arr[i])) ;
      }
      free(gt_arr);
      RETVAL = newRV_noinc((SV*)av_ref);
  OUTPUT:
      RETVAL

void
vcfrow_DESTROY(row)
    Bio::DB::HTS::VCF::Row row
    CODE:
      bcf_destroy(row);

MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::VCF::Sweep PREFIX = vcfs_

Bio::DB::HTS::VCF::Sweep
vcfs_sweep_open(filename)
    char* filename
    PREINIT:
        bcf_sweep_t* sweep;
    CODE:
        sweep = bcf_sweep_init(filename);
        RETVAL = sweep;
    OUTPUT:
        RETVAL

Bio::DB::HTS::VCF::HeaderPtr
vcfs_header_read(sweep)
    Bio::DB::HTS::VCF::Sweep sweep
    PREINIT:
        bcf_hdr_t* h;
    CODE:
        h = bcf_sweep_hdr(sweep);
        RETVAL = h;
    OUTPUT:
        RETVAL

Bio::DB::HTS::VCF::RowPtr
vcfs_sweep_next(sweep)
    Bio::DB::HTS::VCF::Sweep sweep
    PREINIT:
        bcf1_t* line;
    CODE:
        line = bcf_sweep_fwd(sweep);
        if( line )
        {
          RETVAL = line;
        }
        else
        {
          XSRETURN_EMPTY ;
        }
    OUTPUT:
        RETVAL

Bio::DB::HTS::VCF::RowPtr
vcfs_sweep_previous(sweep)
    Bio::DB::HTS::VCF::Sweep sweep
    PREINIT:
        bcf1_t* line;
    CODE:
        line = bcf_sweep_bwd(sweep);
        if( line )
        {
          RETVAL = line;
        }
        else
        {
          XSRETURN_EMPTY ;
        }
    OUTPUT:
        RETVAL

void
vcfs_sweep_close(sweep)
    Bio::DB::HTS::VCF::Sweep sweep
    CODE:
        bcf_sweep_destroy(sweep);
  OUTPUT:


MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::Kseq PREFIX = kseq_

Bio::DB::HTS::Kseq
kseq_new(package, filename, mode="r")
  char *package
  char *filename
  char *mode
  PROTOTYPE: $$$
  CODE:
      RETVAL = gzopen(filename, mode);
  OUTPUT:
      RETVAL

Bio::DB::HTS::Kseq
kseq_newfh(pack, fh, mode="r")
  char *pack
  PerlIO* fh
  char *mode
  PROTOTYPE: $$$
  CODE:
      RETVAL = gzdopen(PerlIO_fileno(fh), mode);
  OUTPUT:
      RETVAL

Bio::DB::HTS::Kseq::Iterator
kseq_iterator(fp)
  Bio::DB::HTS::Kseq fp
  PROTOTYPE: $
  CODE:
      RETVAL = kseq_init(fp);
  OUTPUT:
      RETVAL

void
kseq_DESTROY(fp)
  Bio::DB::HTS::Kseq fp
  PROTOTYPE: $
  CODE:
      gzclose(fp);

MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::Kseq::Kstream   PREFIX=kstream_

Bio::DB::HTS::Kseq::Kstream
kstream_new(package, fh)
  char *package
  Bio::DB::HTS::Kseq fh
  PROTOTYPE: $$
  CODE:
      RETVAL = ks_init(fh);
  OUTPUT:
      RETVAL

int
kstream_begin(kstr)
  Bio::DB::HTS::Kseq::Kstream kstr
  PROTOTYPE: $
  CODE:
      RETVAL = kstr->begin;
  OUTPUT:
      RETVAL

int
kstream_end(kstr)
  Bio::DB::HTS::Kseq::Kstream kstr
  PROTOTYPE: $
  CODE:
      RETVAL = kstr->end;
  OUTPUT:
      RETVAL

int
kstream_is_eof(kstr)
  Bio::DB::HTS::Kseq::Kstream kstr
  PROTOTYPE: $
  CODE:
      RETVAL = kstr->is_eof;
  OUTPUT:
      RETVAL

char *
kstream_buffer(kstr)
  Bio::DB::HTS::Kseq::Kstream kstr
  PROTOTYPE: $
  CODE:
      RETVAL = (char *)kstr->buf;
  OUTPUT:
      RETVAL

Bio::DB::HTS::Kseq
kstream_fh(kstr)
  Bio::DB::HTS::Kseq::Kstream kstr
  PROTOTYPE: $
  CODE:
      RETVAL = kstr->f;
  OUTPUT:
      RETVAL

void
kstream_DESTROY(kstr)
  Bio::DB::HTS::Kseq::Kstream kstr
  PROTOTYPE: $
  CODE:
      ks_destroy(kstr);

MODULE = Bio::DB::HTS PACKAGE = Bio::DB::HTS::Kseq::Iterator   PREFIX=kseqit_

SV *
kseqit_next_seq_hash(it)
  Bio::DB::HTS::Kseq::Iterator it
  PROTOTYPE: $
  INIT:
      HV * results;
  CODE:
      results = (HV *)sv_2mortal((SV *)newHV());
      if (kseq_read(it) >= 0) {
          hv_stores(results, "name", newSVpvn(it->name.s, it->name.l));
          hv_stores(results, "desc", newSVpvn(it->comment.s, it->comment.l));
          hv_stores(results, "seq", newSVpvn(it->seq.s, it->seq.l));
          hv_stores(results, "qual", newSVpvn(it->qual.s, it->qual.l));
          RETVAL = newRV((SV *)results);
      } else {
          XSRETURN_UNDEF;
      }
  OUTPUT:
      RETVAL

SV *
kseqit_next_seq(it)
Bio::DB::HTS::Kseq::Iterator it
PROTOTYPE: $
INIT:
    HV * results;
    HV * class_stash;
    SV * ref;
CODE:
    results = (HV *)sv_2mortal((SV *)newHV());
    class_stash = gv_stashpv("Bio::DB::HTS::Kseq::Record", 0);
    if (kseq_read(it) >= 0) {
        hv_stores(results, "name", newSVpvn(it->name.s, it->name.l));
        hv_stores(results, "desc", newSVpvn(it->comment.s, it->comment.l));
        hv_stores(results, "seq", newSVpvn(it->seq.s, it->seq.l));
        hv_stores(results, "qual", newSVpvn(it->qual.s, it->qual.l));
        ref = newRV((SV *)results);
        sv_bless(ref, class_stash);
        RETVAL = ref;
    } else {
        XSRETURN_UNDEF;
    }
OUTPUT:
    RETVAL

int
kseqit_read(it)
  Bio::DB::HTS::Kseq::Iterator it
  PROTOTYPE: $
  INIT:
  CODE:
      RETVAL = kseq_read(it);
  OUTPUT:
      RETVAL

void
kseqit_rewind(it)
  Bio::DB::HTS::Kseq::Iterator it
  PROTOTYPE: $
  CODE:
      /* kseq_rewind() doesn't completely rewind the file,
        just resets markers */
      kseq_rewind(it);
      /* use zlib to do so */
      gzrewind(it->f->f);

Bio::DB::HTS::Kseq::Kstream
kseqit_kstream(it)
  Bio::DB::HTS::Kseq::Iterator it
  PROTOTYPE: $
  CODE:
      RETVAL = it->f;
  OUTPUT:
      RETVAL

char *
kseqit_name(it)
  Bio::DB::HTS::Kseq::Iterator it
  PROTOTYPE: $
  CODE:
      RETVAL = it->name.s;
  OUTPUT:
      RETVAL

char *
kseqit_comment(it)
  Bio::DB::HTS::Kseq::Iterator it
  PROTOTYPE: $
  CODE:
      RETVAL = it->comment.s;
  OUTPUT:
      RETVAL

char *
kseqit_seq(it)
  Bio::DB::HTS::Kseq::Iterator it
  PROTOTYPE: $
  CODE:
      RETVAL = it->seq.s;
  OUTPUT:
      RETVAL

char *
kseqit_qual(it)
  Bio::DB::HTS::Kseq::Iterator it
  PROTOTYPE: $
  CODE:
      RETVAL = it->qual.s;
  OUTPUT:
      RETVAL

int
kseqit_last_char(it)
  Bio::DB::HTS::Kseq::Iterator it
  PROTOTYPE: $
  CODE:
      RETVAL = it->last_char;
  OUTPUT:
      RETVAL

void
kseqit_DESTROY(it)
  Bio::DB::HTS::Kseq::Iterator it
  PROTOTYPE: $
  CODE:
      kseq_destroy(it);
