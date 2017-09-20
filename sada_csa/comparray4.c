/* comparray4.c
   Copyright (C) 2005, K. Sadakane, all rights reserved.

   This file contains an implementation of CSA.
   For more information, see

   K. Sadakane. Compressed text databases with efficient query
     algorithms based on the compressed suffix array.
     In Proceedings 11th Annual International Symposium on Algorithms
     and Computation (ISAAC)}, LNCS v. 1969, pages 410--421, 2000.

   K. Sadakane. Succinct representations of lcp information and
     improvements in the compressed suffix arrays.
     In Proceedings 13th Annual ACM-SIAM Symposium on Discrete
     Algorithms (SODA), 2002.

   K. Sadakane. New text indexing functionalities of the compressed
     suffix arrays. Journal of Algorithms, 48(2):294--313, 2003.


   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#define USE_MMAP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comparray4.h"

#define inline
#define ENCODENUM encodegamma
#define DECODENUM decodegamma
//#define TWO rankb_w
//#define TWO2 (rankb_w*16)
//#define L rankb_w2

#define D 16
#define TBLSIZE (1<<D)
long R3[D][TBLSIZE];
long R4[TBLSIZE];
long R5[D][TBLSIZE];
long R5n[TBLSIZE],R5b[TBLSIZE],R5x[TBLSIZE];
long R6b[TBLSIZE],R6x[TBLSIZE];

//long rankb_w,rankb_m,rankb_w2;

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif

#define dprintf

inline
long getbitD(unsigned short *B, long i)
{
  long j,l,x;
  i--;
  //j = i / D;
  //l = i % D;
  j = i >> 4;
  l = i & (D-1);
  x = (B[j]<<D)+B[j+1];
  return (x >> (D-l)) & 0xffff;
}

long getbit(unsigned short *B, long i)
{
  long j,l;
  i--;
  //j = i / D;
  //l = i % D;
  j = i >> 4;
  l = i & (D-1);
  return (B[j] >> (D-1-l)) & 1;
}

long setbit(unsigned short *B, long i,long x)
{
  long j,l;
  i--;
  j = i / D;
  l = i % D;
  if (x==0) B[j] &= (~(1<<(D-1-l)));
  else if (x==1) B[j] |= (1<<(D-1-l));
  else {
    printf("error setbit x=%ld\n",x);
    exit(1);
  }
  return x;
}

long initranktables(void)
{
  unsigned short B;
  long i,j,m,r;
  long b;
#if D!=16
  error
#endif
  for (i = 0; i < TBLSIZE; i++) { /* D==16 以外では動かない */
    B = i;
    r = 0;
    for (m = 0; m < D; m++) {
      b = getbit(&B, m+1);
      r += b;
      R3[m][i] = r;
    }
    for (m = 1; m <= D; m++) {
      r = 0;
      for (j = 1; j <= D; j++) {
	b = getbit(&B, j);
	if (b == 1) {
	  r += b;
	  if (r == m) R5[m-1][i] = j-1;
	}
      }
    }
  }
  for (i = 0; i < D; i++) {
    for (j = (1<<i); j < (2<<i); j++) {
      R4[j] = D-1-i;
    }
  }
  R4[0] = D;

  return 0;
}

long blog(long x)
{
long l;
  l = 0;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

long encodegamma(unsigned short *B,long p,long x) /* x は1以上 */
{
long j,w;
  if (x<=0) {
    fprintf(stderr,"encodegamma %ld\n",x);  exit(1);
  }
  w = blog(x);
  for (j=0;j<w-1;j++) setbit(B,1+p+j,0);
  //  setbit(B,1+p+w,1);
  for (j=w-1;j>=0;j--) setbit(B,1+p+(w-1)+(w-1)-j,(x >> j)&1);
  return 2*w-1;
}

#ifndef DEBUG
inline
#endif
long getzerorun(unsigned short *B,long p)
{
  long w,w2;
#if 0
  w = 0;
  while (getbit(B,1+p+w)==0) w++;
#else
  w = 0;
  while (1) {
    w2 = R4[getbitD(B,1+p)];
    w += w2;
    if (w2 < D) break;
    p += D;
  }
#endif
  return w;
}

long decodegamma(unsigned short *B,long p,long *ans)
{
long w,x;
long w2;
#if 0
  x = getbitD(B,1+p);
  b = R6b[x];
  if (b>0) {
    *ans = R6x[x];
    return b;
  }
#endif
  w = getzerorun(B,p);
#if 0
  x = 1;
  for (i=0;i<w;i++) {
    x <<= 1;
    x += getbit(B,1+p+w+1+i);
  }
#else
  /* こっちだとpsiをgamma符号にするとダメ */
  p += w+1;
  x = 1;
  w2 = w;
  while (w2 > D) {
    x <<= D;
    x += getbitD(B,1+p);
    p += D;
    w2 -= D; /* w を変えてたので return value がおかしかった */
  }
  x <<= w2;
  x += (getbitD(B,1+p)>>(D-w2));
#endif
  *ans = x;
  return 2*w+1;
}

void mkdecodetable(void)
{
  unsigned short B[256];
  long i,j,b,b2,d,x;

  //printf("L %d\n",L);
  for (i=0; i<256; i++) B[i] = 0xffff;
  for (i = 0; i < TBLSIZE; i++) {
    B[0] = i;
    R6b[i] = 0;  R6x[i] = 0;
    b = 0;  j = 0;  x = 0;
    while (1) {
      b2 = DECODENUM(B,b,&d);
      if (b+b2 > D) break;
      b += b2;
      x += d;
      j++;
      if (j==1) {R6b[i] = b2;  R6x[i] = d;}
      //printf("i %d b %d\n",i,b);
    }
    R5n[i] = j;  R5b[i] = b;  R5x[i] = x;
    //printf("%d ",j);
  }
}

inline
long psi_list(CSA *SA,long i)
{
  long j,l,r,m;
#ifdef DEBUG
  if (i > SA->n || i < 1) {
    printf("error psi_get i=%d n=%d\n",i,SA->n);
    exit(1);
  }
#endif
  l = 1; r = SA->m;
  while (l < r) {
    m = (l + r) / 2;
    //printf("i %d m %d K %d\n",i,m,SA->K[m+1]);
    if (SA->K[m+1] <= i) {
      l = m + 1;
    } else {
      r = m;
    }
  }
  //printf("l %d r %d m %d\n",l,r,m);
  j = r;
  return j;
}

void psisort2(long *p,long *I,unsigned char *s,long n)
{
  long i,sum;
  long C[SIGMA];
  long x,c;
  //  long *J;
  for (i = 0; i < SIGMA; i++) C[i] = 0;
  for (i = 1; i <= n; i++) {
    c = s[i];
    C[c]++;
  }
  sum = 0;
  for (i = 0; i < SIGMA; i++) {
    sum = sum + C[i];
    C[i] = sum - C[i];
  }

  for (i = 0; i <= n; i++) {
    x = p[i]-1;
    if (x==0) continue;
    c = s[x];
    //printf("%d c %d C %d x\n",i,c,C[c]);
    I[1+C[c]++] = i;
  }
}

void writelong(long x,FILE *f)
{
  long tmp;
  tmp = x;
  fwrite(&tmp,sizeof(long),1,f);
}

void csa_new(long n, long *p, unsigned char *s, char *fname1, char *fname2, long rankb_w, long rankb_w2)
{
  long i,v,b,x,b2,d,w,m;
  long *I,*J;
  long K[SIGMA+2],C[SIGMA+1],C2[SIGMA+1];
  unsigned short Btmp[64];
  FILE *f1,*f2;
  long psize,isize;

  f1 = fopen(fname1,"wb"); /* psi */
  f2 = fopen(fname2,"wb"); /* directory */
  if (f1 == NULL || f2 == NULL) {
    perror("csa2_new1: ");
    exit(1);
  }

  for (i=0; i<SIGMA; ++i) {
    C[i] = 0;
  }
  for (i=0; i<n; ++i) {
    C[s[i]]++;
  }

  for (m=0,v=1,i=0; i<SIGMA; i++) {
    if (C[i]>0) {
      m++;
      C2[m] = i;
      K[m] = v;
      v += C[i];
    }
  }
  K[m+1] = v;

  for (v=0,i=0; i<SIGMA; i++) {
    v = v + C[i];
    C[i] = v;
  }

  psize = isize = 0;

  writelong(n,f2);   /* テキスト長 */
  writelong(rankb_w2,f2); /* psiを何個おきに格納するか */
  writelong(rankb_w,f2); /* SAを何個おきに格納するか */
  writelong((rankb_w*16),f2); /* ISAを何個おきに格納するか */
  writelong(SIGMA,f2);   /* アルファベットサイズ */
  writelong(m,f2);   /* 実際に現れた文字の数 */
  isize += 6*sizeof(long);

  for (i = 0; i < SIGMA; i++) {
    writelong(C[i],f2); /* 文字->辞書順 */
    //printf("C[%d] %d\n",i,C[i]);
  }
  isize += SIGMA*sizeof(long);
  for (i = 1; i <= m+1; i++) {
    writelong(K[i],f2); /* 現れた文字の累積頻度 */
    //printf("K[%d] %d\n",i,K[i]);
  }
  isize += (m+1)*sizeof(long);
  for (i = 1; i <= m; i++) {
    writelong(C2[i],f2); /* 文字の順序->文字コード */
    //printf("C2[%d] %d\n",i,C2[i]);
  }
  isize += m*sizeof(long);

  I=malloc((n+2) * sizeof(*I));
  if (I==NULL) {
    fprintf(stderr, "psi_new2 malloc I failed\n");
    exit(1);
  }

  psisort2(p,I,s-1,n);

  writelong(-1,f2); /* R[0] */
  writelong(0,f2); /* P[0] */
  isize += 2*sizeof(long);

  x = -1;  b = b2 = 0;
  for (i=1; i<=n; i++) {
    //printf("%d I %d x %d\n",i,I[i],x);
    if (I[i] < x) {
      d = (n+65536) - x;
      //printf("%d I %d d %d\n",i,I[i],d);
    } else {
      d = I[i] - x;
    }
    w = ENCODENUM(Btmp,b2,d);
    b += w;  b2 += w;
    if (b2 >= 16) {
      fwrite(Btmp,b2 / 16,sizeof(short),f1);
      psize += (b2/16)*sizeof(short);
      Btmp[0] = Btmp[b2 / 16];
      b2 = b2 % 16;
    };
    if (I[i] < x) {
      x = -1;
      i--;
    } else {
      x = I[i];
      if (i % rankb_w2 == 0) {
	writelong(I[i],f2); /* R[i / L] */
	writelong(b,f2);    /* P[i / L] */
	isize += 2*sizeof(long);
      }
    }
  }
  if (b2 > 0) {
    fwrite(Btmp,(b2+15) / 16,sizeof(short),f1);
    psize += ((b2+15)/16)*sizeof(short);
  };

  writelong(n+1,f2); /* SA[0] */
  isize += sizeof(long);
  for (i=rankb_w; i<=n; i+=rankb_w) {
    //printf("SA[%d] %d\n",i,p[i]);
    writelong(p[i],f2);
    isize += sizeof(long);
  }
  J = malloc(((n-1)/(rankb_w*16)+1)*sizeof(*J));
  if (J==NULL) {
    perror("csa2_new\n");
    exit(1);
  }
  for (i=1; i<=n; i++) {
    if ((p[i]-1) % (rankb_w*16) == 0) {
      J[(p[i]-1) / (rankb_w*16)] = i;
    }
  }
  for (i = 0; i <= (n-1)/(rankb_w*16); i++) {
    writelong(J[i],f2);
    isize += sizeof(long);
  }
  fclose(f1);
  fclose(f2);

  free(I);
  free(J);

//  printf("Psi   %d bytes (%1.3f bpc)\n",psize,(double)psize*8/n);
//  printf("Total %d bytes (%1.3f bpc)\n",psize+isize,(double)(psize+isize)*8/n);

}

long readlong(FILE *f)
{
  long tmp;
  fread(&tmp,sizeof(long),1,f);
  return tmp;
}

long csa_read(CSA *SA,char *fname1,char *fname2)
{
  long i,n,m;
  FILE *f;
  long psize,isize;
  unsigned char *ptr;

#ifndef USE_MMAP
  f = fopen(fname1,"rb");
  if (f == NULL) {
    perror("csa2_read1: ");
    exit(1);
  }
  fseek(f,0,SEEK_END);
  psize = ftell(f);
  fseek(f,0,0);
  SA->B = malloc(psize+1);
  if (SA->B == NULL) {
    perror("csa2_read2: ");
    exit(1);
  }
  fread(SA->B,psize+1,1,f);
  fclose(f);
#else
  SA->mapp = mymmap(fname1);
  if (SA->mapp->addr==NULL) {
    perror("mmap1\n");
    exit(1);
  }
  SA->B = (unsigned short *)SA->mapp->addr;
  psize = SA->mapp->len;
#endif

  f = fopen(fname2,"rb");
  if (f == NULL) {
    perror("csa2_read3: ");
    exit(1);
  }
  fseek(f,0,SEEK_END);
  isize = ftell(f);
  fseek(f,0,0);
  SA->n = n = readlong(f);   /* テキスト長 */
  SA->l = readlong(f); /* psiを何個おきに格納するか */
  SA->two = readlong(f); /* SAを何個おきに格納するか */
  SA->two2 = readlong(f); /* ISAを何個おきに格納するか */

//  printf("D=%d (stores SA for every D)\n",SA->two);
//  printf("L=%d (directory for Psi)\n",SA->l);
//  printf("Psi   %d bytes (%1.3f bpc)\n",psize,(double)psize*8/n);
//  printf("Total %d bytes (%1.3f bpc)\n",psize+isize,(double)(psize+isize)*8/n);

  if ((m=readlong(f)) != SIGMA) {   /* アルファベットサイズ */
    printf("error sigma=%ld\n",m);
  }
  SA->m = m = readlong(f);   /* 実際に現れた文字の数 */
  isize = 6*sizeof(long);

  for (i = 0; i < SIGMA; i++) {
    SA->C[i] = readlong(f); /* 文字->辞書順 */
    //printf("C[%d] %d\n",i,SA->C[i]);
  }
  isize += SIGMA*sizeof(long);
  for (i = 1; i <= m+1; i++) {
    SA->K[i] = readlong(f); /* 現れた文字の累積頻度 */
    //printf("K[%d] %d\n",i,SA->K[i]);
  }
  isize += (m+1)*sizeof(long);
  for (i = 1; i <= m; i++) {
    SA->C2[i] = readlong(f); /* 文字の順序->文字コード */
    //printf("C2[%d] %d\n",i,SA->C2[i]);
  }
  isize += m*sizeof(long);

#ifndef USE_MMAP
  SA->R = malloc((n / SA->l + 1)*2*sizeof(long));
  if (SA->R == NULL) {
    perror("csa2_read4: ");
    exit(1);
  }
  for (i = 0; i <= n / SA->l; i++) {
    SA->R[i*2] = readlong(f); /* psiの値 */
    SA->R[i*2+1] = readlong(f); /* psiへのポインタ */
  }

  SA->SA = malloc((n / SA->two + 1)*sizeof(long));
  if (SA->SA == NULL) {
    perror("csa2_read6: ");
    exit(1);
  }
  for (i = 0; i <= (n / SA->two); i++) {
    SA->SA[i] = readlong(f);
    //printf("SA[%d] %d\n",i,SA->SA[i]);
  }
  SA->ISA = malloc((n / SA->two2 + 1)*sizeof(long));
  if (SA->ISA == NULL) {
    perror("csa2_read7: ");
    exit(1);
  }
  for (i = 0; i <= (n-1) / SA->two2; i++) {
    SA->ISA[i] = readlong(f);
  }
  fclose(f);
#else
  fclose(f);

  SA->mapi = mymmap(fname2);
  if (SA->mapi->addr==NULL) {
    perror("mmap2\n");
    exit(1);
  }
  ptr = (unsigned char *)SA->mapi->addr + isize;
  SA->R = (long *)ptr;

  isize += (n / SA->l+1)*2*sizeof(long);
  ptr = (unsigned char *)SA->mapi->addr + isize;
  SA->SA = (long *)ptr;

  isize += (n / SA->two+1)*sizeof(long);
  ptr = (unsigned char *)SA->mapi->addr + isize;
  SA->ISA = (long *)ptr;
#endif
  return 0;
}

inline
long csa_psi(CSA *SA, long i)
{
  long j,k,b,d,x;
  long k2,p,n;
  long l;
  unsigned short *B;
#ifdef DEBUG
  if (i > SA->n || i < 1) {
    printf("error csa2_psi i=%d n=%d\n",i,SA->n);
    exit(1);
  }
#endif

  l = SA->l;
  x = SA->R[(i / l)*2];
  b = SA->R[(i / l)*2+1];
  j = i % l;
  //j = i & (L-1);

  n = SA->n;
  B = SA->B;

#if 0
  for (k=0; k<j; k++) {
    b += DECODENUM(B,b,&d);
    x += d;
    if (x > n) {
      //printf("i %d k %d d %d x %d n %d\n",i,k,d,x,n);
      x = -1;
      k--;
    }
    //printf("k %d j %d b %d \n",k,j,b);
  }
#else

  k = 0;
  while (k < j) {
    p = getbitD(B,1+b);
    k2 = R5n[p];
    //printf("k %d k2 %d j %d b %d\n",k,k2,j,b);
    if (k2 == 0) {
      //if (k == j) break;
      b += DECODENUM(B,b,&d);
      x += d;
      k++;
      if (x > n) {
	x = -1;
	k--;
      }
    } else {
      if (k+k2 > j) break;
      k += k2;
      b += R5b[p];
      x += R5x[p];
    }
  }

  for (; k<j; k++) {
    b += DECODENUM(B,b,&d);
    x += d;
    if (x > n) {
      x = -1;
      k--;
    }
  }
#endif
#ifdef DEBUG
  if (x < 0 || x > SA->n) {
    printf("error csa2_psi(%d) %d\n",i,x);
  }
#endif
  return x;
}

inline
long csa_T(CSA *SA,long i)
{
  long c;
  c = psi_list(SA,i);
  return SA->C2[c];
}

void csa_decode(unsigned char *p,CSA *SA,long suf,long len)
{
  long pos;
  long i;
  pos = csa_inverse(SA,suf);
  i = 0;
  while (i < len) {
    *p++ = csa_T(SA,pos);
    pos = csa_psi(SA,pos);
    i++;
  }
}

void csa_decode2(unsigned char *p,CSA *SA,long pos,long len)
{
  long i;
  i = 0;
  while (i < len) {
    *p++ = csa_T(SA,pos);
    pos = csa_psi(SA,pos);
    i++;
  }
}

void csa_decode1line(unsigned char *p,CSA *SA,long suf,long maxlen)
{
  long i,k,m,pos;
  unsigned char *tmp;

  m = maxlen*2;
  tmp = malloc(m+1);
  if (tmp==NULL) {perror("csa_decode1line");  exit(1);}

  k = suf - maxlen;  if (k <= 0) k = 1;
  pos = csa_inverse(SA,k);

  i = 0;
  while (i < m) {
    tmp[i] = csa_T(SA,pos);
    pos = csa_psi(SA,pos);
    i++;
  }
  //tmp[m] = 0;  printf("%s\n",tmp);
  for (i = suf-k;  i < m;  i++) {
    if (tmp[i] == 0x0a) {i--;  break;}
  }
  m = i;
  for (i = suf-k;  i >= 0;  i--) {
    if (tmp[i] == 0x0a) {i++;  break;}
  }
  if (m-i > maxlen) i = m-maxlen;
  while (i < m) *p++ = tmp[i++];
  *p = 0;
  free(tmp);
}

void csa_decodeall(unsigned char *p,CSA *SA)
{
  long *I;
  long i,n,pos;
  long x,b,d;
  unsigned short *B;
  n = SA->n;
  I = malloc((n+1)*sizeof(*I));
  if (I == NULL) {perror("decodeall");  exit(1);}

  B = SA->B;
  x = -1;  b = 0;
  for (i=1; i<=n; i++) {
    b += DECODENUM(B,b,&d);
    x += d;
    if (x > n) {
      x = -1;  i--;
    } else {
      I[i] = x;
    }
  }
  pos = csa_inverse(SA,1);
  for (i=1; i<=n; i++) {
    if (pos < 1 || pos > n) {
      printf("i %ld pos %ld\n",i,pos);
    }
    *p++ = csa_T(SA,pos);
    pos = I[pos];
  }
}

long csa_lookup(CSA *SA, long i)
{
  long v,two;
  v = 0;  two = SA->two;
  while (i % two !=0) {
    i = csa_psi(SA,i);
    v++;
  }
  i = i / two;
  return SA->SA[i]-v;
}

long np;
long csa_lookup2(CSA *SA, long i)
{
  long v,two;
  v = 0;  two = SA->two;
  while (i % two !=0) {
    i = csa_psi(SA,i);
    np++;
    v++;
  }
  i = i / two;
  return SA->SA[i]-v;
}

long csa_inverse(CSA *SA, long suf)
{
  long p,pos;
  long two2;

  two2 = SA->two2;

  p = ((suf-1)/two2)*two2+1;
  pos = SA->ISA[(suf-1)/two2];

  while (p < suf) {
    pos = csa_psi(SA,pos);
    p++;
  }
  return pos;
}

int longcompare(const void *pi, const void *pj)
{
  const long i = *((const long *)pi);
  const long j = *((const long *)pj);
  if (i > j)
    return 1;
  if (i < j)
    return -1;
  return 0;
}

long *csa_batchlookup(CSA *SA, long l, long r)
{
  long *I;
  long j;
  I = malloc((r-l+1+1)*sizeof(*I));
  np = 0;
  for (j=0; j<r-l+1; j++) I[1+j] = csa_lookup2(SA,l+j);
  printf("#psi %ld (%1.3f)\n",np,(double)np/(r-l+1));
  qsort(I+1, r-l+1, sizeof(long), longcompare);
  I[0] = r-l+1;
  return I;
}

unsigned long *csa_batchlookup2(CSA *SA,long l, long r)
{
  unsigned long *I; /* 答えを入れる配列 */
  long *V; /* vを入れる配列 */
  long *J; /* Iの逆を入れる配列 */
  long v;  /* 反復の深さ */
  long m;  /* psiを計算した回数(test用) */
  long q;
  long i,j;
  long two;
  long *sa;
  long f,s;

  two = SA->two;
  sa = SA->SA;

  I = malloc((r-l+1)*sizeof(*I));
  V = malloc((r-l+1+1)*sizeof(*V));
  J = malloc((r-l+1+1)*sizeof(*J));

  for (j=l; j<=r; j++) J[j-l] = -1;
  for (j=l; j<=r; j++) I[j-l] = 0;
  for (m=0,j=l; j<=r; j++) {
    //printf("%d ",j-l);
    f = 0;
    i = j;  v = 0;
    while (i % two !=0) {
      i = csa_psi(SA,i);
      v++;
      m++;
      if (l <= i && i <= r) {
        V[j-l] = v;
        J[i-l] = j;
        f = 1;
        break;
      }
    }
    if (f==0) {
      i = i / two;
      I[j-l] = sa[i]-v;
    }
  }
  for (j=l; j<=r; j++) {
    //printf("%d ",j-l);
    if (I[j-l] != 0) {
      q = j;
      while (J[q-l] != -1) {
	s = I[q-l];
	i = J[q-l];
	v = V[i-l];
	I[i-l] = s - v;
	J[q-l] = -1;
	q = i;
      }
    }
  }

  for (j=l; j<=r; j++)
    I[j-l]--;

  free(V);  free(J);
  return I;
}

long *csa_batchlookup3(CSA *SA, long l, long r, long len)
{
  long *I; /* 答えを入れる配列 */
  long *P; /* 途中の i を入れる配列 */
  long v;  /* 反復の深さ */
  long m;  /* すでに求まったSAの数 */
  long q;
  long i,j;
  long two;
  long *sa;
  long k,b,d,x,n,w;
  unsigned short *B;

  n = SA->n;
  B = SA->B;
  two = SA->two;
  sa = SA->SA;
  w = SA->l;

  I = malloc((r-l+1+1)*sizeof(*I));
  P = malloc((r-l+1+1)*sizeof(*I));
#if 1
  x = SA->R[(l / w)*2];
  b = SA->R[(l / w)*2+1];
  j = l % w;
  for (k=0; k<j; k++) {
    //printf("l %d psi %d\n",l/w+k,x);
    b += DECODENUM(B,b,&d);
    x += d;
    if (x > n) {x = -1;  k--;}
  }
  for (m = 0, q = 0, i = l; i <= r; i++) {
    //printf("i %d psi %d\n",i,x);
    if (i % two == 0) {
      I[1+m] = sa[i / two];
      m++;
    } else {
      P[q++] = x;
    }
    b += DECODENUM(B,b,&d);
    x += d;
    if (x > n) {
      x = -1;
      b += DECODENUM(B,b,&d);
      x += d;
    }
  }
  v = 1;
#else
  for (q = 0, i = l; i <= r; i++) {
    P[q++] = i;
  }
  v = 0;
  m = 0;
#endif
  while (q > 0 && v <= len) {
    //printf("v %d q %d %d\n",v,q,P[0]);
    for (k = 0, j = 0; j < q; j++) {
      i = P[j];
      //printf("j %d psi %d\n",j,i);
      if (i % two == 0) {
	I[1+m] = sa[i / two] - v;
	//printf("I %d\n",I[1+m]);
	m++;
      } else {
	P[k++] = csa_psi(SA,i);
      }
    }
    q = k;
    v++;
  }
  for (j = 0; j < q; j++) {
    I[1+m] = csa_lookup(SA,P[j])-v;
    m++;
  }
  qsort(I+1, r-l+1, sizeof(long), longcompare);
  I[0] = r-l+1;
  free(P);
  return I;
}

/* backward search */
long csa_bsearch(const unsigned char *key,long keylen,CSA *SA,long *li,long *ri)
{
  long c,h,l,r,m,ll,rr,pl,pr;
  long x,b,w,d,n,*R;
  unsigned short *B;
  long len;

  c = key[keylen-1];
  r = SA->C[c];  if (c>0) l = SA->C[c-1]+1; else l = 1;
  len = 0;
  if (l > r) goto end;
  len++;
  for (h = keylen-2; h >= 0; h--) {
    pl = l;  pr = r;
    c = key[h];
    r = SA->C[c];  if (c>0) l = SA->C[c-1]+1; else l = 1;
    if (l > r) goto end;
#if 0
    while (1) { // find maximum r such that Psi[r] <= pr
      j = csa_psi(SA,r);
      if (j <= pr) break;
      r--;
      //if (l > r) goto end;
    }
#else
#if 0
    ll = l;  rr = r;
    while (ll <= rr) {
      m = (ll + rr) / 2;
      if (csa_psi(SA,m) <= pr) ll = m+1; else rr = m-1;
    }
    r = ll-1;
#else
    R = SA->R;  B = SA->B;  w = SA->l;  n = SA->n;
    ll = l / w + 1;
    rr = r / w;
    while (ll <= rr) {
      m = (ll + rr) / 2;
      if (R[m*2] <= pr) ll = m+1; else rr = m-1;
    }
    m = (ll-1)*w;
    x = R[(m / w)*2];
    b = R[(m / w)*2+1];

#if 1
    while (m < l) {
      b += DECODENUM(B,b,&d);
      x += d;
      //if (x > n) printf("??? \n");
      if (x > n) {x = -1;  m--;}
      m++;
    }
#endif
    while (x <= pr && m <= r) {
      b += DECODENUM(B,b,&d);
      x += d;
      //if (x > n) printf("??? \n");
      m++;
    }
    r = m-1;
#endif
#endif
    //if (l > r) goto end;
#if 0
    while (1) { // find minimum l such that Psi[l] >= pl
      j = csa_psi(SA,l);
      if (j >= pl) break;
      l++;
      //if (l > r) goto end;
    }
#else
#if 0
    ll = l;  rr = r;
    while (ll <= rr) {
      m = (ll + rr) / 2;
      if (csa_psi(SA,m) >= pl) rr = m-1; else ll = m+1;
    }
    l = rr+1;
#else
    //ll = l / w + 1;
    ll = l / w;
    rr = r / w;
    while (ll <= rr) {
      m = (ll + rr) / 2;
      if (R[m*2] >= pl) rr = m-1; else ll = m+1;
    }
    m = rr*w;
    x = R[(m / w)*2];
    b = R[(m / w)*2+1];

    while (m < l) {
      b += DECODENUM(B,b,&d);
      x += d;
      if (x > n) {x = -1;  m--;}
      m++;
    }
    while (x < pl && m <= r) {
      b += DECODENUM(B,b,&d);
      x += d;
      m++;
    }
    l = m;
#endif
#endif
    if (l > r) goto end;
    len++;
  }
 end:
  *li = l;  *ri = r;
  return len;
}
