#include <immintrin.h>
#include <stdint.h>
#include <assert.h>
#include "fips202.h"

#define NROUNDS 24
#define ROL(a, offset) ((a << offset) ^ (a >> (64-offset)))

static uint64_t load64(const unsigned char *x)
{
  unsigned long long r = 0, i;

  for (i = 0; i < 8; ++i) {
    r |= (unsigned long long)x[i] << 8 * i;
  }
  return r;
}

static void store64(uint8_t *x, uint64_t u)
{
  unsigned int i;

  for(i=0; i<8; ++i) {
    x[i] = u;
    u >>= 8;
  }
}

static const __m256i rconsts[24] = {
  {0x0000000000000001ULL, 0x0000000000000001ULL, 0x0000000000000001ULL, 0x0000000000000001ULL},
  {0x0000000000008082ULL, 0x0000000000008082ULL, 0x0000000000008082ULL, 0x0000000000008082ULL},
  {0x800000000000808aULL, 0x800000000000808aULL, 0x800000000000808aULL, 0x800000000000808aULL},
  {0x8000000080008000ULL, 0x8000000080008000ULL, 0x8000000080008000ULL, 0x8000000080008000ULL},
  {0x000000000000808bULL, 0x000000000000808bULL, 0x000000000000808bULL, 0x000000000000808bULL}, 
  {0x0000000080000001ULL, 0x0000000080000001ULL, 0x0000000080000001ULL, 0x0000000080000001ULL},
  {0x8000000080008081ULL, 0x8000000080008081ULL, 0x8000000080008081ULL, 0x8000000080008081ULL}, 
  {0x8000000000008009ULL, 0x8000000000008009ULL, 0x8000000000008009ULL, 0x8000000000008009ULL},
  {0x000000000000008aULL, 0x000000000000008aULL, 0x000000000000008aULL, 0x000000000000008aULL}, 
  {0x0000000000000088ULL, 0x0000000000000088ULL, 0x0000000000000088ULL, 0x0000000000000088ULL},
  {0x0000000080008009ULL, 0x0000000080008009ULL, 0x0000000080008009ULL, 0x0000000080008009ULL}, 
  {0x000000008000000aULL, 0x000000008000000aULL, 0x000000008000000aULL, 0x000000008000000aULL},
  {0x000000008000808bULL, 0x000000008000808bULL, 0x000000008000808bULL, 0x000000008000808bULL}, 
  {0x800000000000008bULL, 0x800000000000008bULL, 0x800000000000008bULL, 0x800000000000008bULL},
  {0x8000000000008089ULL, 0x8000000000008089ULL, 0x8000000000008089ULL, 0x8000000000008089ULL}, 
  {0x8000000000008003ULL, 0x8000000000008003ULL, 0x8000000000008003ULL, 0x8000000000008003ULL},
  {0x8000000000008002ULL, 0x8000000000008002ULL, 0x8000000000008002ULL, 0x8000000000008002ULL}, 
  {0x8000000000000080ULL, 0x8000000000000080ULL, 0x8000000000000080ULL, 0x8000000000000080ULL},
  {0x000000000000800aULL, 0x000000000000800aULL, 0x000000000000800aULL, 0x000000000000800aULL}, 
  {0x800000008000000aULL, 0x800000008000000aULL, 0x800000008000000aULL, 0x800000008000000aULL},
  {0x8000000080008081ULL, 0x8000000080008081ULL, 0x8000000080008081ULL, 0x8000000080008081ULL}, 
  {0x8000000000008080ULL, 0x8000000000008080ULL, 0x8000000000008080ULL, 0x8000000000008080ULL},
  {0x0000000080000001ULL, 0x0000000080000001ULL, 0x0000000080000001ULL, 0x0000000080000001ULL}, 
  {0x8000000080008008ULL, 0x8000000080008008ULL, 0x8000000080008008ULL, 0x8000000080008008ULL}};



void KeccakF1600_StatePermute4x(__m256i *s)
{
  unsigned char n;
  __m256i B[5], t,y,d;


  for (n = 0; n < 24; ++n) 
  {
    B[0] = _mm256_xor_si256(s[0],s[5]);
    B[0] = _mm256_xor_si256(B[0],s[10]);
    B[0] = _mm256_xor_si256(B[0],s[15]);
    B[0] = _mm256_xor_si256(B[0],s[20]);


    B[1] = _mm256_xor_si256(s[1],s[ 6]);
    B[1] = _mm256_xor_si256(B[1],s[11]);
    B[1] = _mm256_xor_si256(B[1],s[16]);
    B[1] = _mm256_xor_si256(B[1],s[21]);

    B[2] = _mm256_xor_si256(s[2],s[ 7]);
    B[2] = _mm256_xor_si256(B[2],s[12]);
    B[2] = _mm256_xor_si256(B[2],s[17]);
    B[2] = _mm256_xor_si256(B[2],s[22]);

    B[3] = _mm256_xor_si256(s[3],s[ 8]);
    B[3] = _mm256_xor_si256(B[3],s[13]);
    B[3] = _mm256_xor_si256(B[3],s[18]);
    B[3] = _mm256_xor_si256(B[3],s[23]);

    B[4] = _mm256_xor_si256(s[4],s[ 9]);
    B[4] = _mm256_xor_si256(B[4],s[14]);
    B[4] = _mm256_xor_si256(B[4],s[19]);
    B[4] = _mm256_xor_si256(B[4],s[24]);
    
    /* t = B[4] ^ ROL(B[1],1); */
    t = _mm256_slli_epi64(B[1],1);
    y = _mm256_srli_epi64(B[1],63);
    t = _mm256_or_si256(t, y);
    t = _mm256_xor_si256(B[4], t);

    s[ 0] = _mm256_xor_si256(s[ 0],t);
    s[ 5] = _mm256_xor_si256(s[ 5],t);
    s[10] = _mm256_xor_si256(s[10],t);
    s[15] = _mm256_xor_si256(s[15],t);
    s[20] = _mm256_xor_si256(s[20],t);

    /* t = B[0] ^ ROL(B[2],1); */
    t = _mm256_slli_epi64(B[2],1);
    y = _mm256_srli_epi64(B[2],63);
    t = _mm256_xor_si256(t, y);
    t = _mm256_xor_si256(B[0], t);

    s[ 1] = _mm256_xor_si256(s[ 1],t);
    s[ 6] = _mm256_xor_si256(s[ 6],t);
    s[11] = _mm256_xor_si256(s[11],t);
    s[16] = _mm256_xor_si256(s[16],t);
    s[21] = _mm256_xor_si256(s[21],t);

    /* t = B[1] ^ ROL(B[3],1); */
    t = _mm256_slli_epi64(B[3],1);
    y = _mm256_srli_epi64(B[3],63);
    t = _mm256_xor_si256(t, y);
    t = _mm256_xor_si256(B[1], t);

    s[ 2] = _mm256_xor_si256(s[ 2],t);
    s[ 7] = _mm256_xor_si256(s[ 7],t);
    s[12] = _mm256_xor_si256(s[12],t);
    s[17] = _mm256_xor_si256(s[17],t);
    s[22] = _mm256_xor_si256(s[22],t);

    /* t = B[2] ^ ROL(B[4],1); */
    t = _mm256_slli_epi64(B[4],1);
    y = _mm256_srli_epi64(B[4],63);
    t = _mm256_xor_si256(t, y);
    t = _mm256_xor_si256(B[2], t);

    s[ 3] = _mm256_xor_si256(s[ 3],t);
    s[ 8] = _mm256_xor_si256(s[ 8],t);
    s[13] = _mm256_xor_si256(s[13],t);
    s[18] = _mm256_xor_si256(s[18],t);
    s[23] = _mm256_xor_si256(s[23],t);

    /* t = B[3] ^ ROL(B[0],1); */
    t = _mm256_slli_epi64(B[0],1);
    y = _mm256_srli_epi64(B[0],63);
    t = _mm256_xor_si256(t, y);
    t = _mm256_xor_si256(B[3], t);

    s[ 4] = _mm256_xor_si256(s[ 4],t);
    s[ 9] = _mm256_xor_si256(s[ 9],t);
    s[14] = _mm256_xor_si256(s[14],t);
    s[19] = _mm256_xor_si256(s[19],t);
    s[24] = _mm256_xor_si256(s[24],t);



    t = s[10];
    //s[10] = ROL(s[1], 1);
    d     = _mm256_slli_epi64(s[1],1);
    s[10] = _mm256_srli_epi64(s[1],63);
    s[10] = _mm256_xor_si256(s[10], d);

    y = s[7];
    //s[7] = ROL(t, 3);
    d     = _mm256_slli_epi64(t,3);
    s[ 7] = _mm256_srli_epi64(t,61);
    s[ 7] = _mm256_xor_si256(s[ 7], d);

    t = s[11];
    //s[11] = ROL(y, 6);
    d     = _mm256_slli_epi64(y,6);
    s[11] = _mm256_srli_epi64(y,58);
    s[11] = _mm256_xor_si256(s[11], d);

    y = s[17];
    //s[17] = ROL(t, 10);
    d     = _mm256_slli_epi64(t,10);
    s[17] = _mm256_srli_epi64(t,54);
    s[17] = _mm256_xor_si256(s[17], d);

    t = s[18];
    //s[18] = ROL(y, 15);
    d     = _mm256_slli_epi64(y,15);
    s[18] = _mm256_srli_epi64(y,49);
    s[18] = _mm256_xor_si256(s[18], d);

    y = s[3];
    //s[3] = ROL(t, 21);
    d     = _mm256_slli_epi64(t,21);
    s[ 3] = _mm256_srli_epi64(t,43);
    s[ 3] = _mm256_xor_si256(s[ 3], d);

    t = s[5];
    //s[5] = ROL(y, 28);
    d     = _mm256_slli_epi64(y,28);
    s[ 5] = _mm256_srli_epi64(y,36);
    s[ 5] = _mm256_xor_si256(s[ 5], d);

    y = s[16];
    //s[16] = ROL(t, 36);
    d     = _mm256_slli_epi64(t,36);
    s[16] = _mm256_srli_epi64(t,28);
    s[16] = _mm256_xor_si256(s[16], d);

    t = s[8];
    //s[8] = ROL(y, 45);
    d     = _mm256_slli_epi64(y,45);
    s[ 8] = _mm256_srli_epi64(y,19);
    s[ 8] = _mm256_xor_si256(s[ 8], d);

    y = s[21];
    //s[21] = ROL(t, 55);
    d     = _mm256_slli_epi64(t,55);
    s[21] = _mm256_srli_epi64(t, 9);
    s[21] = _mm256_xor_si256(s[21], d);

    t = s[24];
    //s[24] = ROL(y, 2);
    d     = _mm256_slli_epi64(y, 2);
    s[24] = _mm256_srli_epi64(y,62);
    s[24] = _mm256_xor_si256(s[24], d);

    y = s[4];
    //s[4] = ROL(t, 14);
    d     = _mm256_slli_epi64(t,14);
    s[ 4] = _mm256_srli_epi64(t,50);
    s[ 4] = _mm256_xor_si256(s[ 4], d);
    
    t = s[15];
    //s[15] = ROL(y, 27);
    d     = _mm256_slli_epi64(y,27);
    s[15] = _mm256_srli_epi64(y,37);
    s[15] = _mm256_xor_si256(s[15], d);

    y = s[23];
    //s[23] = ROL(t, 41);
    d     = _mm256_slli_epi64(t,41);
    s[23] = _mm256_srli_epi64(t,23);
    s[23] = _mm256_xor_si256(s[23], d);

    t = s[19];
    //s[19] = ROL(y, 56);
    d     = _mm256_slli_epi64(y,56);
    s[19] = _mm256_srli_epi64(y, 8);
    s[19] = _mm256_xor_si256(s[19], d);

    y = s[13];
    //s[13] = ROL(t, 8);
    d     = _mm256_slli_epi64(t, 8);
    s[13] = _mm256_srli_epi64(t,56);
    s[13] = _mm256_xor_si256(s[13], d);
    
    t = s[12];
    //s[12] = ROL(y, 25);
    d     = _mm256_slli_epi64(y,25);
    s[12] = _mm256_srli_epi64(y,39);
    s[12] = _mm256_xor_si256(s[12], d);

    y = s[2];
    //s[2] = ROL(t, 43);
    d     = _mm256_slli_epi64(t,43);
    s[ 2] = _mm256_srli_epi64(t,21);
    s[ 2] = _mm256_xor_si256(s[ 2], d);

    t = s[20];
    //s[20] = ROL(y, 62);
    d     = _mm256_slli_epi64(y,62);
    s[20] = _mm256_srli_epi64(y, 2);
    s[20] = _mm256_xor_si256(s[20], d);

    y = s[14];
    //s[14] = ROL(t, 18);
    d     = _mm256_slli_epi64(t,18);
    s[14] = _mm256_srli_epi64(t,46);
    s[14] = _mm256_xor_si256(s[14], d);
    
    t = s[22];
    //s[22] = ROL(y, 39);
    d     = _mm256_slli_epi64(y,39);
    s[22] = _mm256_srli_epi64(y,25);
    s[22] = _mm256_xor_si256(s[22], d);

    y = s[9];
    //s[9] = ROL(t, 61);
    d     = _mm256_slli_epi64(t,61);
    s[ 9] = _mm256_srli_epi64(t, 3);
    s[ 9] = _mm256_xor_si256(s[ 9], d);

    t = s[6];
    //s[6] = ROL(y, 20);
    d     = _mm256_slli_epi64(y,20);
    s[ 6] = _mm256_srli_epi64(y,44);
    s[ 6] = _mm256_xor_si256(s[ 6], d);

    //s[1] = ROL(t, 44);
    d     = _mm256_slli_epi64(t,44);
    s[ 1] = _mm256_srli_epi64(t,20);
    s[ 1] = _mm256_xor_si256(s[ 1], d);


    B[0] = s[0];
    B[1] = s[1];
    B[2] = s[2];
    B[3] = s[3];
    B[4] = s[4];

    t = _mm256_andnot_si256(B[1], B[2]);
    s[0] = _mm256_xor_si256(B[0], t);
    t = _mm256_andnot_si256(B[2], B[3]);
    s[1] = _mm256_xor_si256(B[1], t);
    t = _mm256_andnot_si256(B[3], B[4]);
    s[2] = _mm256_xor_si256(B[2], t);
    t = _mm256_andnot_si256(B[4], B[0]);
    s[3] = _mm256_xor_si256(B[3], t);
    t = _mm256_andnot_si256(B[0], B[1]);
    s[4] = _mm256_xor_si256(B[4], t);

    B[0] = s[5];
    B[1] = s[6];
    B[2] = s[7];
    B[3] = s[8];
    B[4] = s[9];

    t = _mm256_andnot_si256(B[1], B[2]);
    s[5] = _mm256_xor_si256(B[0], t);
    t = _mm256_andnot_si256(B[2], B[3]);
    s[6] = _mm256_xor_si256(B[1], t);
    t = _mm256_andnot_si256(B[3], B[4]);
    s[7] = _mm256_xor_si256(B[2], t);
    t = _mm256_andnot_si256(B[4], B[0]);
    s[8] = _mm256_xor_si256(B[3], t);
    t = _mm256_andnot_si256(B[0], B[1]);
    s[9] = _mm256_xor_si256(B[4], t);

    B[0] = s[10];
    B[1] = s[11];
    B[2] = s[12];
    B[3] = s[13];
    B[4] = s[14];
 
    t = _mm256_andnot_si256(B[1], B[2]);
    s[10] = _mm256_xor_si256(B[0], t);
    t = _mm256_andnot_si256(B[2], B[3]);
    s[11] = _mm256_xor_si256(B[1], t);
    t = _mm256_andnot_si256(B[3], B[4]);
    s[12] = _mm256_xor_si256(B[2], t);
    t = _mm256_andnot_si256(B[4], B[0]);
    s[13] = _mm256_xor_si256(B[3], t);
    t = _mm256_andnot_si256(B[0], B[1]);
    s[14] = _mm256_xor_si256(B[4], t);

    B[0] = s[15];
    B[1] = s[16];
    B[2] = s[17];
    B[3] = s[18];
    B[4] = s[19];
 
    t = _mm256_andnot_si256(B[1], B[2]);
    s[15] = _mm256_xor_si256(B[0], t);
    t = _mm256_andnot_si256(B[2], B[3]);
    s[16] = _mm256_xor_si256(B[1], t);
    t = _mm256_andnot_si256(B[3], B[4]);
    s[17] = _mm256_xor_si256(B[2], t);
    t = _mm256_andnot_si256(B[4], B[0]);
    s[18] = _mm256_xor_si256(B[3], t);
    t = _mm256_andnot_si256(B[0], B[1]);
    s[19] = _mm256_xor_si256(B[4], t);

    B[0] = s[20];
    B[1] = s[21];
    B[2] = s[22];
    B[3] = s[23];
    B[4] = s[24];

    t = _mm256_andnot_si256(B[1], B[2]);
    s[20] = _mm256_xor_si256(B[0], t);
    t = _mm256_andnot_si256(B[2], B[3]);
    s[21] = _mm256_xor_si256(B[1], t);
    t = _mm256_andnot_si256(B[3], B[4]);
    s[22] = _mm256_xor_si256(B[2], t);
    t = _mm256_andnot_si256(B[4], B[0]);
    s[23] = _mm256_xor_si256(B[3], t);
    t = _mm256_andnot_si256(B[0], B[1]);
    s[24] = _mm256_xor_si256(B[4], t);


    s[ 0] = _mm256_xor_si256(s[0], rconsts[n]);
  }
}


static void keccak_absorb4x(__m256i *s,
                          unsigned int r,
                          const unsigned char *m0, 
                          const unsigned char *m1, 
                          const unsigned char *m2, 
                          const unsigned char *m3, 
                          unsigned long long int mlen,
                          unsigned char p)
{
  unsigned long long i;
  unsigned char t0[200];
  unsigned char t1[200];
  unsigned char t2[200];
  unsigned char t3[200];

  unsigned long long *ss = (unsigned long long *)s;

 
  while (mlen >= r) 
  {
    for (i = 0; i < r / 8; ++i)
    {
      ss[4*i+0] ^= load64(m0 + 8 * i);
      ss[4*i+1] ^= load64(m1 + 8 * i);
      ss[4*i+2] ^= load64(m2 + 8 * i);
      ss[4*i+3] ^= load64(m3 + 8 * i);
    }
    
    KeccakF1600_StatePermute4x(s);
    mlen -= r;
    m0 += r;
    m1 += r;
    m2 += r;
    m3 += r;
  }

  for (i = 0; i < r; ++i)
  {
    t0[i] = 0;
    t1[i] = 0;
    t2[i] = 0;
    t3[i] = 0;
  }
  for (i = 0; i < mlen; ++i)
  {
    t0[i] = m0[i];
    t1[i] = m1[i];
    t2[i] = m2[i];
    t3[i] = m3[i];
  }

  t0[i] = p;
  t1[i] = p;
  t2[i] = p;
  t3[i] = p;

  t0[r - 1] |= 128;
  t1[r - 1] |= 128;
  t2[r - 1] |= 128;
  t3[r - 1] |= 128;

  for (i = 0; i < r / 8; ++i)
  {
    ss[4*i+0] ^= load64(t0 + 8 * i);
    ss[4*i+1] ^= load64(t1 + 8 * i);
    ss[4*i+2] ^= load64(t2 + 8 * i);
    ss[4*i+3] ^= load64(t3 + 8 * i);
  }
}


static void keccak_squeezeblocks4x(unsigned char *h0, 
                                   unsigned char *h1, 
                                   unsigned char *h2, 
                                   unsigned char *h3, 
                                   unsigned long long int nblocks,
                                   __m256i *s, 
                                   unsigned int r)
{
  unsigned int i;

  unsigned long long *ss = (unsigned long long *)s;

  while(nblocks > 0) 
  {
    KeccakF1600_StatePermute4x(s);
    for(i=0;i<(r>>3);i++)
    {
      store64(h0+8*i, ss[4*i+0]);
      store64(h1+8*i, ss[4*i+1]);
      store64(h2+8*i, ss[4*i+2]);
      store64(h3+8*i, ss[4*i+3]);
    }
    h0 += r;
    h1 += r;
    h2 += r;
    h3 += r;
    nblocks--;
  }
}


/* N is assumed to be empty; S is assumed to be exactly one character */
void cshake128_simple4x(unsigned char *output0, 
                        unsigned char *output1,
                        unsigned char *output2,
                        unsigned char *output3, unsigned long long outlen, 
                        uint16_t cstm0, 
                        uint16_t cstm1, 
                        uint16_t cstm2, 
                        uint16_t cstm3,
                        const unsigned char *in, unsigned long long inlen)
{
  __m256i s[25];
  unsigned char *sep = (unsigned char *)s;
  unsigned char t0[SHAKE128_RATE];
  unsigned char t1[SHAKE128_RATE];
  unsigned char t2[SHAKE128_RATE];
  unsigned char t3[SHAKE128_RATE];
  unsigned int i;

  for(i=0;i<25;i++)
    s[i] = _mm256_setzero_si256();

  /* Absorb customization (domain-separation) string */
  for(i=0;i<4;i++)
  {
    sep[8*i+0] = 0x01;
    sep[8*i+1] = 0xa8;
    sep[8*i+2] = 0x01;
    sep[8*i+3] = 0x00;
    sep[8*i+4] = 0x01;
    sep[8*i+5] = 16;
  }
  sep[ 6] = cstm0 & 0xff;
  sep[ 7] = cstm0 >> 8;
  sep[14] = cstm1 & 0xff;
  sep[15] = cstm1 >> 8;
  sep[22] = cstm2 & 0xff;
  sep[23] = cstm2 >> 8;
  sep[30] = cstm3 & 0xff;
  sep[31] = cstm3 >> 8;

  KeccakF1600_StatePermute4x(s);

  /* Absorb input */
  keccak_absorb4x(s, SHAKE128_RATE, in, in, in, in, inlen, 0x04);

  /* Squeeze output */
  keccak_squeezeblocks4x(output0, output1, output2, output3, outlen/SHAKE128_RATE, s, SHAKE128_RATE);
  output0 += (outlen/SHAKE128_RATE)*SHAKE128_RATE;
  output1 += (outlen/SHAKE128_RATE)*SHAKE128_RATE;
  output2 += (outlen/SHAKE128_RATE)*SHAKE128_RATE;
  output3 += (outlen/SHAKE128_RATE)*SHAKE128_RATE;

  if(outlen%SHAKE128_RATE)
  {
    keccak_squeezeblocks4x(t0, t1, t2, t3, 1, s, SHAKE128_RATE);
    for(i=0;i<outlen%SHAKE128_RATE;i++)
    {
      output0[i] = t0[i];
      output1[i] = t1[i];
      output2[i] = t2[i];
      output3[i] = t3[i];
    }
  }
}

void shake1284x(unsigned char *output0,
                unsigned char *output1,
                unsigned char *output2,
                unsigned char *output3, unsigned long long outputByteLen,
                const unsigned char *in0,
                const unsigned char *in1,
                const unsigned char *in2,
                const unsigned char *in3, unsigned long long inputByteLen)
{
  unsigned long long i;
  __m256i s[25];
  unsigned char t0[SHAKE128_RATE];
  unsigned char t1[SHAKE128_RATE];
  unsigned char t2[SHAKE128_RATE];
  unsigned char t3[SHAKE128_RATE];

  for(i=0;i<25;i++)
    s[i] = _mm256_setzero_si256();

  keccak_absorb4x(s, SHAKE128_RATE, in0, in1, in2, in3, inputByteLen, 0x1F);

  keccak_squeezeblocks4x(output0, output1, output2, output3, outputByteLen/SHAKE128_RATE, s, SHAKE128_RATE);
  output0 += (outputByteLen/SHAKE128_RATE)*SHAKE128_RATE;
  output1 += (outputByteLen/SHAKE128_RATE)*SHAKE128_RATE;
  output2 += (outputByteLen/SHAKE128_RATE)*SHAKE128_RATE;
  output3 += (outputByteLen/SHAKE128_RATE)*SHAKE128_RATE;

  if(outputByteLen%SHAKE128_RATE)
  {
    keccak_squeezeblocks4x(t0, t1, t2, t3, 1, s, SHAKE128_RATE);
    for(i=0;i<outputByteLen%SHAKE128_RATE;i++)
    {
      output0[i] = t0[i];
      output1[i] = t1[i];
      output2[i] = t2[i];
      output3[i] = t3[i];
    }
  }
}