#include "../params.h"
#include "../randombytes.h"
#include "../mq.h"
#include "../api.h"
#include "../sofia.h"
#include "../F.h"
#include "../sampling.h"
#ifdef SOFIA_MEASURE_ALTMQ4
#include "../mq_gf4_vertical.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#ifdef SOFIA_EXPOSE_STATICS
void GF4_scalarmul_vector(uint64_t *r, const uint64_t *a, unsigned char alpha);
void hash_transcript(unsigned char *hash, const unsigned char *transcript);
#ifdef SOFIA_MEASURE_4X_PERM_COMMIT
void commit4x(unsigned char *out0, unsigned char *out1, unsigned char *out2, unsigned char *out3,
              const uint64_t *in0, const uint64_t *in1, const uint64_t *in2, const uint64_t *in3, int len);
void permute4x(unsigned char *out0, unsigned char *out1, unsigned char *out2, unsigned char *out3,
               const uint64_t *in0, const uint64_t *in1, const uint64_t *in2, const uint64_t *in3, int len);
#else
void commit(unsigned char *out, const uint64_t *in, int inlen);
void permute(unsigned char *out, const uint64_t *in, int len);
#endif
void sample_rte(uint64_t *out, const unsigned char *seed);
#endif

#define NTESTS 100
#define MLEN 32

static unsigned long long cpucycles(void)
{
  unsigned long long result;
  __asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
    : "=a" (result) ::  "%rdx");
  return result;
}

static int cmp_llu(const void *a, const void*b)
{
    if (*(unsigned long long *)a < *(unsigned long long *)b) return -1;
    if (*(unsigned long long *)a > *(unsigned long long *)b) return 1;
    return 0;
}

static unsigned long long median(unsigned long long *l, size_t llen)
{
    qsort(l, llen, sizeof(unsigned long long), cmp_llu);

    if (llen % 2) return l[llen / 2];
    else return (l[llen/2 - 1] + l[llen/2]) / 2;
}

static unsigned long long average(unsigned long long *t, size_t tlen)
{
    unsigned long long acc=0;
    size_t i;
    for(i = 0; i < tlen; i++) {
        acc += t[i];
    }
    return acc/(tlen);
}

static void print_results(const char *s, unsigned long long *t, size_t tlen, int mult)
{
  size_t i;
  printf("%s", s);
  for (i = 0; i < tlen-1; i++) {
    t[i] = t[i+1] - t[i];
  }
  printf("\n");
  printf("median        : %llu\n", median(t, tlen));
  printf("average       : %llu\n", average(t, tlen-1));
  if (mult > 1) {
    printf("median  (%3dx): %llu\n", mult, mult*median(t, tlen));
    printf("average (%3dx): %llu\n", mult, mult*average(t, tlen-1));
  }
  printf("\n");
}

int main()
{
    unsigned long long t[NTESTS];
    uint64_t x[5*SOFIA_N64] __attribute__ ((aligned (32)));
    uint64_t y[5*SOFIA_N64] __attribute__ ((aligned (32)));
    uint64_t fx[5*SOFIA_M64] __attribute__ ((aligned (32)));
    uint64_t F[SOFIA_F64_GAPPED] __attribute__ ((aligned (32)));
    uint64_t noise[2*SOFIA_N64 + SOFIA_M64];

    unsigned char sk[CRYPTO_SECRETKEYBYTES] __attribute__ ((aligned (32)));
    unsigned char pk[CRYPTO_PUBLICKEYBYTES] __attribute__ ((aligned (32)));

    unsigned char m[MLEN];
    unsigned char sm[MLEN+CRYPTO_BYTES];
    unsigned char buf[(2*SOFIA_NBYTES + SOFIA_MBYTES) * SOFIA_ROUNDS];
    unsigned long long mlen;
    unsigned long long smlen;

    unsigned char seed[SOFIA_SEEDBYTES];
    unsigned char transcript[SOFIA_TRANSCRIPTBYTES];

    int i;

    int signing_total = 0;
    int verifying_total = 0;
    int mq_cost, g_cost, mq3_cost, g3_cost;

    randombytes(m, MLEN);
    randombytes(seed, SOFIA_SEEDBYTES);
    randombytes((unsigned char *)noise, 2*SOFIA_NBYTES + SOFIA_MBYTES);
    randombytes(transcript, SOFIA_TRANSCRIPTBYTES);

    randombytes((unsigned char *)x, 5*SOFIA_N64);
    randombytes((unsigned char *)y, 5*SOFIA_N64);
    randombytes((unsigned char *)fx, 5*SOFIA_M64);

#ifdef SOFIA_GAP_F
    expand_F_gapped(F, seed);
#else
    expand_F(F, seed);
#endif

    printf("-- api --\n\n");

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        crypto_sign_keypair(pk, sk);
    }
    print_results("sofia_keypair: ", t, NTESTS, 1);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        crypto_sign(sm, &smlen, m, MLEN, sk);
    }
    print_results("sofia_sign: ", t, NTESTS, 1);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        crypto_sign_open(m, &mlen, sm, smlen, pk);
    }
    print_results("sofia_verify: ", t, NTESTS, 1);

    printf("-- internals --\n\n");

    randombytes((unsigned char *)x, SOFIA_NBYTES);
    randombytes((unsigned char *)y, SOFIA_NBYTES);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        MQ(fx, x, F);
    }
    print_results("MQ: ", t, NTESTS, SOFIA_ROUNDS);
    mq_cost = median(t, NTESTS);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        G(fx, x, y, F);
    }
    print_results("G: ", t, NTESTS, SOFIA_ROUNDS);
    g_cost = median(t, NTESTS);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        MQ3(fx, x, F);
    }
    print_results("MQ3: ", t, NTESTS, SOFIA_ROUNDS / 3);
    mq3_cost = median(t, NTESTS);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        G3(fx, x, y, F);
    }
    print_results("G3: ", t, NTESTS, SOFIA_ROUNDS / 3);
    g3_cost = median(t, NTESTS);

#ifdef SOFIA_MEASURE_ALTMQ4
    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        mq_gf4_n128_m128(fx, x, F);
    }
    print_results("MQ1 (vert): ", t, NTESTS, SOFIA_ROUNDS);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        G_gf4_n128_m128(fx, x, y, F);
    }
    print_results("G1 (vert): ", t, NTESTS, SOFIA_ROUNDS);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        mq_gf4_n128_m128_x3(fx, x, F);
    }
    print_results("MQ3 (vert): ", t, NTESTS, SOFIA_ROUNDS / 3);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        G_gf4_n128_m128_x3(fx, x, y, F);
    }
    print_results("G3 (vert): ", t, NTESTS, SOFIA_ROUNDS / 3 );

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        mq_gf4_n128_m128_x4(fx, x, F);
    }
    print_results("MQ4 (vert): ", t, NTESTS, SOFIA_ROUNDS / 4);
#endif

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        expand_F(F, seed);
    }
    print_results("expand F: ", t, NTESTS, 1);
#ifndef SOFIA_GAP_F
    signing_total += median(t, NTESTS);
    verifying_total += median(t, NTESTS);
#else
    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        expand_F_gapped(F, seed);
    }
    print_results("expand F gapped: ", t, NTESTS, 1);
    signing_total += median(t, NTESTS);
    verifying_total += median(t, NTESTS);
#endif

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        sample_challenges(buf, SOFIA_ROUNDS, seed, SOFIA_SEEDBYTES);
    }
    print_results("sample challenges: ", t, NTESTS, 1);
    signing_total += median(t, NTESTS);
    verifying_total += median(t, NTESTS);

#ifdef SOFIA_EXPOSE_STATICS
    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        hash_transcript(buf, transcript);
    }
    print_results("hash transcript: ", t, NTESTS, 1);
    signing_total += median(t, NTESTS);
    verifying_total += median(t, NTESTS);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        sample_rte((uint64_t *)buf, seed);
    }
    print_results("sample r, t, and e: ", t, NTESTS, 1);
    signing_total += median(t, NTESTS);

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
#ifdef SOFIA_MEASURE_4X_PERM_COMMIT
        commit4x(buf, buf, buf, buf,
                 noise, noise, noise, noise, 2*SOFIA_N64 + SOFIA_M64);
        }
    print_results("commit4x(N, N, M): ", t, NTESTS, (SOFIA_ROUNDS + 3) / 4);
    signing_total += ((SOFIA_ROUNDS + 3) / 4) * median(t, NTESTS);
    verifying_total += 0.5*((SOFIA_ROUNDS + 3) / 4) * median(t, NTESTS);
#else
        commit(buf, noise, 2*SOFIA_N64 + SOFIA_M64);
    }
    print_results("commit(N, N, M): ", t, NTESTS, SOFIA_ROUNDS);
    signing_total += SOFIA_ROUNDS * median(t, NTESTS);
    verifying_total += 0.5*SOFIA_ROUNDS * median(t, NTESTS);
#endif

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
#ifdef SOFIA_MEASURE_4X_PERM_COMMIT
        commit4x(buf, buf, buf, buf,
                 noise, noise, noise, noise, SOFIA_N64 + SOFIA_M64);
    }
    print_results("commit4x(N, M): ", t, NTESTS, (SOFIA_ROUNDS + 3) / 4);
    signing_total += ((SOFIA_ROUNDS + 3) / 4) * median(t, NTESTS);
    verifying_total += 0.5*((SOFIA_ROUNDS + 3) / 4) * median(t, NTESTS);
#else
        commit(buf, noise, SOFIA_N64 + SOFIA_M64);
    }
    print_results("commit(N, N, M): ", t, NTESTS, SOFIA_ROUNDS);
    signing_total += SOFIA_ROUNDS * median(t, NTESTS);
    verifying_total += 0.5*SOFIA_ROUNDS * median(t, NTESTS);
#endif

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
#ifdef SOFIA_MEASURE_4X_PERM_COMMIT
        permute4x(buf, buf, buf, buf,
                  noise, noise, noise, noise, SOFIA_N64 + SOFIA_M64);
    }
    print_results("permute4x(N, M): ", t, NTESTS, (SOFIA_T*SOFIA_ROUNDS + 3) / 4);
    signing_total += ((SOFIA_T*SOFIA_ROUNDS + 3) / 4) * median(t, NTESTS);
    verifying_total += ((SOFIA_ROUNDS + 3) / 4) * median(t, NTESTS);
#else
        permute(buf, noise, SOFIA_N64 + SOFIA_M64);
    }
    print_results("permute(N, M): ", t, NTESTS, SOFIA_T*SOFIA_ROUNDS);
    signing_total += SOFIA_T*SOFIA_ROUNDS * median(t, NTESTS);
    verifying_total += SOFIA_ROUNDS * median(t, NTESTS);
#endif

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
#ifdef SOFIA_MEASURE_4X_PERM_COMMIT
        permute4x(buf, buf, buf, buf,
                  noise, noise, noise, noise, SOFIA_N64);
    }
    print_results("permute4x(N): ", t, NTESTS, (2*SOFIA_ROUNDS + 3) / 4);
    signing_total += ((2*SOFIA_ROUNDS + 3) / 4) * median(t, NTESTS);
    verifying_total += ((SOFIA_ROUNDS + 3) / 4) * median(t, NTESTS);
#else
        permute(buf, noise, SOFIA_N64);
    }
    print_results("permute(N): ", t, NTESTS, 2*SOFIA_ROUNDS);
    signing_total += 2*SOFIA_ROUNDS * median(t, NTESTS);
    verifying_total += SOFIA_ROUNDS * median(t, NTESTS);
#endif

    for(i=0; i<NTESTS; i++) {
        t[i] = cpucycles();
        GF4_scalarmul_vector((uint64_t *)buf, noise, 2);
    }
    print_results("scalar multiplication: ", t, NTESTS, 2*SOFIA_ROUNDS*SOFIA_T);
    signing_total += 2*SOFIA_ROUNDS*SOFIA_T * median(t, NTESTS);
    verifying_total += 1.5*SOFIA_ROUNDS * median(t, NTESTS);
#endif  // exposing statics

    printf("Total accounted for signing  : %d\n",
           (int)(SOFIA_ROUNDS * (g_cost + mq_cost) + signing_total));
    printf("Total accounted for verifying: %d\n",
           (int)(SOFIA_ROUNDS * (0.5 * g_cost + mq_cost) + verifying_total));

    printf("Total accounted for signing   (MQ3): %d\n",
           (int)(SOFIA_ROUNDS / 3 * (g3_cost + mq3_cost) +
                 (SOFIA_ROUNDS % 3) * (g_cost + mq_cost) + signing_total));
    printf("Total accounted for verifying (MQ3): %d\n",
           (int)(SOFIA_ROUNDS / 3 * (0.5 * g3_cost + mq3_cost) +
                 (SOFIA_ROUNDS % 3) * (0.5 * g_cost + mq_cost) + verifying_total));
    return 0;
}
