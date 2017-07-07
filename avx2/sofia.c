#include "sofia.h"
#include "params.h"
#include "mq.h"
#include "F.h"
#include "randombytes.h"
#include "fips202.h"
#include "fips202x4.h"
#include "sampling.h"
#include "prg.h"

#ifndef SOFIA_EXPOSE_STATICS
    #define SOFIA_STATIC static
#else
    #define SOFIA_STATIC
#endif

#ifdef SOFIA_GAP_F
    #define SOFIA_F_EXPANDING_FUNC expand_F_gapped
    #define SOFIA_F64_SIZE SOFIA_F64_GAPPED
#else
    #define SOFIA_F_EXPANDING_FUNC expand_F
    #define SOFIA_F64_SIZE SOFIA_F64
#endif

/* Multiplies two vectors of 128 bitsliced GF4 elements, and reduces.
   a and r are expended to be vectors of 'len' GF4 elements, bitsliced per 128,
    and alpha should be a non-bitsliced GF4 element. */
SOFIA_STATIC void GF4_scalarmul_vector(uint64_t *r, const uint64_t *a, unsigned char alpha)
{
    int i;
    // Create bitmasks of expanded high and low bits
    uint64_t lo = -(alpha & 1);
    uint64_t hi = -((alpha >> 1) & 1);
    uint64_t lo_xor_hi = lo ^ hi;
    uint64_t t;

    for (i = 0; i < 2; i++) {
        // equivalent to r[i+2] = (a[i+2] & lo) ^ (a[i] & hi) ^ (a[i+2] & hi);
        t = a[i];  // To ensure that r and a can point to the same location.
        r[i] = (a[i] & lo) ^ (a[i+2] & hi);
        r[i+2] = (a[i+2] & lo_xor_hi) ^ (t & hi);
    }
}

SOFIA_STATIC void hash_transcript(unsigned char *hash, const unsigned char *transcript)
{
    shake128(hash, SOFIA_HASHBYTES, transcript, SOFIA_TRANSCRIPTBYTES);
}

SOFIA_STATIC void digest(unsigned char *D, const unsigned char *m, unsigned long long mlen)
{
    shake128(D, SOFIA_HASHBYTES, m, mlen);
}

SOFIA_STATIC void commit4x(unsigned char *out0, unsigned char *out1, unsigned char *out2, unsigned char *out3,
                           const uint64_t *in0, const uint64_t *in1, const uint64_t *in2, const uint64_t *in3, int len)
{
    shake1284x(out0, out1, out2, out3, SOFIA_HASHBYTES,
               (unsigned char *)in0,
               (unsigned char *)in1,
               (unsigned char *)in2,
               (unsigned char *)in3, len*8);
}

SOFIA_STATIC void permute4x(unsigned char *out0, unsigned char *out1, unsigned char *out2, unsigned char *out3,
                            const uint64_t *in0, const uint64_t *in1, const uint64_t *in2, const uint64_t *in3, int len)
{
    shake1284x(out0, out1, out2, out3, len*8,
               (unsigned char *)in0,
               (unsigned char *)in1,
               (unsigned char *)in2,
               (unsigned char *)in3, len*8);
}

SOFIA_STATIC void sample_rte(uint64_t *out, const unsigned char *seed, int seedlen)
{
    shake128((unsigned char *)out, SOFIA_SEEDBYTES, seed, seedlen);
    prg((unsigned char *)out, (2*SOFIA_NBYTES + SOFIA_MBYTES) * SOFIA_ROUNDS, (unsigned char *)out);
}

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk)
{
#if SOFIA_Q != 4
    #error "Currently only supports SOFIA_Q == 4"
#endif
    uint64_t F[SOFIA_F64_SIZE] __attribute__ ((aligned (32)));
    uint64_t v[SOFIA_M64] __attribute__ ((aligned (32)));
    unsigned char skbuf[SOFIA_SEEDBYTES * 2] __attribute__ ((aligned (32)));
    uint64_t *s = (uint64_t *)(skbuf + SOFIA_SEEDBYTES);
    int i;

    // Expand sk to obtain a seed for F and the secret input s.
    // We also expand to obtain a value for sampling r0, t0 and e0 during
    //  signature generation, but that is not relevant here.
    randombytes(sk, SOFIA_SEEDBYTES);
    shake128(skbuf, SOFIA_SEEDBYTES + SOFIA_NBYTES, sk, SOFIA_SEEDBYTES);

    for (i = 0; i < SOFIA_SEEDBYTES; i++) {
        pk[i] = skbuf[i];
    }

    SOFIA_F_EXPANDING_FUNC(F, pk);

    MQ(v, s, F);
    // Since we cannot assume that pk is aligned, there's some indirection.
    for (i = 0; i < SOFIA_MBYTES; i++) {
        pk[SOFIA_SEEDBYTES + i] = *((unsigned char *)v + i);
    }
    return 0;
}

int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                const unsigned char *m, unsigned long long mlen,
                const unsigned char *sk)
{
#if SOFIA_Q != 4
    #error "Currently only supports SOFIA_Q == 4"
#endif
    // Seed for F, input s, seed for r, t and e.
    #define SOFIA_SKBUFBYTES 2*SOFIA_SEEDBYTES + SOFIA_NBYTES
    unsigned char skbuf[SOFIA_SKBUFBYTES] __attribute__ ((aligned (32)));
    // Transcript is put into one buffer for convenient hashing.
    unsigned char transcript[SOFIA_TRANSCRIPTBYTES] __attribute__ ((aligned (32)));
    unsigned char *pk = transcript;
    unsigned char *D = pk + SOFIA_PUBLICKEYBYTES;
    unsigned char *commits = D + SOFIA_HASHBYTES;
    unsigned char *h1 = commits + SOFIA_HASHBYTES * SOFIA_ROUNDS * 2;
    unsigned char *h2 = h1 + SOFIA_ROUNDS * (SOFIA_T * (SOFIA_NBYTES + SOFIA_MBYTES));
    // Buffer for challenges I and J, rounded up to closest byte.
    unsigned char indices[(SOFIA_ROUNDS * (2 + 1) + 7) & ~7];
    // Number of bits offset that is left when all 2-bit challenges are done.
    #define SOFIA_INDICES_OFFSET (SOFIA_ROUNDS & 0x03)
    unsigned char sample_seed[SOFIA_SEEDBYTES + SOFIA_HASHBYTES];
    uint64_t *s = (uint64_t *)(skbuf + SOFIA_SEEDBYTES);
    uint64_t rte0[(2*SOFIA_N64 + SOFIA_M64) * SOFIA_ROUNDS] __attribute__ ((aligned (32)));
    uint64_t *r0 = rte0;
    uint64_t *t0 = rte0 + SOFIA_ROUNDS * SOFIA_N64;
    uint64_t *e0 = rte0 + SOFIA_ROUNDS * 2 * SOFIA_N64;
    uint64_t r1[SOFIA_N64 * SOFIA_ROUNDS] __attribute__ ((aligned (32)));
    uint64_t Fr0[SOFIA_ROUNDS * SOFIA_M64] __attribute__ ((aligned (32)));
    uint64_t com0buf[4*(2*SOFIA_N64 + SOFIA_M64)] __attribute__ ((aligned (32)));
    uint64_t G_outbuf[SOFIA_ROUNDS * SOFIA_M64] __attribute__ ((aligned (32)));
    uint64_t com1buf[4*(SOFIA_N64 + SOFIA_M64)] __attribute__ ((aligned (32)));
    uint64_t resp1[SOFIA_ROUNDS * SOFIA_T * (SOFIA_N64 + SOFIA_M64)] __attribute__ ((aligned (32)));
    uint64_t *resp;
    uint64_t F[SOFIA_F64_SIZE] __attribute__ ((aligned (32)));
    unsigned long long m_idx;  // This is an ULL because we're iterating over m.
    int i, j, k;
    int alpha;
    int idx;

    shake128(skbuf, SOFIA_SKBUFBYTES, sk, SOFIA_SEEDBYTES);

    digest(D, m, mlen);

    SOFIA_F_EXPANDING_FUNC(F, skbuf);

    // Copy the seed of F into pk
    for (i = 0; i < SOFIA_SEEDBYTES; i++) {
        pk[i] = skbuf[i];
    }
    // Compute v = F(s)
    MQ((uint64_t *)(pk + SOFIA_SEEDBYTES), s, F);

    for (i = 0; i < SOFIA_SEEDBYTES; i++) {
        sample_seed[i] = skbuf[SOFIA_SEEDBYTES + SOFIA_NBYTES + i];
    }
    for (i = 0; i < SOFIA_HASHBYTES; i++) {
        sample_seed[i + SOFIA_SEEDBYTES] = D[i];
    }

    // Sample all the random values for r, t and e.
    // Interpret these as bitsliced elements in GF4, in 128-bit chunks.
    sample_rte(rte0, sample_seed, SOFIA_SEEDBYTES + SOFIA_HASHBYTES);

    for (i = 0; i < SOFIA_ROUNDS; i++) {
        for (j = 0; j < SOFIA_N64; j++) {
            r1[i*SOFIA_N64 + j] = s[j] ^ r0[i*SOFIA_N64 + j];
        }
    }

    for (i = 0; i < SOFIA_ROUNDS - 2; i += 3) {
        // Compute G(t_0, r_1) + e_0.
        G3(G_outbuf + i*SOFIA_M64, t0 + i*SOFIA_N64, r1 + i*SOFIA_N64, F);
        // Compute F(r_0).
        MQ3(Fr0 + i*SOFIA_M64, r0 + i*SOFIA_N64, F);
    }

    // If the number of rounds is not a multiple of 3, we have some MQ / G left
    for (i = (SOFIA_ROUNDS / 3) * 3; i < SOFIA_ROUNDS; i++) {
        // Compute G(t_0, r_1) + e_0.
        G(G_outbuf + i*SOFIA_M64, t0 + i*SOFIA_N64, r1 + i*SOFIA_N64, F);
        // Compute F(r_0).
        MQ(Fr0 + i*SOFIA_M64, r0 + i*SOFIA_N64, F);
    }

    for (i = 0; i < SOFIA_ROUNDS; i++) {
        // Commit to r_0, t_0 and e_0.
        for (j = 0; j < SOFIA_N64; j++) {
            com0buf[(i & 0x03) * (2*SOFIA_N64 + SOFIA_M64) + j              ] = r0[i*SOFIA_N64 + j];
            com0buf[(i & 0x03) * (2*SOFIA_N64 + SOFIA_M64) + j +   SOFIA_N64] = t0[i*SOFIA_N64 + j];
            com0buf[(i & 0x03) * (2*SOFIA_N64 + SOFIA_M64) + j + 2*SOFIA_N64] = e0[i*SOFIA_M64 + j];
        }
        for (j = 0; j < SOFIA_N64; j++) {
            com1buf[(i & 0x03) * (SOFIA_N64 + SOFIA_M64) + j] = r1[i*SOFIA_N64 + j];
        }
        for (j = 0; j < SOFIA_M64; j++) {
            com1buf[(i & 0x03) * (SOFIA_N64 + SOFIA_M64) + SOFIA_N64 + j] = G_outbuf[i*SOFIA_M64 + j] ^ e0[i*SOFIA_M64 + j];
        }
        if ((i & 0x03) == 3) {
            // Commit to r_0, t_0, e_0.
            commit4x(commits + (2*(i-3) + 0)*SOFIA_HASHBYTES,
                     commits + (2*(i-2) + 0)*SOFIA_HASHBYTES,
                     commits + (2*(i-1) + 0)*SOFIA_HASHBYTES,
                     commits + (2*(i-0) + 0)*SOFIA_HASHBYTES,
                     com0buf + 0*(2*SOFIA_N64 + SOFIA_M64),
                     com0buf + 1*(2*SOFIA_N64 + SOFIA_M64),
                     com0buf + 2*(2*SOFIA_N64 + SOFIA_M64),
                     com0buf + 3*(2*SOFIA_N64 + SOFIA_M64),
                     2*SOFIA_N64 + SOFIA_M64);
            // Commit to r_1 and G(t_0, r_1) + e_0.
            commit4x(commits + (2*(i-3) + 1)*SOFIA_HASHBYTES,
                     commits + (2*(i-2) + 1)*SOFIA_HASHBYTES,
                     commits + (2*(i-1) + 1)*SOFIA_HASHBYTES,
                     commits + (2*(i-0) + 1)*SOFIA_HASHBYTES,
                     com1buf + 0*(SOFIA_N64 + SOFIA_M64),
                     com1buf + 1*(SOFIA_N64 + SOFIA_M64),
                     com1buf + 2*(SOFIA_N64 + SOFIA_M64),
                     com1buf + 3*(SOFIA_N64 + SOFIA_M64),
                     SOFIA_N64 + SOFIA_M64);
        }
    }
#if ((SOFIA_ROUNDS / 4) * 4 != SOFIA_ROUNDS)
#define SOFIA_4ROUNDS_OFFSET (SOFIA_ROUNDS & 0x03)
    commit4x(commits + ((SOFIA_ROUNDS > 3)*2*(SOFIA_ROUNDS-4) + 0)*SOFIA_HASHBYTES,
             commits + ((SOFIA_ROUNDS > 2)*2*(SOFIA_ROUNDS-3) + 0)*SOFIA_HASHBYTES,
             commits + ((SOFIA_ROUNDS > 1)*2*(SOFIA_ROUNDS-2) + 0)*SOFIA_HASHBYTES,
             commits + ((SOFIA_ROUNDS > 0)*2*(SOFIA_ROUNDS-1) + 0)*SOFIA_HASHBYTES,
             com0buf + ((SOFIA_ROUNDS > 3)*(SOFIA_4ROUNDS_OFFSET + 0) & 0x03)*(2*SOFIA_N64 + SOFIA_M64),
             com0buf + ((SOFIA_ROUNDS > 2)*(SOFIA_4ROUNDS_OFFSET + 1) & 0x03)*(2*SOFIA_N64 + SOFIA_M64),
             com0buf + ((SOFIA_ROUNDS > 1)*(SOFIA_4ROUNDS_OFFSET + 2) & 0x03)*(2*SOFIA_N64 + SOFIA_M64),
             com0buf + ((SOFIA_ROUNDS > 0)*(SOFIA_4ROUNDS_OFFSET + 3) & 0x03)*(2*SOFIA_N64 + SOFIA_M64),
             2*SOFIA_N64 + SOFIA_M64);
    commit4x(commits + ((SOFIA_ROUNDS > 3)*2*(SOFIA_ROUNDS-4) + 1)*SOFIA_HASHBYTES,
             commits + ((SOFIA_ROUNDS > 2)*2*(SOFIA_ROUNDS-3) + 1)*SOFIA_HASHBYTES,
             commits + ((SOFIA_ROUNDS > 1)*2*(SOFIA_ROUNDS-2) + 1)*SOFIA_HASHBYTES,
             commits + ((SOFIA_ROUNDS > 0)*2*(SOFIA_ROUNDS-1) + 1)*SOFIA_HASHBYTES,
             com1buf + ((SOFIA_ROUNDS > 3)*(SOFIA_4ROUNDS_OFFSET + 0) & 0x03)*(SOFIA_N64 + SOFIA_M64),
             com1buf + ((SOFIA_ROUNDS > 2)*(SOFIA_4ROUNDS_OFFSET + 1) & 0x03)*(SOFIA_N64 + SOFIA_M64),
             com1buf + ((SOFIA_ROUNDS > 1)*(SOFIA_4ROUNDS_OFFSET + 2) & 0x03)*(SOFIA_N64 + SOFIA_M64),
             com1buf + ((SOFIA_ROUNDS > 0)*(SOFIA_4ROUNDS_OFFSET + 3) & 0x03)*(SOFIA_N64 + SOFIA_M64),
             SOFIA_N64 + SOFIA_M64);
#endif

    for (i = 0; i < SOFIA_ROUNDS; i++) {

        for (j = 0, alpha = 0; j < SOFIA_T; j++, alpha++) {
            resp = resp1 + ((i * SOFIA_T + j) * (SOFIA_N64 + SOFIA_M64));

            GF4_scalarmul_vector(resp, r0 + i*SOFIA_N64, alpha);  // alpha * r_0
            for (k = 0; k < SOFIA_M64; k++) {
                resp[k] ^= t0[i*SOFIA_N64 + k];  // minus t_0
            }

            GF4_scalarmul_vector(resp + SOFIA_N64, Fr0 + i*SOFIA_M64, alpha);  // alpha * F(r_0)
            for (k = 0; k < SOFIA_M64; k++) {
                resp[SOFIA_N64 + k] ^= e0[i*SOFIA_M64 + k];  // min e_0
            }
        }
    }

    for (i = 0; i < SOFIA_ROUNDS * SOFIA_T - 3; i += 4) {
        permute4x(h1 + (i+0) * (SOFIA_NBYTES + SOFIA_MBYTES),
                  h1 + (i+1) * (SOFIA_NBYTES + SOFIA_MBYTES),
                  h1 + (i+2) * (SOFIA_NBYTES + SOFIA_MBYTES),
                  h1 + (i+3) * (SOFIA_NBYTES + SOFIA_MBYTES),
                  resp1 + (i+0) * (SOFIA_N64 + SOFIA_M64),
                  resp1 + (i+1) * (SOFIA_N64 + SOFIA_M64),
                  resp1 + (i+2) * (SOFIA_N64 + SOFIA_M64),
                  resp1 + (i+3) * (SOFIA_N64 + SOFIA_M64),
                  SOFIA_N64 + SOFIA_M64);
    }
#if (((SOFIA_ROUNDS * SOFIA_T) / 4) * 4 != (SOFIA_ROUNDS * SOFIA_T))
    permute4x(h1 + (SOFIA_ROUNDS * SOFIA_T > 0)*(SOFIA_ROUNDS * SOFIA_T-1) * (SOFIA_NBYTES + SOFIA_MBYTES),
              h1 + (SOFIA_ROUNDS * SOFIA_T > 1)*(SOFIA_ROUNDS * SOFIA_T-2) * (SOFIA_NBYTES + SOFIA_MBYTES),
              h1 + (SOFIA_ROUNDS * SOFIA_T > 2)*(SOFIA_ROUNDS * SOFIA_T-3) * (SOFIA_NBYTES + SOFIA_MBYTES),
              h1 + (SOFIA_ROUNDS * SOFIA_T > 3)*(SOFIA_ROUNDS * SOFIA_T-4) * (SOFIA_NBYTES + SOFIA_MBYTES),
              resp1 + (SOFIA_ROUNDS * SOFIA_T > 0)*(SOFIA_ROUNDS * SOFIA_T-1) * (SOFIA_N64 + SOFIA_M64),
              resp1 + (SOFIA_ROUNDS * SOFIA_T > 1)*(SOFIA_ROUNDS * SOFIA_T-2) * (SOFIA_N64 + SOFIA_M64),
              resp1 + (SOFIA_ROUNDS * SOFIA_T > 2)*(SOFIA_ROUNDS * SOFIA_T-3) * (SOFIA_N64 + SOFIA_M64),
              resp1 + (SOFIA_ROUNDS * SOFIA_T > 3)*(SOFIA_ROUNDS * SOFIA_T-4) * (SOFIA_N64 + SOFIA_M64),
              SOFIA_N64 + SOFIA_M64);
#endif
    for (i = 0; i < SOFIA_ROUNDS - 1; i += 2) {
        // Commit to r_0 and r_1.
        permute4x(h2 + (2*i     + 0) * SOFIA_NBYTES,
                  h2 + (2*(i+1) + 0) * SOFIA_NBYTES,
                  h2 + (2*i     + 1) * SOFIA_NBYTES,
                  h2 + (2*(i+1) + 1) * SOFIA_NBYTES,
                  r0 + i    *SOFIA_N64,
                  r0 + (i+1)*SOFIA_N64,
                  r1 + i    *SOFIA_N64,
                  r1 + (i+1)*SOFIA_N64, SOFIA_N64);
    }
    // .. and if we have an odd number of rounds;
#if ((SOFIA_ROUNDS / 2) * 2 != SOFIA_ROUNDS)
    // Commit to r_0 and r_1.
    permute4x(h2 + (2*(SOFIA_ROUNDS-1) + 0) * SOFIA_NBYTES,
              h2 + (2*(SOFIA_ROUNDS-1) + 0) * SOFIA_NBYTES,
              h2 + (2*(SOFIA_ROUNDS-1) + 1) * SOFIA_NBYTES,
              h2 + (2*(SOFIA_ROUNDS-1) + 1) * SOFIA_NBYTES,
              r0 + (SOFIA_ROUNDS-1)*SOFIA_N64,
              r0 + (SOFIA_ROUNDS-1)*SOFIA_N64,
              r1 + (SOFIA_ROUNDS-1)*SOFIA_N64,
              r1 + (SOFIA_ROUNDS-1)*SOFIA_N64, SOFIA_N64);
#endif

    // Compute hash over transcript to include in signature, so that the
    //  verifier can compute the challenges. Without this hash, the verifier
    //  does not know which responses were opened and which were permuted.
    // Since the commitments to randomness are parts of the transcript, a hash
    //  over all commits does not need to be included separately for the
    //  verifier to be able to check. They can recreate half the commits, and
    //  the other half are included as part of the signature.
    hash_transcript(sm, transcript);

    // Derive 'reveal indices' from hash over transcript.
    sample_challenges(indices, SOFIA_ROUNDS, sm, SOFIA_HASHBYTES);
    sm += SOFIA_HASHBYTES;

    for (i = 0; i < SOFIA_ROUNDS; i++) {
        // Include challenged resp1 in signature, and h1 for others.
        idx = (indices[i >> 2] >> ((i & 0x03) << 1)) & 0x03;
        for (j = 0; j < SOFIA_T; j++) {
            if (idx == j) {
                for (k = 0; k < SOFIA_N64 + SOFIA_M64; k++) {
                    *((uint64_t *)sm + k) = resp1[(i * SOFIA_T + j) * (SOFIA_N64 + SOFIA_M64) + k];
                }
            }
            else {
                for (k = 0; k < SOFIA_NBYTES + SOFIA_MBYTES; k++) {
                    sm[k] = h1[(i * SOFIA_T + j) * (SOFIA_NBYTES + SOFIA_MBYTES) + k];
                }
            }
            sm += SOFIA_NBYTES + SOFIA_MBYTES;
        }
        // Include challenged resp2 in signature, and h2 otherwise.
        idx = (indices[(SOFIA_ROUNDS >> 2) + ((i + 2*SOFIA_INDICES_OFFSET) >> 3)] >> ((i + 2*SOFIA_INDICES_OFFSET) & 0x7)) & 0x1;
        if (idx == 0) {
            resp = r0;
        }
        else {
            resp = r1;
        }
        for (j = 0; j < 2; j++) {
            if (idx == j) {
                for (k = 0; k < SOFIA_N64; k++) {
                    *((uint64_t *)sm + k) = resp[i*SOFIA_N64 + k];
                }
            }
            else {
                for (k = 0; k < SOFIA_NBYTES; k++) {
                    sm[k] = h2[(2*i + j) * SOFIA_NBYTES + k];
                }
            }
            sm += SOFIA_NBYTES;
        }
        // Include the commits that the verifier cannot reconstruct.
        for (j = 0; j < SOFIA_HASHBYTES; j++) {
            sm[j] = commits[(1-idx)*SOFIA_HASHBYTES + j];
        }
        commits += 2*SOFIA_HASHBYTES;
        sm += SOFIA_HASHBYTES;
    }

    *smlen = SOFIA_BYTES + mlen;
    for (m_idx = mlen; m_idx > 0; m_idx--) {
        sm[m_idx-1] = m[m_idx-1];
    }

    return 0;
}

int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                     const unsigned char *sm, unsigned long long smlen,
                     const unsigned char *pk)
{
    const unsigned char *transcript_hash = sm;
    unsigned char transcript[SOFIA_TRANSCRIPTBYTES];
    unsigned char *trans_pk = transcript;
    unsigned char *D = trans_pk + SOFIA_PUBLICKEYBYTES;
    unsigned char *commits = D + SOFIA_HASHBYTES;
    unsigned char *h1 = commits + SOFIA_HASHBYTES * SOFIA_ROUNDS * 2;
    unsigned char *h2 = h1 + SOFIA_ROUNDS * (SOFIA_T * (SOFIA_NBYTES + SOFIA_MBYTES));
    // Buffer for challenges I and J, rounded up to closest byte.
    unsigned char indices[(SOFIA_ROUNDS * (2 + 1) + 7) & ~7];
    unsigned long long m_idx;
    uint64_t F[SOFIA_F64_SIZE] __attribute__ ((aligned (32)));
    // Initialized purely to silence -Wmaybe-uninitialized;
    const uint64_t *r = (uint64_t *)sm;
    const uint64_t *t = (uint64_t *)sm;
    const uint64_t *e = (uint64_t *)sm;
    const uint64_t *v = (uint64_t *)(pk + SOFIA_SEEDBYTES);
    uint64_t commit0buf[4*(2*SOFIA_N64 + SOFIA_M64)] __attribute__ ((aligned (32)));
    uint64_t commit1buf[4*(SOFIA_N64 + SOFIA_M64)] __attribute__ ((aligned (32)));
    uint64_t Fr_inbuf[3 * SOFIA_N64] __attribute__ ((aligned (32)));
    uint64_t Gt_inbuf[3 * SOFIA_N64] __attribute__ ((aligned (32)));
    uint64_t Gr_inbuf[3 * SOFIA_M64] __attribute__ ((aligned (32)));
    uint64_t G_outbuf[SOFIA_ROUNDS * SOFIA_M64] __attribute__ ((aligned (32)));
    uint64_t Fr[SOFIA_ROUNDS * SOFIA_M64] __attribute__ ((aligned (32)));
    const unsigned char *sm_precomp;
    // Store pointers to delay permutes until we have collected some.
    // We only need to store 3, as the 4th can immediately be used to permute.
    uint64_t *permute_h1_i_ptrs[3];
    unsigned char *permute_h1_o_ptrs[3];
    uint64_t *permute_h2_i_ptrs[3];
    unsigned char *permute_h2_o_ptrs[3];
    unsigned char *commit0_o_ptrs[3];
    unsigned char *commit1_o_ptrs[3];
    int G_offset = 0;
    int alpha = 0;
    int i, j, k;
    int ch;
    int c0_o_ptrs = 0;
    int c1_o_ptrs = 0;

    digest(D, sm + SOFIA_BYTES, smlen - SOFIA_BYTES);

    SOFIA_F_EXPANDING_FUNC(F, pk);

    for (i = 0; i < SOFIA_PUBLICKEYBYTES; i++) {
        trans_pk[i] = pk[i];
    }

    sample_challenges(indices, SOFIA_ROUNDS, transcript_hash, SOFIA_HASHBYTES);
    sm += SOFIA_HASHBYTES;  // Skip over H(transcript) for now.

    sm_precomp = sm;
    // To more easily use MQ3 and G3, we precompute, depending on challenges.
    for (i = 0; i < SOFIA_ROUNDS; i++) {
        ch = (indices[i >> 2] >> ((i & 0x03) << 1)) & 0x03;
        for (j = 0; j < SOFIA_T; j++) {
            if (ch == j) {
                t = (uint64_t *)sm_precomp;
            }
            sm_precomp += SOFIA_NBYTES + SOFIA_MBYTES;
        }
        ch = (indices[(SOFIA_ROUNDS >> 2) + ((i + 2*SOFIA_INDICES_OFFSET) >> 3)] >> ((i + 2*SOFIA_INDICES_OFFSET) & 0x7)) & 0x1;
        for (j = 0; j < 2; j++) {
            if (ch == j) {
                r = (uint64_t *)sm_precomp;
            }
            sm_precomp += SOFIA_NBYTES;
        }
        sm_precomp += SOFIA_HASHBYTES;

        for (j = 0; j < SOFIA_N64; j++) {
            Fr_inbuf[(i % 3) * SOFIA_N64 + j] = r[j];
        }
        // If we've collected 3 inputs to the MQ function;
        if (i % 3 == 2) {
            MQ3(Fr + (i-2)*SOFIA_M64, Fr_inbuf, F);
        }

        if (ch == 1) {
            for (j = 0; j < SOFIA_N64; j++) {
                Gt_inbuf[(G_offset % 3) * SOFIA_N64 + j] = t[j];
                Gr_inbuf[(G_offset % 3) * SOFIA_N64 + j] = r[j];
            }
            if (G_offset % 3 == 2) {
                G3(G_outbuf + (G_offset-2)*SOFIA_M64, Gt_inbuf, Gr_inbuf, F);
            }
            G_offset++;
        }
    }
    // If number of rounds is not a multiple of 3, we have some incomplete MQ
    for (i = 0; i < SOFIA_ROUNDS % 3; i++) {
        MQ(Fr + (i+((SOFIA_ROUNDS / 3) * 3))*SOFIA_M64, Fr_inbuf + i*SOFIA_N64, F);
    }
    // If there are incomplete instances of G3 left..
    for (i = 0; i < G_offset % 3; i++) {
        G(G_outbuf + (3*(G_offset/3) + i)*SOFIA_M64,
          Gt_inbuf + i*SOFIA_N64, Gr_inbuf + i*SOFIA_N64, F);
    }

    G_offset = 0;

    for (i = 0; i < SOFIA_ROUNDS; i++) {
        // Copy or compute commitments over resp1 into the transcript.
        ch = (indices[i >> 2] >> ((i & 0x03) << 1)) & 0x03;
        for (j = 0; j < SOFIA_T; j++) {
            if (ch == j) {
                if ((i & 0x03) == 3) {
                    permute4x(permute_h1_o_ptrs[0],
                              permute_h1_o_ptrs[1],
                              permute_h1_o_ptrs[2],
                              h1,
                              permute_h1_i_ptrs[0],
                              permute_h1_i_ptrs[1],
                              permute_h1_i_ptrs[2],
                              (uint64_t *)sm, SOFIA_N64 + SOFIA_M64);
                }
                else {
                    permute_h1_o_ptrs[i & 0x03] = h1;
                    permute_h1_i_ptrs[i & 0x03] = (uint64_t *)sm;
                }
                alpha = j;
                t = (uint64_t *)sm;
                e = (uint64_t *)(sm + SOFIA_NBYTES);
            }
            else {
                for (k = 0; k < SOFIA_NBYTES + SOFIA_MBYTES; k++) {
                    h1[k] = sm[k];
                }
            }
            h1 += SOFIA_NBYTES + SOFIA_MBYTES;
            sm += SOFIA_NBYTES + SOFIA_MBYTES;
        }

        // Copy or compute commitments over resp2 into the transcript.
        ch = (indices[(SOFIA_ROUNDS >> 2) + ((i + 2*SOFIA_INDICES_OFFSET) >> 3)] >> ((i + 2*SOFIA_INDICES_OFFSET) & 0x7)) & 0x1;
        for (j = 0; j < 2; j++) {
            if (ch == j) {
                if ((i & 0x03) == 3) {
                    permute4x(permute_h2_o_ptrs[0],
                              permute_h2_o_ptrs[1],
                              permute_h2_o_ptrs[2],
                              h2,
                              permute_h2_i_ptrs[0],
                              permute_h2_i_ptrs[1],
                              permute_h2_i_ptrs[2],
                              (uint64_t *)sm, SOFIA_N64);
                }
                else {
                    permute_h2_o_ptrs[i & 0x03] = h2;
                    permute_h2_i_ptrs[i & 0x03] = (uint64_t *)sm;
                }
                r = (uint64_t *)sm;
            }
            else {
                for (k = 0; k < SOFIA_NBYTES; k++) {
                    h2[k] = sm[k];
                }
            }
            h2 += SOFIA_NBYTES;
            sm += SOFIA_NBYTES;
        }
        // Add the commit that we cannot reconstruct.
        for (j = 0; j < SOFIA_HASHBYTES; j++) {
            commits[(1-ch)*SOFIA_HASHBYTES + j] = sm[j];
        }
        sm += SOFIA_HASHBYTES;
        // Reconstruct the other commit
        if (ch == 0) {
            for (j = 0; j < SOFIA_N64; j++) {
                commit0buf[j + c0_o_ptrs*(2*SOFIA_N64 + SOFIA_M64)] = r[j];
            }
            GF4_scalarmul_vector(commit0buf + SOFIA_N64 + c0_o_ptrs*(2*SOFIA_N64 + SOFIA_M64), r, alpha);
            for (j = 0; j < SOFIA_N64; j++) {
                commit0buf[SOFIA_N64 + j + c0_o_ptrs*(2*SOFIA_N64 + SOFIA_M64)] ^= t[j];
            }
            GF4_scalarmul_vector(commit0buf + 2*SOFIA_N64 + c0_o_ptrs*(2*SOFIA_N64 + SOFIA_M64),
                                 Fr + i*SOFIA_M64, alpha);
            for (j = 0; j < SOFIA_M64; j++) {
                commit0buf[2*SOFIA_N64 + j + c0_o_ptrs*(2*SOFIA_N64 + SOFIA_M64)] ^= e[j];  // a * F(r_0) - e_1
            }
            if (c0_o_ptrs < 3) {
                commit0_o_ptrs[c0_o_ptrs] = commits + ch*SOFIA_HASHBYTES;
                c0_o_ptrs++;
            }
            else {
                commit4x(commit0_o_ptrs[0],
                         commit0_o_ptrs[1],
                         commit0_o_ptrs[2],
                         commits + ch*SOFIA_HASHBYTES,
                         commit0buf + 0*(2*SOFIA_N64 + SOFIA_M64),
                         commit0buf + 1*(2*SOFIA_N64 + SOFIA_M64),
                         commit0buf + 2*(2*SOFIA_N64 + SOFIA_M64),
                         commit0buf + 3*(2*SOFIA_N64 + SOFIA_M64),
                         2*SOFIA_N64 + SOFIA_M64);
                c0_o_ptrs = 0;
            }
        }
        else if (ch == 1) {
            for (j = 0; j < SOFIA_N64; j++) {
                commit1buf[j + c1_o_ptrs*(SOFIA_N64 + SOFIA_M64)] = r[j];
            }
            for (j = 0; j < SOFIA_M64; j++) {
                commit1buf[SOFIA_N64 + j + c1_o_ptrs*(SOFIA_N64 + SOFIA_M64)] = Fr[i*SOFIA_M64 + j] ^ v[j];  // v - F(r_1)
            }
            GF4_scalarmul_vector(commit1buf + SOFIA_N64 + c1_o_ptrs*(SOFIA_N64 + SOFIA_M64),
                                 commit1buf + SOFIA_N64 + c1_o_ptrs*(SOFIA_N64 + SOFIA_M64), alpha);
            for (j = 0; j < SOFIA_M64; j++) {
                commit1buf[SOFIA_N64 + j + c1_o_ptrs*(SOFIA_N64 + SOFIA_M64)] ^= G_outbuf[G_offset*SOFIA_M64 + j] ^ e[j];
            }
            G_offset++;
            if (c1_o_ptrs < 3) {
                commit1_o_ptrs[c1_o_ptrs] = commits + ch*SOFIA_HASHBYTES;
                c1_o_ptrs++;
            }
            else {
                commit4x(commit1_o_ptrs[0],
                         commit1_o_ptrs[1],
                         commit1_o_ptrs[2],
                         commits + ch*SOFIA_HASHBYTES,
                         commit1buf + 0*(SOFIA_N64 + SOFIA_M64),
                         commit1buf + 1*(SOFIA_N64 + SOFIA_M64),
                         commit1buf + 2*(SOFIA_N64 + SOFIA_M64),
                         commit1buf + 3*(SOFIA_N64 + SOFIA_M64),
                         SOFIA_N64 + SOFIA_M64);
                c1_o_ptrs = 0;
            }
        }
        commits += 2*SOFIA_HASHBYTES;
    }

// If the number of rounds is not a multiple of 4
// We have at most 3 permute pointer pairs ready
#if SOFIA_ROUNDS & 0x03
    permute4x(permute_h1_o_ptrs[0],
              permute_h1_o_ptrs[(SOFIA_ROUNDS > 1) * 1],
              permute_h1_o_ptrs[(SOFIA_ROUNDS > 2) * 2],
              permute_h1_o_ptrs[0],  // this one is useless, but 3x == 4x
              permute_h1_i_ptrs[0],
              permute_h1_i_ptrs[(SOFIA_ROUNDS > 1) * 1],
              permute_h1_i_ptrs[(SOFIA_ROUNDS > 2) * 2],
              permute_h1_i_ptrs[0], SOFIA_N64 + SOFIA_M64);
    permute4x(permute_h2_o_ptrs[0],
              permute_h2_o_ptrs[(SOFIA_ROUNDS > 1) * 1],
              permute_h2_o_ptrs[(SOFIA_ROUNDS > 2) * 2],
              permute_h2_o_ptrs[0],  // this one is useless, but 3x == 4x
              permute_h2_i_ptrs[0],
              permute_h2_i_ptrs[(SOFIA_ROUNDS > 1) * 1],
              permute_h2_i_ptrs[(SOFIA_ROUNDS > 2) * 2],
              permute_h2_i_ptrs[0], SOFIA_N64);
#endif

    if (c0_o_ptrs > 0) {  // If we have at least one more commit to do
        commit4x(commit0_o_ptrs[0],
                 commit0_o_ptrs[(c0_o_ptrs > 1) * 1],
                 commit0_o_ptrs[(c0_o_ptrs > 2) * 2],
                 commit0_o_ptrs[0],  // this one is useless, but 3x == 4x
                 commit0buf + 0*(2*SOFIA_N64 + SOFIA_M64),
                 commit0buf + (c0_o_ptrs > 1) * 1*(2*SOFIA_N64 + SOFIA_M64),
                 commit0buf + (c0_o_ptrs > 2) * 2*(2*SOFIA_N64 + SOFIA_M64),
                 commit0buf + 0*(2*SOFIA_N64 + SOFIA_M64),
                 2*SOFIA_N64 + SOFIA_M64);
    }
    if (c1_o_ptrs > 0) {
        commit4x(commit1_o_ptrs[0],
                 commit1_o_ptrs[(c1_o_ptrs > 1) * 1],
                 commit1_o_ptrs[(c1_o_ptrs > 2) * 2],
                 commit1_o_ptrs[0],  // this one is useless, but 3x == 4x
                 commit1buf + 0*(SOFIA_N64 + SOFIA_M64),
                 commit1buf + (c1_o_ptrs > 1) * 1*(SOFIA_N64 + SOFIA_M64),
                 commit1buf + (c1_o_ptrs > 2) * 2*(SOFIA_N64 + SOFIA_M64),
                 commit1buf + 0*(SOFIA_N64 + SOFIA_M64),
                 SOFIA_N64 + SOFIA_M64);
    }

    hash_transcript(transcript, transcript);

    for (i = 0; i < SOFIA_HASHBYTES; i++) {
        if (transcript[i] != transcript_hash[i]) {
            for (m_idx = 0; m_idx < smlen - SOFIA_BYTES; m_idx++) {
                m[m_idx] = 0;
            }
            *mlen = 0;
            return -1;
        }
    }

    for (m_idx = 0; m_idx < smlen - SOFIA_BYTES; m_idx++) {
        m[m_idx] = sm[m_idx];
    }
    *mlen = smlen - SOFIA_BYTES;

    return 0;
}
