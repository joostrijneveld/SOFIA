#include "sofia.h"
#include "params.h"
#include "mq.h"
#include "F.h"
#include "randombytes.h"
#include "fips202.h"
#include "sampling.h"

#ifndef SOFIA_EXPOSE_STATICS
    #define SOFIA_STATIC
#else
    #define SOFIA_STATIC static
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

SOFIA_STATIC void hash_transcript(unsigned char *hash, const unsigned char *transcript,
                                  const unsigned char *m, unsigned long long mlen)
{
    uint64_t transcript_shakestate[25];
    unsigned long long absorbed_bytes;
    int i;

    for(i = 0; i < 25; i++) {
        transcript_shakestate[i] = 0;
    }
    absorbed_bytes = 0;
    shake128_partial_absorb(transcript_shakestate, transcript, SOFIA_TRANSCRIPTBYTES, &absorbed_bytes);
    shake128_partial_absorb(transcript_shakestate, m, mlen, &absorbed_bytes);
    shake128_close_absorb(transcript_shakestate, &absorbed_bytes);

    shake128_squeezebytes(hash, SOFIA_HASHBYTES, transcript_shakestate);
}

SOFIA_STATIC void commit(unsigned char *out, const uint64_t *in, int inlen)
{
    shake128(out, SOFIA_HASHBYTES, (unsigned char *)in, inlen*8);
}

SOFIA_STATIC void permute(unsigned char *out, const uint64_t *in, int len)
{
    shake128(out, 8*len, (unsigned char *)in, len*8);
}

SOFIA_STATIC void sample_rte(uint64_t *out, uint64_t *shakestate)
{
    shake128_squeezebytes((unsigned char *)out, (2*SOFIA_NBYTES + SOFIA_MBYTES) * SOFIA_ROUNDS, shakestate);
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
    unsigned char *commits = pk + SOFIA_PUBLICKEYBYTES;
    unsigned char *h1 = commits + SOFIA_HASHBYTES * SOFIA_ROUNDS * 2;
    unsigned char *h2 = h1 + SOFIA_ROUNDS * (SOFIA_T * (SOFIA_NBYTES + SOFIA_MBYTES));
    // Buffer for challenges I and J, rounded up to closest byte.
    unsigned char indices[(SOFIA_ROUNDS * (2 + 1) + 7) & ~7];
    // Number of bits offset that is left when all 2-bit challenges are done.
    #define SOFIA_INDICES_OFFSET (SOFIA_ROUNDS & 0x03)
    uint64_t sample_shakestate[25];
    unsigned long long absorbed_bytes;
    uint64_t *s = (uint64_t *)(skbuf + SOFIA_SEEDBYTES);
    uint64_t rte0[(2*SOFIA_N64 + SOFIA_M64) * SOFIA_ROUNDS] __attribute__ ((aligned (32)));
    uint64_t *r0 = rte0;
    uint64_t *t0 = rte0 + SOFIA_ROUNDS * SOFIA_N64;
    uint64_t *e0 = rte0 + SOFIA_ROUNDS * 2 * SOFIA_N64;
    uint64_t r1[SOFIA_N64 * SOFIA_ROUNDS] __attribute__ ((aligned (32)));
    uint64_t Fr0[SOFIA_M64] __attribute__ ((aligned (32)));
    uint64_t com0buf[2*SOFIA_N64 + SOFIA_M64] __attribute__ ((aligned (32)));
    uint64_t com1buf[SOFIA_N64 + SOFIA_M64] __attribute__ ((aligned (32)));
    uint64_t resp1[SOFIA_ROUNDS * SOFIA_T * (SOFIA_N64 + SOFIA_M64)] __attribute__ ((aligned (32)));
    uint64_t *resp;
    uint64_t F[SOFIA_F64_SIZE] __attribute__ ((aligned (32)));
    unsigned long long m_idx;  // This is an ULL because we're iterating over m.
    int i, j, k;
    int alpha;
    int idx;

    shake128(skbuf, SOFIA_SKBUFBYTES, sk, SOFIA_SEEDBYTES);

    SOFIA_F_EXPANDING_FUNC(F, skbuf);

    // Copy the seed of F into pk
    for (i = 0; i < SOFIA_SEEDBYTES; i++) {
        pk[i] = skbuf[i];
    }
    // Compute v = F(s)
    MQ((uint64_t *)(pk + SOFIA_SEEDBYTES), s, F);

    // Absorb S_rte and the message into a shake128 state
    // S_rte also serves to protect against internal collisions
    for(i = 0; i < 25; i++) {
        sample_shakestate[i] = 0;
    }
    absorbed_bytes = 0;
    shake128_partial_absorb(sample_shakestate, skbuf + SOFIA_SEEDBYTES + SOFIA_NBYTES, SOFIA_SEEDBYTES, &absorbed_bytes);
    shake128_partial_absorb(sample_shakestate, m, mlen, &absorbed_bytes);
    shake128_close_absorb(sample_shakestate, &absorbed_bytes);

    // Sample all the random values for r, t and e.
    // Interpret these as bitsliced elements in GF4, in 128-bit chunks.
    sample_rte(rte0, sample_shakestate);

    for (i = 0; i < SOFIA_ROUNDS; i++) {
        // Commit to r_0, t_0 and e_0.
        for (j = 0; j < SOFIA_N64; j++) {
            com0buf[j              ] = r0[i*SOFIA_N64 + j];
            com0buf[j +   SOFIA_N64] = t0[i*SOFIA_N64 + j];
            com0buf[j + 2*SOFIA_N64] = e0[i*SOFIA_M64 + j];
        }
        commit(commits + (2*i + 0)*SOFIA_HASHBYTES,
               com0buf, 2*SOFIA_N64 + SOFIA_M64);

        for (j = 0; j < SOFIA_N64; j++) {
            r1[i*SOFIA_N64 + j] = s[j] ^ r0[i*SOFIA_N64 + j];
        }
        for (j = 0; j < SOFIA_N64; j++) {
            com1buf[j] = r1[i*SOFIA_N64 + j];
        }
        // Compute G(t_0, r_1) + e_0.
        G(com1buf + SOFIA_N64, t0 + i*SOFIA_N64, r1 + i*SOFIA_N64, F);
        for (j = 0; j < SOFIA_M64; j++) {
            com1buf[SOFIA_N64 + j] ^= e0[i*SOFIA_M64 + j];
        }

        // Commit to r_1 and G(t_0, r_1) + e_0.
        commit(commits + (2*i + 1)*SOFIA_HASHBYTES,
               com1buf, SOFIA_N64 + SOFIA_M64);

        MQ(Fr0, r0 + i*SOFIA_N64, F);  // F(r_0)
        for (j = 0, alpha = 0; j < SOFIA_T; j++, alpha++) {
            resp = resp1 + ((i * SOFIA_T + j) * (SOFIA_N64 + SOFIA_M64));

            GF4_scalarmul_vector(resp, r0 + i*SOFIA_N64, alpha);  // alpha * r_0
            for (k = 0; k < SOFIA_M64; k++) {
                resp[k] ^= t0[i*SOFIA_N64 + k];  // minus t_0
            }

            GF4_scalarmul_vector(resp + SOFIA_N64, Fr0, alpha);  // alpha * F(r_0)
            for (k = 0; k < SOFIA_M64; k++) {
                resp[SOFIA_N64 + k] ^= e0[i*SOFIA_M64 + k];  // min e_0
            }
            // Commit to the response.
            permute(h1 + (i * SOFIA_T + j) * (SOFIA_NBYTES + SOFIA_MBYTES),
                    resp, SOFIA_N64 + SOFIA_M64);
        }
        // Commit to r_0.
        permute(h2 + (2*i + 0) * SOFIA_NBYTES, r0 + i*SOFIA_N64, SOFIA_N64);
        // Commit to r_1.
        permute(h2 + (2*i + 1) * SOFIA_NBYTES, r1 + i*SOFIA_N64, SOFIA_N64);
    }

    // Compute hash over transcript to include in signature, so that the
    //  verifier can compute the challenges. Without this hash, the verifier
    //  does not know which responses were opened and which were permuted.
    // Since the commitments to randomness are parts of the transcript, a hash
    //  over all commits does not need to be included separately for the
    //  verifier to be able to check. They can recreate half the commits, and
    //  the other half are included as part of the signature.
    // The randomness included in the transcript is necessary and sufficient
    //  to guarantee collision resilience.
    hash_transcript(sm, transcript, m, mlen);

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
#if SOFIA_Q != 4
    #error "Currently only supports SOFIA_Q == 4"
#endif
    const unsigned char *transcript_hash = sm;
    unsigned char transcript[SOFIA_TRANSCRIPTBYTES];
    unsigned char *trans_pk = transcript;
    unsigned char *commits = trans_pk + SOFIA_PUBLICKEYBYTES;
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
    uint64_t commitbuf[2*SOFIA_N64 + SOFIA_M64] __attribute__ ((aligned (32)));
    uint64_t G_buf[SOFIA_M64 + SOFIA_N64] __attribute__ ((aligned (32)));
    int alpha = 0;
    int i, j, k;
    int ch;

    SOFIA_F_EXPANDING_FUNC(F, pk);

    for (i = 0; i < SOFIA_PUBLICKEYBYTES; i++) {
        trans_pk[i] = pk[i];
    }

    sample_challenges(indices, SOFIA_ROUNDS, transcript_hash, SOFIA_HASHBYTES);
    sm += SOFIA_HASHBYTES;  // Skip over H(transcript) for now.

    for (i = 0; i < SOFIA_ROUNDS; i++) {
        // Copy or compute commitments over resp1 into the transcript.
        ch = (indices[i >> 2] >> ((i & 0x03) << 1)) & 0x03;
        for (j = 0; j < SOFIA_T; j++) {
            if (ch == j) {
                permute(h1, (uint64_t *)sm, SOFIA_N64 + SOFIA_M64);
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
                permute(h2, (uint64_t *)sm, SOFIA_N64);
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
        for (j = 0; j < SOFIA_N64; j++) {
            commitbuf[j] = r[j];
        }
        if (ch == 0) {
            GF4_scalarmul_vector(commitbuf + SOFIA_N64, r, alpha);
            for (j = 0; j < SOFIA_N64; j++) {
                commitbuf[SOFIA_N64 + j] ^= t[j];
            }
            MQ(commitbuf + 2*SOFIA_N64, commitbuf, F);  // F(r_0)
            GF4_scalarmul_vector(commitbuf + 2*SOFIA_N64,
                                 commitbuf + 2*SOFIA_N64, alpha);
            for (j = 0; j < SOFIA_M64; j++) {
                commitbuf[2*SOFIA_N64 + j] ^= e[j];  // a * F(r_0) - e_1
            }
            commit(commits + ch*SOFIA_HASHBYTES, commitbuf, 2*SOFIA_N64 + SOFIA_M64);
        }
        else if (ch == 1) {
            MQ(commitbuf + SOFIA_N64, commitbuf, F);  // F(r_1)
            for (j = 0; j < SOFIA_M64; j++) {
                commitbuf[SOFIA_N64 + j] ^= v[j];  // v - F(r_1)
            }
            GF4_scalarmul_vector(commitbuf + SOFIA_N64,
                                 commitbuf + SOFIA_N64, alpha);
            for (j = 0; j < SOFIA_N64; j++) {
                G_buf[SOFIA_M64 + j] = t[j];  // Ensure alignment of t.
            }
            G(G_buf, G_buf + SOFIA_M64, commitbuf, F);  // G(t1, r1)
            for (j = 0; j < SOFIA_M64; j++) {
                commitbuf[SOFIA_N64 + j] ^= G_buf[j] ^ e[j];
            }
            commit(commits + ch*SOFIA_HASHBYTES, commitbuf, SOFIA_N64 + SOFIA_M64);
        }
        commits += 2*SOFIA_HASHBYTES;
    }

    hash_transcript(transcript, transcript, sm, smlen - SOFIA_BYTES);

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
