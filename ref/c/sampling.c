#include "params.h"
#include "sampling.h"
#include "fips202.h"

/* Samples a sequence of values in {0, .., SOFIA_T}.
   Uses rejection sampling for SOFIA_T == 3. */
void sample_challenges(unsigned char *out, int rounds,
                       const unsigned char *seed, int seedlen)
{
#if SOFIA_T == 4
    shake128(out, ((rounds * (2 + 1) + 7) & ~7) >> 3, seed, seedlen);
#elif SOFIA_T == 3
    uint64_t s[25];
    unsigned char buf[SHAKE128_RATE];
    unsigned char idx;
    int sampled = 0, offset = 0, i, j;

    for (i = 0; i < 25; i++) {
        s[i] = 0;
    }

    shake128_absorb(s, seed, seedlen);

    for (i = 0; i < ((rounds * (2 + 1) + 7) & ~7) >> 3; i++) {
        out[i] = 0;
    }

    // Note that the typical number of rounds would be far below the SHAKE128
    //  rate, making it unlikely that more than one call to squeezeblocks would
    //  ever be needed.

    while (sampled < rounds) {
        shake128_squeezeblocks(buf, 1, s);
        for (i = 0; i < SHAKE128_RATE && sampled < rounds; i++) {
            for (j = 0; j < 4 && sampled < rounds; j++) {
                offset++;
                idx = (buf[i] >> (j << 1)) & 0x03;
                if (idx == 0x03) {
                    continue;
                }
                out[sampled >> 2] |= idx << ((sampled & 0x03) << 1);
                sampled++;
            }
        }
    }

    // If we have ran through exactly one block of SHAKE128 output;
    if (offset % (SHAKE128_RATE << 2) == 0) {
        shake128_squeezeblocks(buf, 1, s);
    }

    i = (offset >> 2) % SHAKE128_RATE;
    j = (offset & 0x03) << 1;
    sampled = 0;
    while (sampled < rounds) {
        for (; i < SHAKE128_RATE; i++) {
            // TODO this can be more efficient by operating on bytes instead of
            //  bits, but then the conditionals get a bit more iffy.
            for (; j < 8; j++) {
                out[(2*rounds + sampled) >> 3] |= ((buf[i] >> j) & 0x01) << ((2*rounds + sampled) & 0x07);
                sampled++;
                if (sampled == rounds) {
                    return;
                }
            }
            j = 0;
        }
        shake128_squeezeblocks(buf, 1, s);
        i = 0;
    }
#else
    #error "Sampling currently only supports SOFIA_T == 4 or SOFIA_T == 3"
#endif
}