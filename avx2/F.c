#include <stdint.h>
#include "mq.h"
#include "fips202.h"
#include "fips202x4.h"
#include "params.h"

void expand_F(uint64_t *F, const unsigned char *seed)
{
#if SOFIA_N != 128
    #error "Dividing F in blocks of 128 depends on N = 128"
#endif
    cshake128_simple4x((unsigned char *)F + 0 * (SOFIA_FBYTES >> 2),
                       (unsigned char *)F + 1 * (SOFIA_FBYTES >> 2),
                       (unsigned char *)F + 2 * (SOFIA_FBYTES >> 2),
                       (unsigned char *)F + 3 * (SOFIA_FBYTES >> 2),
                       SOFIA_FBYTES >> 2,
                       0, 1, 2, 3,
                       seed, SOFIA_SEEDBYTES);
}

/* Expands F, but leaves a 16-byte gap after every 64.5 * 32 bytes.
   This makes it easier to evaluate the MQ function using 32-byte registers.
   Resulting structure of F:
   For all 128 accumulators;
       2x 128 bits = 32 bytes for the linear terms
       2x (128 x 129 / 2) bits = 2064 = 64x 32 + 16 bytes for the quadratic terms
       16 bytes zeroed out. */
void expand_F_gapped(uint64_t *F, const unsigned char *seed)
{
#if SOFIA_N != 128
    #error "Expanding F to align to 32 bytes depends on N = 128"
#endif
    int i, j;

    expand_F(F, seed);

    #define SOFIA_FPART (SOFIA_N + SOFIA_N * (SOFIA_N+1) / 2) / 4 / 8

    // TODO no. of moves can be reduced by splitting into multiple SHAKE calls
    // TODO  lcm of 64.5 * 32 and 168 (shake128 rate) is 14448 (7 x 64.5 x 32)
    // but that may not combine nicely with the multi-stream shake API
    for (i = SOFIA_M - 1; i >= 0; i--) {
        F[i*(SOFIA_FPART + 2) + SOFIA_FPART + 1] = 0;
        F[i*(SOFIA_FPART + 2) + SOFIA_FPART] = F[i*SOFIA_FPART + SOFIA_FPART - 1];
        F[i*(SOFIA_FPART + 2) + SOFIA_FPART - 1] = 0;
        for (j = SOFIA_FPART - 2; j >= 0; j--) {
            F[i*(SOFIA_FPART + 2) + j] = F[i*SOFIA_FPART + j];
        }
    }
}
