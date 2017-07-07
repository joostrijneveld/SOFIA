#include <stdio.h>
#include <stdint.h>
#include "../randombytes.h"
#include "../mq.h"
#include "../params.h"

#define TRIALS 100

int test_consistency()
{
    uint64_t x[SOFIA_N64]  __attribute__ ((aligned (32)));
    uint64_t y[SOFIA_N64]  __attribute__ ((aligned (32)));
    uint64_t fx[SOFIA_M64]  __attribute__ ((aligned (32)));
    uint64_t fx2[SOFIA_M64]  __attribute__ ((aligned (32)));
    uint64_t F[SOFIA_F64_GAPPED]  __attribute__ ((aligned (32)));
    int i, j;

    randombytes((unsigned char *)x, SOFIA_NBYTES);
    randombytes((unsigned char *)y, SOFIA_NBYTES);
    randombytes((unsigned char *)F, SOFIA_FBYTES_GAPPED);

    randombytes((unsigned char *)fx, SOFIA_MBYTES);

    MQ(fx, x, F);
    for (i = 0; i < TRIALS; i++) {
        randombytes((unsigned char *)fx2, SOFIA_MBYTES);
        MQ(fx2, x, F);
        for (j = 0; j < SOFIA_M64; j++) {
            if (fx2[j] != fx[j]) {
                return 1;
            }
        }
    }

    G(fx, x, y, F);
    for (i = 0; i < TRIALS; i++) {
        randombytes((unsigned char *)fx2, SOFIA_MBYTES);
        G(fx2, x, y, F);
        for (j = 0; j < SOFIA_M64; j++) {
            if (fx2[j] != fx[j]) {
                return 1;
            }
        }
    }

    return 0;
}

int main()
{
    int r = test_consistency();
    printf("Testing if MQ and G are deterministic.. ");
    printf(r ? "FAIL!" : "Success.");
    printf("\n");
    return r;
}
