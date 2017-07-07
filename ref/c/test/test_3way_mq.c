#include <stdio.h>
#include <stdint.h>
#include "../randombytes.h"
#include "../mq.h"
#include "../F.h"
#include "../params.h"

#define TRIALS 100

int main()
{
    uint64_t x[3*SOFIA_N64]  __attribute__ ((aligned (32)));
    uint64_t y[3*SOFIA_N64]  __attribute__ ((aligned (32)));
    uint64_t fx[3*SOFIA_M64]  __attribute__ ((aligned (32)));
    uint64_t fx2[3*SOFIA_M64]  __attribute__ ((aligned (32)));
    uint64_t F[SOFIA_F64_GAPPED]  __attribute__ ((aligned (32)));
    uint64_t Fnoop_gapped[SOFIA_F64_GAPPED] __attribute__ ((aligned (32)));
    unsigned char *Fnoop_gapped_bytes = (unsigned char *)Fnoop_gapped;
    unsigned char seed[SOFIA_SEEDBYTES];
    int i, j;
    int c[8] = {0};

    // Set F to contain all-ones (bitsliced) to make it a no-op.
    for (i = 0; i < SOFIA_FBYTES_GAPPED / 32; i++) {
        for (j = 0; j < 16; j++) {
            Fnoop_gapped_bytes[i * 32 + j] = 0xFF;
            Fnoop_gapped_bytes[i * 32 + j + 16] = 0x0;
        }
    }

    for (i = 0; i < TRIALS; i++) {
        randombytes((unsigned char *)x, 3*SOFIA_NBYTES);
        randombytes((unsigned char *)y, 3*SOFIA_NBYTES);
        randombytes(seed, SOFIA_SEEDBYTES);

        expand_F_gapped(F, seed);

        randombytes((unsigned char *)fx, 3*SOFIA_MBYTES);
        randombytes((unsigned char *)fx2, 3*SOFIA_MBYTES);

        for (j = 0; j < 3; j++) {
            MQ(fx + j*SOFIA_M64, x + j * SOFIA_N64, F);
        }
        MQ3(fx2, x, F);
        for (j = 0; j < SOFIA_M64; j++) {
            if (fx2[j] != fx[j]) {
                c[0]++;
                break;
            }
        }
        for (j = 0; j < 3*SOFIA_M64; j++) {
            if (fx2[j] != fx[j]) {
                c[1]++;
                break;
            }
        }

        for (j = 0; j < 3; j++) {
            G(fx + j*SOFIA_M64, x + j * SOFIA_N64, y + j * SOFIA_N64, F);
        }
        G3(fx2, x, y, F);
        for (j = 0; j < SOFIA_M64; j++) {
            if (fx2[j] != fx[j]) {
                c[2]++;
                break;
            }
        }
        for (j = 0; j < 3*SOFIA_M64; j++) {
            if (fx2[j] != fx[j]) {
                c[3]++;
                break;
            }
        }

        for (j = 0; j < 3; j++) {
            MQ(fx + j*SOFIA_M64, x + j * SOFIA_N64, Fnoop_gapped);
        }
        MQ3(fx2, x, Fnoop_gapped);
        for (j = 0; j < SOFIA_M64; j++) {
            if (fx2[j] != fx[j]) {
                c[4]++;
                break;
            }
        }
        for (j = 0; j < 3*SOFIA_M64; j++) {
            if (fx2[j] != fx[j]) {
                c[5]++;
                break;
            }
        }

        for (j = 0; j < 3; j++) {
            G(fx + j*SOFIA_M64, x + j * SOFIA_N64, y + j * SOFIA_N64, Fnoop_gapped);
        }
        G3(fx2, x, y, Fnoop_gapped);
        for (j = 0; j < SOFIA_M64; j++) {
            if (fx2[j] != fx[j]) {
                c[6]++;
                break;
            }
        }
        for (j = 0; j < 3*SOFIA_M64; j++) {
            if (fx2[j] != fx[j]) {
                c[7]++;
                break;
            }
        }
    }

    printf("MQ == MQ3[0]:          ERR %3d / %d\n", c[0], TRIALS);
    printf("3x MQ == MQ3:          ERR %3d / %d\n", c[1], TRIALS);
    printf("G == G3[0]:            ERR %3d / %d\n", c[2], TRIALS);
    printf("3x G == G3:            ERR %3d / %d\n", c[3], TRIALS);
    printf("MQ == MQ3[0] (noop F): ERR %3d / %d\n", c[4], TRIALS);
    printf("3x MQ == MQ3 (noop F): ERR %3d / %d\n", c[5], TRIALS);
    printf("G == G3[0] (noop F):   ERR %3d / %d\n", c[6], TRIALS);
    printf("3x G == G3 (noop F):   ERR %3d / %d\n", c[7], TRIALS);

    return 0;
}