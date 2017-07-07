#include <stdio.h>
#include <stdint.h>
#include "../params.h"
#include "../randombytes.h"
#include "../F.h"

int main() {
    uint64_t F[SOFIA_F64];
    uint64_t F_gapped[SOFIA_F64_GAPPED];
    unsigned char seed[SOFIA_SEEDBYTES];
    int i, j;

    randombytes(seed, SOFIA_SEEDBYTES);

    expand_F(F, seed);
    expand_F_gapped(F_gapped, seed);

    #define F_part (SOFIA_N + SOFIA_N * (SOFIA_N+1) / 2) / 4 / 8

    // Structure of F-gapped should be 65 full 32-byte lines, and one line that
    //  has only half the bits set; since it's bitsliced, this means the low
    //  64 bits in each of the 16-byte half-lines should be set.
    for (i = 0; i < SOFIA_M; i++) {
        for (j = 0; j < 4*(SOFIA_M / 2 + 1) + 1; j++) {
            if (F[i*F_part + j] != F_gapped[i*(F_part + 2) + j]) {
                printf("expand_F and expand_F_gapped differ at "
                       "quad %d / %d.\n", i*F_part + j, i*(F_part + 2) + j);
                return 1;
            }
        }
        // now skip one;
        if (F_gapped[j] != 0) {
            printf("Expected a zeroed quadword at quad %d.\n", i*F_part + j);
            return 1;
        }
        j++;
        // now expect F = F_gapped again, at an offset
        if (F[i*F_part + j - 1] != F_gapped[i*(F_part + 2) + j]) {
            printf("expand_F and expand_F_gapped differ at "
                   "quad %d / %d.\n", i*F_part + j - 1, i*(F_part + 2) + j);
            return 1;
        }
        j++;
        // and a final quadword of zeros
        if (F_gapped[j] != 0) {
            printf("Expected a zeroed quadword at quad %d.\n", i*F_part + j);
            return 1;
        }
    }

    printf("expand_F and expand_F_gapped are consistent.\n");
    return 0;
}
