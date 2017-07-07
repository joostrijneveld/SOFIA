#include <stdio.h>
#include "../sofia.h"
#include "../randombytes.h"
#include "../fips202.h"
#include "../sampling.h"

#define INDICES 10000

unsigned char get_doublebit(const unsigned char *buf, int i)
{
    return (buf[i >> 2] >> ((i & 0x03) << 1)) & 0x03;
}

unsigned char get_singlebit(const unsigned char *buf, int i)
{
    return (buf[i >> 3] >> (i & 0x07)) & 0x01;
}

int main()
{
    unsigned char seed[SOFIA_SEEDBYTES];
    unsigned char buf[10*INDICES];
    unsigned char indices[((INDICES * (2 + 1) + 7) & ~7) >> 3];
    int buckets[SOFIA_T] = {0};
    int k, i, buf_idx;
    int ralpha = 0, r2 = 0;

    randombytes(seed, SOFIA_SEEDBYTES);
    shake128(buf, 10*INDICES, seed, SOFIA_SEEDBYTES);

    printf("Testing if sampling challenges is correct..\n");
    for (k = 0; k < INDICES; k++) {
        sample_challenges(indices, k, seed, SOFIA_SEEDBYTES);

#if SOFIA_T == 3
        for (i = 0, buf_idx = 0; i < k; i++, buf_idx++) {
            while (get_doublebit(buf, buf_idx) == 0x03) {
                buf_idx++;
            }
            if (buf_idx >= 10*k) {
                printf("Test inconclusive. Try again.\n");
                return 1;
            }
            if (get_doublebit(buf, buf_idx) != get_doublebit(indices, i)) {
                if (ralpha == 0) {
                    ralpha = k;
                }
            }
        }
#elif SOFIA_T == 4
        for (i = 0, buf_idx = 0; i < k; i++, buf_idx++) {
            if (get_doublebit(buf, buf_idx) != get_doublebit(indices, i)) {
                if (ralpha == 0) {
                    ralpha = k;
                }
            }
        }
#endif
        for (i = 0; i < k; i++) {
            if (get_singlebit(buf, 2*buf_idx + i) != get_singlebit(indices, 2*k + i)) {
                if (r2 == 0) {
                    r2 = k;
                }
            }
        }
    }
    printf("  - alpha challenges.. ");
    if (ralpha == 0) {
        printf("Success.\n");
    }
    else {
        printf("FAIL (at %d indices)!\n", ralpha);
    }

    printf("  - binary challenges.. ");
    if (r2 == 0) {
        printf("Success.\n");
    }
    else {
        printf("FAIL (at %d indices)!\n", r2);
    }

    printf("Testing sampling distribution (should be T times ~1/T).. ");
    for (i = 0; i < INDICES; i++) {
        buckets[get_doublebit(indices, i)]++;
    }
    for (i = 0; i < SOFIA_T; i++) {
        printf("[%d] ", buckets[i]);
    }
    printf("\n");

    buckets[0] = 0;
    buckets[1] = 0;

    printf("Testing sampling distribution (should be 2 times ~1/2).. ");
    for (i = 0; i < INDICES; i++) {
        buckets[get_singlebit(indices, i + INDICES*2)]++;
    }
    for (i = 0; i < 2; i++) {
        printf("[%d] ", buckets[i]);
    }
    printf("\n");

    return ralpha | r2;
}
