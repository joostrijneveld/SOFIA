#include <stdio.h>
#include "../sofia.h"
#include "../randombytes.h"

#define TRIALS 100
#define MLEN 3742UL

int main(void)
{
    int c[8] = {0};
    int i;
    unsigned char pk[SOFIA_PUBLICKEYBYTES];
    unsigned char sk[SOFIA_SECRETKEYBYTES];
    unsigned char sm[MLEN + SOFIA_BYTES];
    unsigned char sm2[MLEN + SOFIA_BYTES];
    unsigned char mi[MLEN];
    unsigned char mo[MLEN + SOFIA_BYTES];
    unsigned long long smlen;
    unsigned long long mlen;
    unsigned long long idx;

    for(i = 0; i < TRIALS; i++) {
        (void) crypto_sign_keypair(pk, sk);

        for (idx = 0; idx < MLEN + SOFIA_BYTES; idx++) {
            sm[idx] = 0x37;
            mo[idx] = 0x42;
        }
        randombytes(mi, MLEN);

        (void) crypto_sign(sm, &smlen, mi, MLEN, sk);
        (void) crypto_sign(sm2, &smlen, mi, MLEN, sk);
        c[0] += smlen != SOFIA_BYTES + MLEN;
        c[1] += crypto_sign_open(mo, &mlen, sm, smlen, pk) != 0;
        c[2] += mlen != MLEN;
        if (mlen == 0 && mlen != MLEN) {
            c[3]++;  // Avoid a false positive for zeroed message.
        }
        for (idx = 0; idx < mlen; idx++) {
            if (mi[idx] != mo[idx]) {
                c[3]++;
                break;
            }
        }
        for (idx = 0; idx < smlen; idx++) {
            if (sm[i] != sm2[i]) {
                c[7]++;
                break;
            }
        }

        // Test what happens when a random bit in the signature is flipped.
        for (idx = 0; idx < MLEN + SOFIA_BYTES; idx++) {
            mo[idx] = 0x42;
        }
        randombytes((unsigned char*) &idx, sizeof(unsigned long long));
        idx %= SOFIA_BYTES*8;  // select which bit to flip
        sm[idx >> 3] ^= 1 << (idx & 0x07);
        c[4] += 1 - (crypto_sign_open(mo, &mlen, sm, smlen, pk) != 0);

        // Test with overlapping signing buffers.
        for (idx = 0; idx < MLEN; idx++) {
            sm[idx + SOFIA_BYTES] = mi[idx];
        }
        (void) crypto_sign(sm, &smlen, sm + SOFIA_BYTES, MLEN, sk);
        c[5] += crypto_sign_open(mo, &mlen, sm, smlen, pk) != 0;
        if (mlen == 0 && mlen != MLEN) {
            c[6]++;  // Avoid a false positive for zeroed message.
        }
        for (idx = 0; idx < mlen; idx++) {
            if (mi[idx] != mo[idx]) {
                c[6]++;
                break;
            }
        }
    }
    printf("Errors in signature length:                %3d/%3d\n", c[0], TRIALS);
    printf("Errors in signature:                       %3d/%3d\n", c[1], TRIALS);
    printf("Errors in message length:                  %3d/%3d\n", c[2], TRIALS);
    printf("Errors in message:                         %3d/%3d\n", c[3], TRIALS);
    printf("Errors in bitflipped signatures:           %3d/%3d\n", c[4], TRIALS);
    printf("Errors in overlapping signing - signature: %3d/%3d\n", c[5], TRIALS);
    printf("Errors in overlapping signing - message:   %3d/%3d\n", c[6], TRIALS);
    printf("Errors in if signature is deterministic:   %3d/%3d\n", c[7], TRIALS);

    return 0;
}
