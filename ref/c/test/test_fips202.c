#include <stdio.h>
#include <stdint.h>
#include "../randombytes.h"
#include "../fips202.h"

int main()
{
    unsigned char inbuf[1024];
    unsigned char outbuf0[1024];
    unsigned char outbuf1[1024];
    unsigned char block[SHAKE128_RATE];
    unsigned long long absorbed_bytes;
    uint64_t s[25];
    int i;

    shake128(outbuf0, 256, (unsigned char *)"", 0);

    printf("Testing SHAKE128 test vectors.. ");

    unsigned char v1[32] = {0x7f, 0x9c, 0x2b, 0xa4, 0xe8, 0x8f, 0x82, 0x7d,
                            0x61, 0x60, 0x45, 0x50, 0x76, 0x05, 0x85, 0x3e,
                            0xd7, 0x3b, 0x80, 0x93, 0xf6, 0xef, 0xbc, 0x88,
                            0xeb, 0x1a, 0x6e, 0xac, 0xfa, 0x66, 0xef, 0x26};

    unsigned char v2[32] = {0xf4, 0x20, 0x2e, 0x3c, 0x58, 0x52, 0xf9, 0x18,
                            0x2a, 0x04, 0x30, 0xfd, 0x81, 0x44, 0xf0, 0xa7,
                            0x4b, 0x95, 0xe7, 0x41, 0x7e, 0xca, 0xe1, 0x7d,
                            0xb0, 0xf8, 0xcf, 0xee, 0xd0, 0xe3, 0xe6, 0x6e};

    for (i = 0; i < 32; i++) {
        if (outbuf0[i] != v1[i]) {
            printf("failed v1!\n");
            break;
        }
    }
    if (i == 32) {
        shake128(outbuf0, 256, (unsigned char *)"The quick brown fox jumps over the lazy dog", 43);

        for (i = 0; i < 32; i++) {
            if (outbuf0[i] != v2[i]) {
                printf("failed v2!\n");
                break;
            }
        }
        printf("succeeded!\n");
    }

    randombytes(inbuf, 1024);
    shake128(outbuf0, 1024, inbuf, 1024);

    printf("Testing absorb / squeeze with single block.. ");

    for (i = 0; i < 25; i++) {
        s[i] = 0;
    }
    shake128_absorb(s, inbuf, 1024);
    shake128_squeezeblocks(outbuf1, 1, s);

    for (i = 0; i < SHAKE128_RATE; i++) {
        if (outbuf0[i] != outbuf1[i]) {
            printf("failed!\n");
            break;
        }
    }
    if (i == SHAKE128_RATE) {
        printf("succeeded!\n");
    }

    printf("Testing absorb / squeeze with arbitrary multiple.. ");

    for (i = 0; i < 25; i++) {
        s[i] = 0;
    }
    shake128_absorb(s, inbuf, 1024);
    shake128_squeezeblocks(outbuf1, 1024 / SHAKE128_RATE, s);
    shake128_squeezeblocks(block, 1, s);
    for (i = 0; i < 1024 % SHAKE128_RATE; i++) {
        outbuf1[1024 / SHAKE128_RATE * SHAKE128_RATE + i] = block[i];
    }

    for (i = 0; i < 1024; i++) {
        if (outbuf0[i] != outbuf1[i]) {
            printf("failed!\n");
            break;
        }
    }
    if (i == 1024) {
        printf("succeeded!\n");
    }

    printf("Testing absorb / squeeze with single partial absorb.. ");

    for (i = 0; i < 25; i++) {
        s[i] = 0;
    }
    absorbed_bytes = 0;
    shake128_partial_absorb(s, inbuf, 1024, &absorbed_bytes);
    shake128_close_absorb(s, &absorbed_bytes);
    shake128_squeezeblocks(outbuf1, 1024 / SHAKE128_RATE, s);
    shake128_squeezeblocks(block, 1, s);
    for (i = 0; i < 1024 % SHAKE128_RATE; i++) {
        outbuf1[1024 / SHAKE128_RATE * SHAKE128_RATE + i] = block[i];
    }

    for (i = 0; i < 1024; i++) {
        if (outbuf0[i] != outbuf1[i]) {
            printf("failed!\n");
            break;
        }
    }
    if (i == 1024) {
        printf("succeeded!\n");
    }

    printf("Testing absorb / squeeze with small partial absorb.. ");

    for (i = 0; i < 25; i++) {
        s[i] = 0;
    }
    absorbed_bytes = 0;
    shake128_partial_absorb(s, inbuf, 32, &absorbed_bytes);
    shake128_partial_absorb(s, inbuf + 32, 1024 - 32, &absorbed_bytes);
    shake128_close_absorb(s, &absorbed_bytes);
    shake128_squeezeblocks(outbuf1, 1024 / SHAKE128_RATE, s);
    shake128_squeezeblocks(block, 1, s);
    for (i = 0; i < 1024 % SHAKE128_RATE; i++) {
        outbuf1[1024 / SHAKE128_RATE * SHAKE128_RATE + i] = block[i];
    }

    for (i = 0; i < 1024; i++) {
        if (outbuf0[i] != outbuf1[i]) {
            printf("failed!\n");
            break;
        }
    }
    if (i == 1024) {
        printf("succeeded!\n");
    }

    printf("Testing absorb / squeeze with arbitrary partial absorb.. ");

    for (i = 0; i < 25; i++) {
        s[i] = 0;
    }
    absorbed_bytes = 0;
    shake128_partial_absorb(s, inbuf, 32, &absorbed_bytes);
    shake128_partial_absorb(s, inbuf + 32, 500, &absorbed_bytes);
    shake128_partial_absorb(s, inbuf + 532, 400, &absorbed_bytes);
    shake128_partial_absorb(s, inbuf + 932, 92, &absorbed_bytes);
    shake128_close_absorb(s, &absorbed_bytes);
    shake128_squeezeblocks(outbuf1, 1024 / SHAKE128_RATE, s);
    shake128_squeezeblocks(block, 1, s);
    for (i = 0; i < 1024 % SHAKE128_RATE; i++) {
        outbuf1[1024 / SHAKE128_RATE * SHAKE128_RATE + i] = block[i];
    }

    for (i = 0; i < 1024; i++) {
        if (outbuf0[i] != outbuf1[i]) {
            printf("failed!\n");
            break;
        }
    }
    if (i == 1024) {
        printf("succeeded!\n");
    }

    printf("Testing absorb / squeeze with single-byte partial absorb.. ");

    for (i = 0; i < 25; i++) {
        s[i] = 0;
    }
    absorbed_bytes = 0;
    for (i = 0; i < 1024; i++) {
        shake128_partial_absorb(s, inbuf + i, 1, &absorbed_bytes);
    }
    shake128_close_absorb(s, &absorbed_bytes);
    shake128_squeezeblocks(outbuf1, 1024 / SHAKE128_RATE, s);
    shake128_squeezeblocks(block, 1, s);
    for (i = 0; i < 1024 % SHAKE128_RATE; i++) {
        outbuf1[1024 / SHAKE128_RATE * SHAKE128_RATE + i] = block[i];
    }

    for (i = 0; i < 1024; i++) {
        if (outbuf0[i] != outbuf1[i]) {
            printf("failed!\n");
            break;
        }
    }
    if (i == 1024) {
        printf("succeeded!\n");
    }
}
