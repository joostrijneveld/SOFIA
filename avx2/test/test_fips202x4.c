#include "../fips202.h"
#include "../fips202x4.h"
#include "../randombytes.h"
#include <stdio.h>

#define LEN 8192

int test_shake128() {
    unsigned char inbuf[LEN*4];
    unsigned char outbuf[LEN];
    unsigned char outbuf_i[LEN*4];
    int i, j;

    randombytes(inbuf, 4*LEN);

    shake1284x(outbuf_i + LEN*0, outbuf_i + LEN*1, outbuf_i + LEN*2, outbuf_i + LEN*3, LEN,
               inbuf    + LEN*0, inbuf    + LEN*1, inbuf    + LEN*2, inbuf    + LEN*3, LEN);

    printf("Testing 4x parallel SHAKE128.. ");

    for (i = 0; i < 4; i++) {
        shake128(outbuf, LEN, inbuf + i * LEN, LEN);
        for (j = 0; j < LEN; j++) {
            if (outbuf_i[i * LEN + j] != outbuf[j]) {
                printf("Failed in instance %d.\n", i);
                return 1;
            }
        }
    }
    printf("Success.\n");
    return 0;
}

int test_cshake128_simple4x() {
    unsigned char inbuf[LEN];
    unsigned char outbuf[LEN];
    unsigned char outbuf_i[LEN*4];
    unsigned char i;
    int j;

    randombytes(inbuf, LEN);

    cshake128_simple4x(outbuf_i + LEN*0, outbuf_i + LEN*1, outbuf_i + LEN*2, outbuf_i + LEN*3, LEN,
                       0, 1, 2, 3, inbuf, LEN);

    printf("Testing 4x parallel CSHAKE128.. ");

    for (i = 0; i < 4; i++) {
        cshake128_simple(outbuf, LEN, &i, 1, inbuf, LEN);
        for (j = 0; j < LEN; j++) {
            if (outbuf_i[i * LEN + j] != outbuf[j]) {
                printf("Failed in instance %d.\n", i);
                return 1;
            }
        }
    }
    printf("Success.\n");
    return 0;
}

int main() {
    return test_shake128() | test_cshake128_simple4x();
}
