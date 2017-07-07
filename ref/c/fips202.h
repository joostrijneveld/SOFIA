#ifndef FIPS202_H
#define FIPS202_H

#include <stdint.h>

#define SHAKE128_RATE 168
#define SHA3_256_RATE 136
#define SHA3_512_RATE 72

#define SHAKE128_STREAM_KEYBYTES 32
#define SHAKE128_STREAM_NONCEBYTES 8

void shake128_absorb(uint64_t *s, const unsigned char *input, unsigned long long inputByteLen);
void shake128_squeezeblocks(unsigned char *output, unsigned long long nblocks, uint64_t *s);
void shake128(unsigned char *output, unsigned long long outputByteLen, const unsigned char *input, unsigned long long inputByteLen);

void cshake128_simple(unsigned char *output, unsigned long long outlen, unsigned char *custom, unsigned long long customlen, const unsigned char *in, unsigned long long inlen);
void cshake128_256simple(unsigned char *output, const char *custom, const unsigned char *in, unsigned long long inlen);

int crypto_stream_shake128(unsigned char *output, unsigned long long outlen, const unsigned char *nonce, const unsigned char *key);

#endif
