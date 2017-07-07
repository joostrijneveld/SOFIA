#include "crypto_stream_aes256ctr.h"
#include "params.h"
#include "prg.h"

static unsigned char nonce[CRYPTO_STREAM_AES256CTR_NONCEBYTES] = {0};

#if CRYPTO_STREAM_AES256CTR_KEYBYTES != SOFIA_SEEDBYTES
  #error "SOFIA_SEEDBYTES needs to match CRYPTO_STREAM_AES256CTR_KEYBYTES"
#endif

void prg(unsigned char *r, unsigned long long rlen, const unsigned char key[SOFIA_SEEDBYTES])
{
    (void) crypto_stream_aes256ctr(r, rlen, nonce, key);
}
