#include <stdint.h>
#include "mq.h"
#include "params.h"

static void poly_F4_mul_red64(uint64_t *r_lo, uint64_t *r_hi,
                            const uint64_t *a_lo, const uint64_t *a_hi,
                            const uint64_t *b_lo, const uint64_t *b_hi)
{
    // input and output are arrays of 2 uint64_t sequences, low and high resp.
    // outputs reduced result of multiplication in GF4
    *r_hi = (*a_hi & *b_lo) ^ (*a_lo & *b_hi) ^ (*a_hi & *b_hi);
    *r_lo = (*a_lo & *b_lo) ^ (*a_hi & *b_hi);

    // without reduction, (*a_hi & *b_hi) is not yet included
}

static void rol(unsigned char *r, const unsigned char *a, const int delta, const int len)
{
    int i = 0;
    int blockdelta = (delta / (sizeof(unsigned char) * 8)) % len;
    int bitdelta = delta % (sizeof(unsigned char) * 8);
    unsigned char a2[len];

    for (i = 0; i < len; i++) {
        unsigned char t = a[(i + len - 1) % len] >> (sizeof(unsigned char) * 8 - bitdelta);
        a2[i] = a[i] << bitdelta;
        a2[i] ^= t;
    }

    for (i = 0; i < len; i++) {
        r[i] = a2[(i - blockdelta + len) % len];
    }
}

static void rol64(uint64_t *r, const uint64_t *a, const int delta, const int len)
{
    int i = 0;
    int blockdelta = (delta / (sizeof(uint64_t) * 8)) % len;
    int bitdelta = delta % (sizeof(uint64_t) * 8);
    uint64_t a2[len];

    for (i = 0; i < len; i++) {
        uint64_t t = a[(i + len - 1) % len] >> (sizeof(uint64_t) * 8 - bitdelta);
        a2[i] = a[i] << bitdelta;
        a2[i] ^= t;
    }

    for (i = 0; i < len; i++) {
        r[i] = a2[(i - blockdelta + len) % len];
    }
}

// these functions assume input and output in bitsliced representation

void MQ(uint64_t *fx, const uint64_t *x, const uint64_t *F)
{
    int i, j, k;
    uint64_t x_rol[SOFIA_N64];
    uint64_t x_quads[65][SOFIA_N64];
    uint64_t accums[SOFIA_M][SOFIA_N64];
    uint64_t t[SOFIA_N64];

    for (i = 0; i < SOFIA_N64; i++) {
        x_rol[i] = x[i];
    }

    for (i = 0; i < SOFIA_N / 8 / 2; i++) {
        for (j = 0; j < SOFIA_N64 / 2; j++) {
            poly_F4_mul_red64(x_quads[i] + j,
                            x_quads[i] + SOFIA_N64/2 + j,
                            x + j, x + SOFIA_N64/2 + j,
                            x_rol + j, x_rol + SOFIA_N64/2 + j);
        }
        rol64(x_rol, x_rol, 8, SOFIA_N64 / 2);
        rol64(x_rol + SOFIA_N64 / 2, x_rol + SOFIA_N64 / 2, 8, SOFIA_N64 / 2);
    }

    for (i = 1; i < 4; i++) {
        for (j = 0; j < SOFIA_NBYTES; j++) {
            rol(((unsigned char*)x_rol)+j, ((unsigned char*)x)+j, i, 1);
        }
        for (j = 0; j < SOFIA_N / 8; j++) {
            for (k = 0; k < SOFIA_N64 / 2; k++) {
                poly_F4_mul_red64(x_quads[8+(i-1)*(SOFIA_N / 8)+j] + k,
                                x_quads[8+(i-1)*(SOFIA_N / 8)+j] + SOFIA_N64/2 + k,
                                x + k, x + SOFIA_N64/2 + k,
                                x_rol + k, x_rol + SOFIA_N64/2 + k);
            }
            rol64(x_rol, x_rol, 8, SOFIA_N64 / 2);
            rol64(x_rol + SOFIA_N64 / 2, x_rol + SOFIA_N64 / 2, 8, SOFIA_N64 / 2);
        }
    }

    rol64(x_rol, x, 4, SOFIA_N64 / 2);
    rol64(x_rol + SOFIA_N64 / 2, x + SOFIA_N64 / 2, 4, SOFIA_N64 / 2);

    for (i = 0; i < SOFIA_N / 8 / 2; i++) {
        for (j = 0; j < SOFIA_N64 / 2; j++) {
            poly_F4_mul_red64(x_quads[8+3*(SOFIA_N / 8)+i] + j,
                            x_quads[8+3*(SOFIA_N / 8)+i] + SOFIA_N64/2 + j,
                            x + j, x + SOFIA_N64/2 + j,
                            x_rol + j, x_rol + SOFIA_N64/2 + j);
        }
        rol64(x_rol, x_rol, 8, SOFIA_N64 / 2);
        rol64(x_rol + SOFIA_N64 / 2, x_rol + SOFIA_N64 / 2, 8, SOFIA_N64 / 2);
    }

    for (i = 0; i < SOFIA_N64 / 4; i++) {
        x_rol[i] = x[i + SOFIA_N64 / 4];
        x_rol[i + SOFIA_N64 / 2] = x[i + SOFIA_N64 / 4 + SOFIA_N64 / 2];
        // mask out top half
        x_rol[i + SOFIA_N64 / 4] = 0;
        x_rol[i + SOFIA_N64 / 4 + SOFIA_N64 / 2] = 0;
    }

    for (j = 0; j < SOFIA_N64 / 2; j++) {
        poly_F4_mul_red64(x_quads[64] + j,
                          x_quads[64] + SOFIA_N64/2 + j,
                          x + j, x + SOFIA_N64/2 + j,
                          x_rol + j, x_rol + SOFIA_N64/2 + j);
    }

    for (i = 0; i < SOFIA_M; i++) {
        for (j = 0; j < SOFIA_N64; j++) {
            accums[i][j] = 0;
        }
    }

    for (k = 0; k < SOFIA_M; k++) {
        // add in the linear terms
        for (j = 0; j < SOFIA_N64 / 2; j++) {
            poly_F4_mul_red64(t + j,
                            t + SOFIA_N64/2 + j,
                            x + j, x + SOFIA_N64/2 + j,
                            &F[(int)(k*65.5 * SOFIA_N64)] + j,
                            &F[(int)(k*65.5 * SOFIA_N64)] + SOFIA_N64/2 + j);
        }
        for (j = 0; j < SOFIA_N64; j++) {
            accums[k][j] ^= t[j];
        }
        for (i = 0; i < 65; i++) {  // and for every series of quads
            for (j = 0; j < SOFIA_N64 / 2; j++) {
                poly_F4_mul_red64(t + j,
                                t + SOFIA_N64/2 + j,
                                x_quads[i] + j, x_quads[i] + SOFIA_N64/2 + j,
                                &F[(int)((k*65.5 + i + 1) * SOFIA_N64)] + j,
                                &F[(int)((k*65.5 + i + 1) * SOFIA_N64)] + SOFIA_N64/2 + j - (i == 64)*(SOFIA_N64/4));
            }
            for (j = 0; j < SOFIA_N64; j++) {
                accums[k][j] ^= t[j];
            }
        }
    }

    for (i = 0; i < SOFIA_M64; i++) {
        fx[i] = 0;
    }

    // for accums into fx
    for (i = 0; i < SOFIA_M; i++) {
        uint64_t f[2];
        for (k = 0; k < 2; k++) {
            f[k] = 0;
            for (j = 0; j < SOFIA_N64 / 2; j++) {
                f[k] ^= accums[i][j + k * (SOFIA_N64 / 2)];
            }
            f[k] ^= f[k] >> 32;
            f[k] ^= f[k] >> 16;
            f[k] ^= f[k] >> 8;
            f[k] ^= f[k] >> 4;
            f[k] ^= f[k] >> 2;
            f[k] ^= f[k] >> 1;
            f[k] &= 1;
        }
        fx[(i >> 6)] ^= f[0] << (i % 64);
        fx[(i >> 6) + SOFIA_M64 / 2] ^= f[1] << (i % 64);
    }
}

void G(uint64_t *gx, const uint64_t *x, const uint64_t *y, const uint64_t *F)
{
    int i, j, k;
    uint64_t x_rol[SOFIA_N64];
    uint64_t y_rol[SOFIA_N64];
    uint64_t x_quads[65][SOFIA_N64];
    uint64_t xiyj_lo;
    uint64_t xjyi_lo;
    uint64_t xiyj_hi;
    uint64_t xjyi_hi;
    uint64_t accums[SOFIA_M][SOFIA_N64];
    uint64_t t[SOFIA_N64];

    for (i = 0; i < SOFIA_N64; i++) {
        x_rol[i] = x[i];
    }
    for (i = 0; i < SOFIA_N64; i++) {
        y_rol[i] = y[i];
    }

    for (i = 0; i < SOFIA_N / 8 / 2; i++) {
        for (j = 0; j < SOFIA_N64 / 2; j++) {
            poly_F4_mul_red64(&xiyj_lo, &xiyj_hi,
                            x + j, x + SOFIA_N64/2 + j,
                            y_rol + j, y_rol + SOFIA_N64/2 + j);
            poly_F4_mul_red64(&xjyi_lo, &xjyi_hi,
                            y + j, y + SOFIA_N64/2 + j,
                            x_rol + j, x_rol + SOFIA_N64/2 + j);
            x_quads[i][j] = xiyj_lo ^ xjyi_lo;
            x_quads[i][j + SOFIA_N64/2] = xiyj_hi ^ xjyi_hi;
        }
        rol64(x_rol, x_rol, 8, SOFIA_N64 / 2);
        rol64(y_rol, y_rol, 8, SOFIA_N64 / 2);
        rol64(x_rol + SOFIA_N64 / 2, x_rol + SOFIA_N64 / 2, 8, SOFIA_N64 / 2);
        rol64(y_rol + SOFIA_N64 / 2, y_rol + SOFIA_N64 / 2, 8, SOFIA_N64 / 2);
    }

    for (i = 1; i < 4; i++) {
        for (j = 0; j < SOFIA_NBYTES; j++) {
            rol(((unsigned char*)x_rol)+j, ((unsigned char*)x)+j, i, 1);
            rol(((unsigned char*)y_rol)+j, ((unsigned char*)y)+j, i, 1);
        }
        for (j = 0; j < SOFIA_N / 8; j++) {
            for (k = 0; k < SOFIA_N64 / 2; k++) {
                poly_F4_mul_red64(&xiyj_lo, &xiyj_hi,
                                x + k, x + SOFIA_N64/2 + k,
                                y_rol + k, y_rol + SOFIA_N64/2 + k);
                poly_F4_mul_red64(&xjyi_lo, &xjyi_hi,
                                y + k, y + SOFIA_N64/2 + k,
                                x_rol + k, x_rol + SOFIA_N64/2 + k);
                x_quads[8+(i-1)*(SOFIA_N / 8)+j][k] = xiyj_lo ^ xjyi_lo;
                x_quads[8+(i-1)*(SOFIA_N / 8)+j][SOFIA_N64/2 + k] = xiyj_hi ^ xjyi_hi;
            }
            rol64(x_rol, x_rol, 8, SOFIA_N64 / 2);
            rol64(y_rol, y_rol, 8, SOFIA_N64 / 2);
            rol64(x_rol + SOFIA_N64 / 2, x_rol + SOFIA_N64 / 2, 8, SOFIA_N64 / 2);
            rol64(y_rol + SOFIA_N64 / 2, y_rol + SOFIA_N64 / 2, 8, SOFIA_N64 / 2);
        }
    }

    rol64(x_rol, x, 4, SOFIA_N64 / 2);
    rol64(y_rol, y, 4, SOFIA_N64 / 2);
    rol64(x_rol + SOFIA_N64 / 2, x + SOFIA_N64 / 2, 4, SOFIA_N64 / 2);
    rol64(y_rol + SOFIA_N64 / 2, y + SOFIA_N64 / 2, 4, SOFIA_N64 / 2);

    for (i = 0; i < SOFIA_N / 8 / 2; i++) {
        for (j = 0; j < SOFIA_N64 / 2; j++) {
            poly_F4_mul_red64(&xiyj_lo, &xiyj_hi,
                            x + j, x + SOFIA_N64/2 + j,
                            y_rol + j, y_rol + SOFIA_N64/2 + j);
            poly_F4_mul_red64(&xjyi_lo, &xjyi_hi,
                            y + j, y + SOFIA_N64/2 + j,
                            x_rol + j, x_rol + SOFIA_N64/2 + j);
            x_quads[8+3*(SOFIA_N / 8)+i][j] = xiyj_lo ^ xjyi_lo;
            x_quads[8+3*(SOFIA_N / 8)+i][SOFIA_N64/2 + j] = xiyj_hi ^ xjyi_hi;
        }
        rol64(x_rol, x_rol, 8, SOFIA_N64 / 2);
        rol64(y_rol, y_rol, 8, SOFIA_N64 / 2);
        rol64(x_rol + SOFIA_N64 / 2, x_rol + SOFIA_N64 / 2, 8, SOFIA_N64 / 2);
        rol64(y_rol + SOFIA_N64 / 2, y_rol + SOFIA_N64 / 2, 8, SOFIA_N64 / 2);
    }

    for (i = 0; i < SOFIA_N64 / 4; i++) {
        x_rol[i] = x[i + SOFIA_N64 / 4];
        y_rol[i] = y[i + SOFIA_N64 / 4];
        x_rol[i + SOFIA_N64 / 2] = x[i + SOFIA_N64 / 4 + SOFIA_N64 / 2];
        y_rol[i + SOFIA_N64 / 2] = y[i + SOFIA_N64 / 4 + SOFIA_N64 / 2];
        // mask out bottom half (can be combined with shufbytes in asm)
        x_rol[i + SOFIA_N64 / 4] = 0;
        y_rol[i + SOFIA_N64 / 4] = 0;
        x_rol[i + SOFIA_N64 / 4 + SOFIA_N64 / 2] = 0;
        y_rol[i + SOFIA_N64 / 4 + SOFIA_N64 / 2] = 0;
    }

    for (j = 0; j < SOFIA_N64 / 2; j++) {
        poly_F4_mul_red64(&xiyj_lo, &xiyj_hi,
                        x + j, x + SOFIA_N64/2 + j,
                        y_rol + j, y_rol + SOFIA_N64/2 + j);
        poly_F4_mul_red64(&xjyi_lo, &xjyi_hi,
                        y + j, y + SOFIA_N64/2 + j,
                        x_rol + j, x_rol + SOFIA_N64/2 + j);
        x_quads[64][j] = xiyj_lo ^ xjyi_lo;
        x_quads[64][j + SOFIA_N64/2] = xiyj_hi ^ xjyi_hi;
    }

    for (i = 0; i < SOFIA_M; i++) {
        for (j = 0; j < SOFIA_N64; j++) {
            accums[i][j] = 0;
        }
    }

    for (k = 0; k < SOFIA_M; k++) {
        // there are no linear terms in G; they cancel
        for (i = 0; i < 65; i++) {  // for every series of quads
            for (j = 0; j < SOFIA_N64 / 2; j++) {
                poly_F4_mul_red64(t + j,
                                t + SOFIA_N64/2 + j,
                                x_quads[i] + j, x_quads[i] + SOFIA_N64/2 + j,
                                &F[(int)((k*65.5 + i + 1) * SOFIA_N64)] + j,
                                &F[(int)((k*65.5 + i + 1) * SOFIA_N64)] + SOFIA_N64/2 + j - (i == 64)*(SOFIA_N64/4));
            }
            for (j = 0; j < SOFIA_N64; j++) {
                accums[k][j] ^= t[j];
            }
        }
    }

    for (i = 0; i < SOFIA_M64; i++) {
        gx[i] = 0;
    }

    // for accums into gx
    for (i = 0; i < SOFIA_M; i++) {
        uint64_t g[2];
        for (k = 0; k < 2; k++) {
            g[k] = 0;
            for (j = 0; j < SOFIA_N64 / 2; j++) {
                g[k] ^= accums[i][j + k * (SOFIA_N64 / 2)];
            }
            g[k] ^= g[k] >> 32;
            g[k] ^= g[k] >> 16;
            g[k] ^= g[k] >> 8;
            g[k] ^= g[k] >> 4;
            g[k] ^= g[k] >> 2;
            g[k] ^= g[k] >> 1;
            g[k] &= 1;
        }
        gx[(i >> 6)] ^= g[0] << (i % 64);
        gx[(i >> 6) + SOFIA_M64 / 2] ^= g[1] << (i % 64);
    }
}


void MQ3(uint64_t *fx, const uint64_t *x, const uint64_t *F)
{
    int i;
    for (i = 0; i < 3; i++) {
        MQ(fx + i*SOFIA_M64, x + i*SOFIA_N64, F);
    }
}

void G3(uint64_t *gx, const uint64_t *x, const uint64_t *y, const uint64_t *F)
{
    int i;
    for (i = 0; i < 3; i++) {
        G(gx + i*SOFIA_M64, x + i*SOFIA_N64, y + i*SOFIA_N64, F);
    }
}
