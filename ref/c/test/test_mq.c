#include <stdio.h>
#include "../mq.h"
#include "../F.h"
#include "../params.h"
#include "../randombytes.h"

#define TRIALS 100

static void poly_F4_mul_red(unsigned char *r_lo, unsigned char *r_hi,
                            const unsigned char *a_lo, const unsigned char *a_hi,
                            const unsigned char *b_lo, const unsigned char *b_hi)
{
    // input and output are arrays of 2 unsigned char sequences, low and high resp.
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

static void G_quad_mono(unsigned char *gx, const unsigned char *x, const unsigned char *y)
{
    int i, j;
    unsigned char xi[2];
    unsigned char xj[2];
    unsigned char yi[2];
    unsigned char yj[2];
    unsigned char xiyj[2];
    unsigned char xjyi[2];
    unsigned char xiyj_xjyi[2];

    for (i = 0; i < SOFIA_MBYTES; i++) {
        gx[i] = 0;
    }

    for (i = 0; i < SOFIA_N; i++) {
        xi[0] = (x[i >> 3] >> (i & 0x7)) & 0x1;
        xi[1] = (x[SOFIA_NBYTES / 2 + (i >> 3)] >> (i & 0x7)) & 0x1;
        yi[0] = (y[i >> 3] >> (i & 0x7)) & 0x1;
        yi[1] = (y[SOFIA_NBYTES / 2 + (i >> 3)] >> (i & 0x7)) & 0x1;
        for (j = 0; j <= i; j++) {
            xj[0] = (x[j >> 3] >> (j & 0x7)) & 0x1;
            xj[1] = (x[SOFIA_NBYTES / 2 + (j >> 3)] >> (j & 0x7)) & 0x1;
            yj[0] = (y[j >> 3] >> (j & 0x7)) & 0x1;
            yj[1] = (y[SOFIA_NBYTES / 2 + (j >> 3)] >> (j & 0x7)) & 0x1;
            poly_F4_mul_red(&xiyj[0], &xiyj[1], &xi[0], &xi[1], &yj[0], &yj[1]);
            poly_F4_mul_red(&xjyi[0], &xjyi[1], &xj[0], &xj[1], &yi[0], &yi[1]);
            xiyj_xjyi[0] = xiyj[0] ^ xjyi[0];
            xiyj_xjyi[1] = xiyj[1] ^ xjyi[1];
            gx[0] ^= xiyj_xjyi[0];
            gx[SOFIA_MBYTES / 2] ^= xiyj_xjyi[1];
        }
    }
}

static void MQ_quad_mono(unsigned char *fx, const unsigned char *x)
{
    int i, j;
    unsigned char xi[2];
    unsigned char xj[2];
    unsigned char xixj[2];

    for (i = 0; i < SOFIA_MBYTES; i++) {
        fx[i] = 0;
    }

    for (i = 0; i < SOFIA_N; i++) {
        xi[0] = (x[i >> 3] >> (i & 0x7)) & 0x1;
        xi[1] = (x[SOFIA_NBYTES / 2 + (i >> 3)] >> (i & 0x7)) & 0x1;
        for (j = 0; j <= i; j++) {
            xj[0] = (x[j >> 3] >> (j & 0x7)) & 0x1;
            xj[1] = (x[SOFIA_NBYTES / 2 + (j >> 3)] >> (j & 0x7)) & 0x1;
            poly_F4_mul_red(&xixj[0], &xixj[1], &xi[0], &xi[1], &xj[0], &xj[1]);
            fx[0] ^= xixj[0];
            fx[SOFIA_MBYTES / 2] ^= xixj[1];
        }
        fx[0] ^= xi[0];
        fx[SOFIA_MBYTES / 2] ^= xi[1];
    }
}

static void MQ_quad_mono_rol(unsigned char *fx, const unsigned char *x)
{
    int i, j;
    unsigned char x_rol[SOFIA_NBYTES];
    unsigned char xixj[2];

    for (i = 0; i < SOFIA_MBYTES; i++) {
        fx[i] = 0;
    }

    for (i = 0; i < SOFIA_N / 2; i++) {
        // rotate high and low bits separately
        rol(x_rol, x, i, SOFIA_NBYTES / 2);
        rol(x_rol + SOFIA_NBYTES / 2, x + SOFIA_NBYTES / 2, i, SOFIA_NBYTES / 2);
        for (j = 0; j < SOFIA_NBYTES / 2; j++) {
            poly_F4_mul_red(&xixj[0], &xixj[1],
                            &x[j], &x[j + SOFIA_NBYTES / 2],
                            &x_rol[j], &x_rol[j + SOFIA_NBYTES / 2]);
            fx[0] ^= xixj[0];
            fx[SOFIA_MBYTES / 2] ^= xixj[1];
        }
    }
    for (i = 0; i < SOFIA_NBYTES / 4; i++) {
        poly_F4_mul_red(&xixj[0], &xixj[1],
                        &x[i], &x[i + SOFIA_NBYTES / 2],
                        &x[i + SOFIA_NBYTES / 4], &x[i + SOFIA_NBYTES / 2 + SOFIA_NBYTES / 4]);
        fx[0] ^= xixj[0];
        fx[SOFIA_MBYTES / 2] ^= xixj[1];
    }

    for (i = 0; i < SOFIA_NBYTES / 2; i++) {
        fx[0] ^= x[i];
        fx[SOFIA_MBYTES / 2] ^= x[i + SOFIA_MBYTES / 2];
    }

    // fold 8-bit fx[0] and fx[1] onto itself
    for (i = 0; i < 2; i++) {
        fx[i*(SOFIA_MBYTES / 2)] ^= fx[i*(SOFIA_MBYTES / 2)] >> 4;
        fx[i*(SOFIA_MBYTES / 2)] ^= fx[i*(SOFIA_MBYTES / 2)] >> 2;
        fx[i*(SOFIA_MBYTES / 2)] ^= fx[i*(SOFIA_MBYTES / 2)] >> 1;
        fx[i*(SOFIA_MBYTES / 2)] &= 0x01;
    }
}

static void MQ_quad_as_asm(unsigned char *fx, const unsigned char *x, const unsigned char *F)
{
    int i, j, k;
    unsigned char x_rol[SOFIA_NBYTES];
    unsigned char x_quads[65][SOFIA_NBYTES];
    unsigned char accums[SOFIA_M][SOFIA_NBYTES];
    unsigned char t[SOFIA_NBYTES];

    for (i = 0; i < SOFIA_NBYTES; i++) {
        x_rol[i] = x[i];
    }

    for (i = 0; i < SOFIA_N / 8 / 2; i++) {
        for (j = 0; j < SOFIA_NBYTES / 2; j++) {
            poly_F4_mul_red(x_quads[i] + j,
                            x_quads[i] + SOFIA_NBYTES/2 + j,
                            x + j, x + SOFIA_NBYTES/2 + j,
                            x_rol + j, x_rol + SOFIA_NBYTES/2 + j);
        }
        rol(x_rol, x_rol, 8, SOFIA_NBYTES / 2);
        rol(x_rol + SOFIA_NBYTES / 2, x_rol + SOFIA_NBYTES / 2, 8, SOFIA_NBYTES / 2);
    }

    for (i = 1; i < 4; i++) {
        for (j = 0; j < SOFIA_NBYTES; j++) {
            rol(x_rol+j, x+j, i, 1);
        }
        for (j = 0; j < SOFIA_N / 8; j++) {
            for (k = 0; k < SOFIA_NBYTES / 2; k++) {
                poly_F4_mul_red(x_quads[8+(i-1)*(SOFIA_N / 8)+j] + k,
                                x_quads[8+(i-1)*(SOFIA_N / 8)+j] + SOFIA_NBYTES/2 + k,
                                x + k, x + SOFIA_NBYTES/2 + k,
                                x_rol + k, x_rol + SOFIA_NBYTES/2 + k);
            }
            rol(x_rol, x_rol, 8, SOFIA_NBYTES / 2);
            rol(x_rol + SOFIA_NBYTES / 2, x_rol + SOFIA_NBYTES / 2, 8, SOFIA_NBYTES / 2);
        }
    }

    rol(x_rol, x, 4, SOFIA_NBYTES / 2);
    rol(x_rol + SOFIA_NBYTES / 2, x + SOFIA_NBYTES / 2, 4, SOFIA_NBYTES / 2);

    for (i = 0; i < SOFIA_N / 8 / 2; i++) {
        for (j = 0; j < SOFIA_NBYTES / 2; j++) {
            poly_F4_mul_red(x_quads[8+3*(SOFIA_N / 8)+i] + j,
                            x_quads[8+3*(SOFIA_N / 8)+i] + SOFIA_NBYTES/2 + j,
                            x + j, x + SOFIA_NBYTES/2 + j,
                            x_rol + j, x_rol + SOFIA_NBYTES/2 + j);
        }
        rol(x_rol, x_rol, 8, SOFIA_NBYTES / 2);
        rol(x_rol + SOFIA_NBYTES / 2, x_rol + SOFIA_NBYTES / 2, 8, SOFIA_NBYTES / 2);
    }

    for (i = 0; i < SOFIA_NBYTES / 4; i++) {
        x_rol[i] = x[i + SOFIA_NBYTES / 4];
        x_rol[i + SOFIA_NBYTES / 2] = x[i + SOFIA_NBYTES / 2 + SOFIA_NBYTES / 4];
        // mask out top half (can be combined with shufbytes in asm)
        x_rol[i + SOFIA_NBYTES / 4] = 0;
        x_rol[i + SOFIA_NBYTES / 4 + SOFIA_NBYTES / 2] = 0;
    }

    for (j = 0; j < SOFIA_NBYTES / 2; j++) {
        poly_F4_mul_red(x_quads[64] + j,
                        x_quads[64] + SOFIA_NBYTES/2 + j,
                        x + j, x + SOFIA_NBYTES/2 + j,
                        x_rol + j, x_rol + SOFIA_NBYTES/2 + j);
    }

    for (i = 0; i < SOFIA_M; i++) {
        for (j = 0; j < SOFIA_NBYTES; j++) {
            accums[i][j] = 0;
        }
    }

    for (k = 0; k < SOFIA_M; k++) {
        // add in the linear terms
        for (j = 0; j < SOFIA_NBYTES / 2; j++) {
            poly_F4_mul_red(t + j,
                            t + SOFIA_NBYTES/2 + j,
                            x + j, x + SOFIA_NBYTES/2 + j,
                            &F[k*66 * SOFIA_NBYTES] + j, &F[k*66 * SOFIA_NBYTES] + SOFIA_NBYTES/2 + j);
        }
        for (j = 0; j < SOFIA_NBYTES; j++) {
            accums[k][j] ^= t[j];
        }
        for (i = 0; i < 65; i++) {  // and for every series of quads
            for (j = 0; j < SOFIA_NBYTES / 2; j++) {
                poly_F4_mul_red(t + j,
                                t + SOFIA_NBYTES/2 + j,
                                x_quads[i] + j, x_quads[i] + SOFIA_NBYTES/2 + j,
                                &F[(k*66 + i + 1) * SOFIA_NBYTES] + j, &F[(k*66 + i + 1) * SOFIA_NBYTES] + SOFIA_NBYTES/2 + j);
            }
            for (j = 0; j < SOFIA_NBYTES; j++) {
                accums[k][j] ^= t[j];
            }
        }
    }

    for (i = 0; i < SOFIA_MBYTES; i++) {
        fx[i] = 0;
    }

    // for accums into fx
    for (i = 0; i < SOFIA_M; i++) {
        unsigned char f[2];
        for (k = 0; k < 2; k++) {
            f[k] = 0;
            for (j = 0; j < SOFIA_NBYTES / 2; j++) {
                f[k] ^= accums[i][j + k * (SOFIA_NBYTES / 2)];
            }
            f[k] ^= f[k] >> 4;
            f[k] ^= f[k] >> 2;
            f[k] ^= f[k] >> 1;
            f[k] &= 1;
        }
        fx[(i >> 3)] ^= f[0] << (i % 8);
        fx[(i >> 3) + SOFIA_MBYTES / 2] ^= f[1] << (i % 8);
    }
}

static void G_quad_as_asm(unsigned char *fx, const unsigned char *x,  const unsigned char *y, const unsigned char *F)
{
    int i, j, k;
    unsigned char x_rol[SOFIA_NBYTES];
    unsigned char y_rol[SOFIA_NBYTES];
    unsigned char x_quads[65][SOFIA_NBYTES];
    unsigned char xiyj_lo;
    unsigned char xjyi_lo;
    unsigned char xiyj_hi;
    unsigned char xjyi_hi;
    unsigned char accums[SOFIA_M][SOFIA_NBYTES];
    unsigned char t[SOFIA_NBYTES];

    for (i = 0; i < SOFIA_NBYTES; i++) {
        x_rol[i] = x[i];
    }
    for (i = 0; i < SOFIA_NBYTES; i++) {
        y_rol[i] = y[i];
    }

    for (i = 0; i < SOFIA_N / 8 / 2; i++) {
        for (j = 0; j < SOFIA_NBYTES / 2; j++) {
            poly_F4_mul_red(&xiyj_lo, &xiyj_hi,
                            x + j, x + SOFIA_NBYTES/2 + j,
                            y_rol + j, y_rol + SOFIA_NBYTES/2 + j);
            poly_F4_mul_red(&xjyi_lo, &xjyi_hi,
                            y + j, y + SOFIA_NBYTES/2 + j,
                            x_rol + j, x_rol + SOFIA_NBYTES/2 + j);
            x_quads[i][j] = xiyj_lo ^ xjyi_lo;
            x_quads[i][j + SOFIA_NBYTES/2] = xiyj_hi ^ xjyi_hi;
        }
        rol(x_rol, x_rol, 8, SOFIA_NBYTES / 2);
        rol(y_rol, y_rol, 8, SOFIA_NBYTES / 2);
        rol(x_rol + SOFIA_NBYTES / 2, x_rol + SOFIA_NBYTES / 2, 8, SOFIA_NBYTES / 2);
        rol(y_rol + SOFIA_NBYTES / 2, y_rol + SOFIA_NBYTES / 2, 8, SOFIA_NBYTES / 2);
    }

    for (i = 1; i < 4; i++) {
        for (j = 0; j < SOFIA_NBYTES; j++) {
            rol(x_rol+j, x+j, i, 1);
            rol(y_rol+j, y+j, i, 1);
        }
        for (j = 0; j < SOFIA_N / 8; j++) {
            for (k = 0; k < SOFIA_NBYTES / 2; k++) {
                poly_F4_mul_red(&xiyj_lo, &xiyj_hi,
                                x + k, x + SOFIA_NBYTES/2 + k,
                                y_rol + k, y_rol + SOFIA_NBYTES/2 + k);
                poly_F4_mul_red(&xjyi_lo, &xjyi_hi,
                                y + k, y + SOFIA_NBYTES/2 + k,
                                x_rol + k, x_rol + SOFIA_NBYTES/2 + k);
                x_quads[8+(i-1)*(SOFIA_N / 8)+j][k] = xiyj_lo ^ xjyi_lo;
                x_quads[8+(i-1)*(SOFIA_N / 8)+j][SOFIA_NBYTES/2 + k] = xiyj_hi ^ xjyi_hi;
            }
            rol(x_rol, x_rol, 8, SOFIA_NBYTES / 2);
            rol(y_rol, y_rol, 8, SOFIA_NBYTES / 2);
            rol(x_rol + SOFIA_NBYTES / 2, x_rol + SOFIA_NBYTES / 2, 8, SOFIA_NBYTES / 2);
            rol(y_rol + SOFIA_NBYTES / 2, y_rol + SOFIA_NBYTES / 2, 8, SOFIA_NBYTES / 2);
        }
    }

    rol(x_rol, x, 4, SOFIA_NBYTES / 2);
    rol(y_rol, y, 4, SOFIA_NBYTES / 2);
    rol(x_rol + SOFIA_NBYTES / 2, x + SOFIA_NBYTES / 2, 4, SOFIA_NBYTES / 2);
    rol(y_rol + SOFIA_NBYTES / 2, y + SOFIA_NBYTES / 2, 4, SOFIA_NBYTES / 2);

    for (i = 0; i < SOFIA_N / 8 / 2; i++) {
        for (j = 0; j < SOFIA_NBYTES / 2; j++) {
            poly_F4_mul_red(&xiyj_lo, &xiyj_hi,
                            x + j, x + SOFIA_NBYTES/2 + j,
                            y_rol + j, y_rol + SOFIA_NBYTES/2 + j);
            poly_F4_mul_red(&xjyi_lo, &xjyi_hi,
                            y + j, y + SOFIA_NBYTES/2 + j,
                            x_rol + j, x_rol + SOFIA_NBYTES/2 + j);
            x_quads[8+3*(SOFIA_N / 8)+i][j] = xiyj_lo ^ xjyi_lo;
            x_quads[8+3*(SOFIA_N / 8)+i][SOFIA_NBYTES/2 + j] = xiyj_hi ^ xjyi_hi;
        }
        rol(x_rol, x_rol, 8, SOFIA_NBYTES / 2);
        rol(y_rol, y_rol, 8, SOFIA_NBYTES / 2);
        rol(x_rol + SOFIA_NBYTES / 2, x_rol + SOFIA_NBYTES / 2, 8, SOFIA_NBYTES / 2);
        rol(y_rol + SOFIA_NBYTES / 2, y_rol + SOFIA_NBYTES / 2, 8, SOFIA_NBYTES / 2);
    }

    for (i = 0; i < SOFIA_NBYTES / 4; i++) {
        x_rol[i] = x[i + SOFIA_NBYTES / 4];
        y_rol[i] = y[i + SOFIA_NBYTES / 4];
        x_rol[i + SOFIA_NBYTES / 2] = x[i + SOFIA_NBYTES / 4 + SOFIA_NBYTES / 2];
        y_rol[i + SOFIA_NBYTES / 2] = y[i + SOFIA_NBYTES / 4 + SOFIA_NBYTES / 2];
        // mask out top half (can be combined with shufbytes in asm)
        x_rol[i + SOFIA_NBYTES / 4] = 0;
        y_rol[i + SOFIA_NBYTES / 4] = 0;
        x_rol[i + SOFIA_NBYTES / 4 + SOFIA_NBYTES / 2] = 0;
        y_rol[i + SOFIA_NBYTES / 4 + SOFIA_NBYTES / 2] = 0;
    }

    for (j = 0; j < SOFIA_NBYTES / 2; j++) {
        poly_F4_mul_red(&xiyj_lo, &xiyj_hi,
                        x + j, x + SOFIA_NBYTES/2 + j,
                        y_rol + j, y_rol + SOFIA_NBYTES/2 + j);
        poly_F4_mul_red(&xjyi_lo, &xjyi_hi,
                        y + j, y + SOFIA_NBYTES/2 + j,
                        x_rol + j, x_rol + SOFIA_NBYTES/2 + j);
        x_quads[64][j] = xiyj_lo ^ xjyi_lo;
        x_quads[64][j + SOFIA_NBYTES/2] = xiyj_hi ^ xjyi_hi;
    }

    for (i = 0; i < SOFIA_M; i++) {
        for (j = 0; j < SOFIA_NBYTES; j++) {
            accums[i][j] = 0;
        }
    }

    for (k = 0; k < SOFIA_M; k++) {
        // there are no linear terms in G; they cancel
        for (i = 0; i < 65; i++) {  // for every series of quads
            for (j = 0; j < SOFIA_NBYTES / 2; j++) {
                poly_F4_mul_red(t + j,
                                t + SOFIA_NBYTES/2 + j,
                                x_quads[i] + j, x_quads[i] + SOFIA_NBYTES/2 + j,
                                &F[(k*66 + i + 1) * SOFIA_NBYTES] + j, &F[(k*66 + i + 1) * SOFIA_NBYTES] + SOFIA_NBYTES/2 + j);
            }
            for (j = 0; j < SOFIA_NBYTES; j++) {
                accums[k][j] ^= t[j];
            }
        }
    }

    for (i = 0; i < SOFIA_MBYTES; i++) {
        fx[i] = 0;
    }

    // for accums into fx
    for (i = 0; i < SOFIA_M; i++) {
        unsigned char f[2];
        for (k = 0; k < 2; k++) {
            f[k] = 0;
            for (j = 0; j < SOFIA_NBYTES / 2; j++) {
                f[k] ^= accums[i][j + k * (SOFIA_NBYTES / 2)];
            }
            f[k] ^= f[k] >> 4;
            f[k] ^= f[k] >> 2;
            f[k] ^= f[k] >> 1;
            f[k] &= 1;
        }
        fx[(i >> 3)] ^= f[0] << (i % 8);
        fx[(i >> 3) + SOFIA_MBYTES / 2] ^= f[1] << (i % 8);
    }
}

int meta_test_polarform_mono() {
    unsigned char fx[SOFIA_MBYTES];
    unsigned char fy[SOFIA_MBYTES];
    unsigned char fxy[SOFIA_MBYTES];
    unsigned char gxy[SOFIA_MBYTES];
    unsigned char x[SOFIA_NBYTES];
    unsigned char y[SOFIA_NBYTES];
    unsigned char xy[SOFIA_NBYTES];
    int i;

    randombytes(x, SOFIA_NBYTES);
    randombytes(y, SOFIA_NBYTES);

    for (i = 0; i < SOFIA_MBYTES; i++) {
        fxy[i] = 0;
        gxy[i] = 1;
        fx[i] = 2;
        fy[i] = 3;
    }

    for (i = 0; i < SOFIA_NBYTES; i++) {
        xy[i] = x[i] ^ y[i];  // X + Y
    }

    MQ_quad_mono(fx, x);
    MQ_quad_mono(fy, y);
    MQ_quad_mono(fxy, xy);
    G_quad_mono(gxy, x, y);

    for (i = 0; i < SOFIA_MBYTES; i++) {
        fxy[i] ^= fx[i] ^ fy[i];  // F(x + y) - F(x) - F(y)
    }
    // F and G accumulate in the first (bitsliced) element
    return gxy[0] != fxy[0] || gxy[SOFIA_MBYTES / 2] != fxy[SOFIA_MBYTES / 2];
}

int meta_test_polarform_as_asm(const unsigned char *F) {
    unsigned char fx[SOFIA_MBYTES];
    unsigned char fy[SOFIA_MBYTES];
    unsigned char fxy[SOFIA_MBYTES];
    unsigned char gxy[SOFIA_MBYTES];
    unsigned char x[SOFIA_NBYTES];
    unsigned char y[SOFIA_NBYTES];
    unsigned char xy[SOFIA_NBYTES];
    int i;

    randombytes(x, SOFIA_NBYTES);
    randombytes(y, SOFIA_NBYTES);

    for (i = 0; i < SOFIA_MBYTES; i++) {
        fxy[i] = 0;
        gxy[i] = 1;
        fx[i] = 2;
        fy[i] = 3;
    }

    for (i = 0; i < SOFIA_NBYTES; i++) {
        xy[i] = x[i] ^ y[i];  // X + Y
    }

    MQ_quad_as_asm(fx, x, F);
    MQ_quad_as_asm(fy, y, F);
    MQ_quad_as_asm(fxy, xy, F);
    G_quad_as_asm(gxy, x, y, F);

    for (i = 0; i < SOFIA_MBYTES; i++) {
        fxy[i] ^= fx[i] ^ fy[i];  // F(x + y) - F(x) - F(y)
    }
    // F and G accumulate in the first (bitsliced) element
    for (i = 0; i < SOFIA_MBYTES; i++) {
        if (fxy[i] != gxy[i]) {
            return 1;
        }
    }
    return 0;
}

int main(void)
{
    int i, j;
    unsigned char fx_a[SOFIA_MBYTES] __attribute__ ((aligned (32)));
    unsigned char fx_a_noop[SOFIA_MBYTES] __attribute__ ((aligned (32)));
    unsigned char fx_b[SOFIA_MBYTES] __attribute__ ((aligned (32)));
    unsigned char fx_c[SOFIA_MBYTES] __attribute__ ((aligned (32)));
    unsigned char fx_d[SOFIA_MBYTES] __attribute__ ((aligned (32)));
    unsigned char fx_d_noop[SOFIA_MBYTES] __attribute__ ((aligned (32)));
    unsigned char fx_e[SOFIA_MBYTES] __attribute__ ((aligned (32)));
    unsigned char fx_e_noop[SOFIA_MBYTES] __attribute__ ((aligned (32)));
    unsigned char fx_f[SOFIA_MBYTES] __attribute__ ((aligned (32)));
    unsigned char fx_f_noop[SOFIA_MBYTES] __attribute__ ((aligned (32)));
    unsigned char x[SOFIA_NBYTES] __attribute__ ((aligned (32)));
    unsigned char y[SOFIA_NBYTES] __attribute__ ((aligned (32)));
    unsigned char F[SOFIA_FBYTES] __attribute__ ((aligned (32)));
#ifndef SOFIA_GAP_F
    unsigned char Fnoop[SOFIA_FBYTES] __attribute__ ((aligned (32)));
    int k;
#endif
    unsigned char F_gapped[SOFIA_FBYTES_GAPPED] __attribute__ ((aligned (32)));
    unsigned char Fnoop_gapped[SOFIA_FBYTES_GAPPED] __attribute__ ((aligned (32)));
    unsigned char Fseed[SOFIA_SEEDBYTES];

#define SOFIA_FPART ((SOFIA_N + SOFIA_N * (SOFIA_N+1) / 2) / 4)

    // Set F to contain all-ones (bitsliced) to make it a no-op.
    for (i = 0; i < SOFIA_FBYTES_GAPPED / 32; i++) {
        for (j = 0; j < 16; j++) {
            Fnoop_gapped[i * 32 + j] = 0xFF;
            Fnoop_gapped[i * 32 + j + 16] = 0x0;
        }
    }
#ifndef SOFIA_GAP_F
    // For the non-gapped F, the structure is a bit more intricate;
    //  The 65.5th line is half the length.
    for (i = 0; i < SOFIA_M; i++) {
        for (j = 0; j < 65; j++) {
            for (k = 0; k < 16; k++) {
                Fnoop[i*SOFIA_FPART + j*SOFIA_NBYTES + k] = 0xFF;
                Fnoop[i*SOFIA_FPART + j*SOFIA_NBYTES + k + 16] = 0x00;
            }
        }
        for (k = 0; k < 8; k++) {
            Fnoop[i*SOFIA_FPART + 65*SOFIA_NBYTES + k] = 0xFF;
            Fnoop[i*SOFIA_FPART + 65*SOFIA_NBYTES + k + 8] = 0x0;
        }
    }
#endif

    for (i = 0; i < SOFIA_MBYTES; i++) {
        fx_a[i] = 5;
        fx_a_noop[i] = 6;
        fx_b[i] = 2;
        fx_c[i] = 3;
        fx_d[i] = 4;
        fx_d_noop[i] = 7;
        fx_e[i] = 8;
        fx_e_noop[i] = 9;
        fx_f[i] = 10;
        fx_f_noop[i] = 11;
    }

    int c[10] = {0};
    for (i = 0; i < TRIALS; i++) {
        randombytes(x, SOFIA_NBYTES);
        randombytes(y, SOFIA_NBYTES);
        randombytes(Fseed, SOFIA_SEEDBYTES);
        expand_F_gapped((uint64_t *)F_gapped, Fseed);
        expand_F((uint64_t *)F, Fseed);
#ifdef SOFIA_GAP_F
        MQ((uint64_t *)fx_a, (uint64_t *)x, (uint64_t *)F_gapped);
        MQ((uint64_t *)fx_a_noop, (uint64_t *)x, (uint64_t *)Fnoop_gapped);
#else
        MQ((uint64_t *)fx_a, (uint64_t *)x, (uint64_t *)F);
        MQ((uint64_t *)fx_a_noop, (uint64_t *)x, (uint64_t *)Fnoop);
#endif
        MQ_quad_mono(fx_b, x);
        MQ_quad_mono_rol(fx_c, x);
        MQ_quad_as_asm(fx_d, x, F_gapped);
        MQ_quad_as_asm(fx_d_noop, x, Fnoop_gapped);
#ifdef SOFIA_GAP_F
        G((uint64_t *)fx_e, (uint64_t *)x, (uint64_t *)y, (uint64_t *)F_gapped);
        G((uint64_t *)fx_e_noop, (uint64_t *)x, (uint64_t *)y, (uint64_t *)Fnoop_gapped);
#else
        G((uint64_t *)fx_e, (uint64_t *)x, (uint64_t *)y, (uint64_t *)F);
        G((uint64_t *)fx_e_noop, (uint64_t *)x, (uint64_t *)y, (uint64_t *)Fnoop);
#endif
        G_quad_as_asm(fx_f, x, y, F_gapped);
        G_quad_as_asm(fx_f_noop, x, y, Fnoop_gapped);

        if (fx_b[0] != fx_c[0] || fx_b[SOFIA_MBYTES / 2] != fx_c[SOFIA_MBYTES / 2]) {
            c[0]++;
        }
        for (j = 0; j < SOFIA_MBYTES / 2; j++) {
            if ((fx_a_noop[j] & 1) != fx_c[0] || (fx_a_noop[j + SOFIA_MBYTES / 2] & 1) != fx_c[SOFIA_MBYTES / 2]) {
                c[1]++;
                break;
            }
        }
        c[2] += meta_test_polarform_mono();
        for (j = 0; j < SOFIA_MBYTES / 2; j++) {
            if ((fx_d_noop[j] & 1) != fx_c[0] || (fx_d_noop[j + SOFIA_MBYTES / 2] & 1) != fx_c[SOFIA_MBYTES / 2]) {
                c[3]++;
                break;
            }
        }
        for (j = 0; j < SOFIA_MBYTES; j++) {
            if (fx_d[j] != fx_a[j]) {
                c[4]++;
                break;
            }
        }
        for (j = 0; j < SOFIA_MBYTES; j++) {
            if (fx_d_noop[j] != fx_a_noop[j]) {
                c[5]++;
                break;
            }
        }
        for (j = 0; j < SOFIA_MBYTES; j++) {
            if (fx_a_noop[(SOFIA_MBYTES / 2) * (j >= SOFIA_MBYTES / 2)] != fx_a_noop[j]) {
                c[6]++;
                break;
            }
        }
        for (j = 0; j < SOFIA_MBYTES; j++) {
            if (fx_e[j] != fx_f[j]) {
                c[7]++;
                break;
            }
        }
        for (j = 0; j < SOFIA_MBYTES; j++) {
            if (fx_e_noop[j] != fx_f_noop[j]) {
                c[8]++;
                break;
            }
        }
        c[9] += meta_test_polarform_as_asm(F);
    }
    printf("MQ_quad_mono == MQ_quad_mono_rol:     ERR %3d / %d\n", c[0], TRIALS);
    printf("MQ == MQ_quad_mono_rol:               ERR %3d / %d\n", c[1], TRIALS);
    printf("meta-test MQ_quad_mono in polar form: ERR %3d / %d\n", c[2], TRIALS);
    printf("MQ_quad_as_asm == MQ_quad_mono_rol:   ERR %3d / %d\n", c[3], TRIALS);
    printf("MQ_quad_as_asm == MQ:                 ERR %3d / %d\n", c[4], TRIALS);
    printf("MQ_quad_as_asm == MQ (noop F):        ERR %3d / %d\n", c[5], TRIALS);
    printf("MQ (noop F) internally consistent:    ERR %3d / %d\n", c[6], TRIALS);
    printf("G_quad_as_asm == G:                   ERR %3d / %d\n", c[7], TRIALS);
    printf("G_quad_as_asm == G (noop F):          ERR %3d / %d\n", c[8], TRIALS);
    printf("meta-test MQ_asm in polar form:       ERR %3d / %d\n", c[9], TRIALS);

    return 0;
}
