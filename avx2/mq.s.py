import math

p = print
SOFIA_N = 128

if __name__ == '__main__':
    p(".data")
    p(".align 32")

    p("no_1lsb_in_byte:")
    for i in range(32):
        p(".byte 0xFE")

    p("only_1lsb_in_byte:")
    for i in range(32):
        p(".byte 0x01")

    p("no_2lsb_in_byte:")
    for i in range(32):
        p(".byte 0xFC")

    p("only_2lsb_in_byte:")
    for i in range(32):
        p(".byte 0x03")

    p("no_3lsb_in_byte:")
    for i in range(32):
        p(".byte 0xF8")

    p("only_3lsb_in_byte:")
    for i in range(32):
        p(".byte 0x07")

    p("mask_low7bit_in_q:")
    for i in range(4):
        p(".byte 0x3F")
        for j in range(7):
            p(".byte 0")

    p("ones_quad:")
    for i in range(4):
        p(".byte 0x1")
        for j in range(7):
            p(".byte 0")

    p("rol8_shufb:")
    for _ in range(2):
        for i in range(16):
            p(".byte", (i - 1) % 16)

    p("shr64_shufb:")
    for _ in range(2):
        for i in range(8):
            p(".byte", 8+i)
        for i in range(8):
            p(".byte 255")

    p("top_128bits:")
    for i in range(16):
        p(".byte 0")
    for i in range(16):
        p(".byte 255")

    p(".text")
    p(".global MQ")
    p(".att_syntax prefix")

    p("MQ:")

    p("mov %rsp, %r8")  # Use r8 to store the old stack pointer during execution.
    p("andq $-32, %rsp")  # Align rsp to the next 32-byte value, for vmovdqa.
    # allocate for 64 quadratics, and the half-row (32 for alignment)
    p("subq ${}, %rsp".format(32 * 65))

    # rdi: fx, rsi: x, rdx: F
    x = 0
    x_rol = 1
    xixj = 2
    x_hihi = 3
    x_loxhi_lo = 4

    t0 = 15
    t1 = 14
    rol8_shufb = 13
    p("vmovdqa {}(%rsi), %ymm{}".format(0, x))
    p("vmovdqa rol8_shufb, %ymm{}".format(rol8_shufb))

    p("vpermq ${}, %ymm{}, %ymm{}".format(int('11101110', 2), x, x_hihi))

    p("vpxor %ymm{}, %ymm{}, %ymm{}".format(t1, t1, t1))
    p("vinserti128 $1, %xmm{}, %ymm{}, %ymm{}".format(x, t1, x_loxhi_lo))
    p("vpxor %ymm{}, %ymm{}, %ymm{}".format(x_loxhi_lo, x, x_loxhi_lo))

    # This variant relies on a slight precomputation of b being available;
    #  it requires a register with [b_hi / b_hi] and [b_lo ^ b_hi / b_lo].
    def poly_mul_GF4(r, a, b_hihi, b_loxhi_lo, t0):
        # *r_hi = (*a_hi & (*b_lo ^ b_hi)) ^ (*a_lo & *b_hi);
        # *r_lo = (*a_lo & *b_lo) ^          (*a_hi & *b_hi);
        p("vpand %ymm{}, %ymm{}, %ymm{}".format(a, b_loxhi_lo, r))
        p("vpermq ${}, %ymm{}, %ymm{}".format(int('01001110', 2), a, t0))
        p("vpand %ymm{}, %ymm{}, %ymm{}".format(t0, b_hihi, t0))
        p("vpxor %ymm{}, %ymm{}, %ymm{}".format(t0, r, r))

    def poly_F4_hilo_lohi_reduce(r, r_lohi, t1):
        p("vextracti128 $1, %ymm{}, %xmm{}".format(r, t1))
        # now t1 contains NULL / (*a_hi & *b_hi)
        p("vinserti128 $1, %xmm{}, %ymm{}, %ymm{}".format(r_lohi, t1, t1))
        # now t1 contains (*a_lo & *b_hi) / (*a_hi & *b_hi)
        t0 = r_lohi
        p("vpand %ymm{}, %ymm{}, %ymm{}".format(top_128bits, r_lohi, t0))
        # now t0 contains (*a_hi & *b_lo) / NULL
        p("vpxor %ymm{}, %ymm{}, %ymm{}".format(t1, r, r))
        p("vpxor %ymm{}, %ymm{}, %ymm{}".format(t0, r, r))
        # r_hi = (a_hi & b_lo) ^ (a_lo & b_hi) ^ (a_hi & b_hi);
        # r_lo = (a_lo & b_lo) ^ (a_hi & b_hi);

    def store_xixj(i, xixj):
        p("vmovdqa %ymm{}, {}(%rsp)".format(xixj, i * 32))
        # swapping after loading is actually faster than storing swapped values

    for i in range(SOFIA_N // 8 // 2):
        poly_mul_GF4(xixj, x_rol if i > 0 else x, x_hihi, x_loxhi_lo, t0)
        store_xixj(i, xixj)
        if i < SOFIA_N // 8 // 2 - 1:
            p("vpshufb %ymm{}, %ymm{}, %ymm{}".format(rol8_shufb, x_rol if i > 0 else x, x_rol))

    for i in range(1, 4):
        p("vpsllq ${}, %ymm{}, %ymm{}".format(i, x, t0))
        p("vpsrlq ${}, %ymm{}, %ymm{}".format(8-i, x, t1))
        p("vpand no_{}lsb_in_byte, %ymm{}, %ymm{}".format(i, t0, t0))
        p("vpand only_{}lsb_in_byte, %ymm{}, %ymm{}".format(i, t1, t1))
        p("vpor %ymm{}, %ymm{}, %ymm{}".format(t0, t1, t1))
        for j in range(SOFIA_N // 8):  # for each byterotation
            poly_mul_GF4(xixj, t1 if j == 0 else x_rol, x_hihi, x_loxhi_lo, t0)
            store_xixj((8+(i-1)*(SOFIA_N // 8)+j), xixj)
            if j < SOFIA_N // 8 - 1:
                p("vpshufb %ymm{}, %ymm{}, %ymm{}".format(rol8_shufb, t1 if j == 0 else x_rol, x_rol))

    p("vpsrlq ${}, %ymm{}, %ymm{}".format(60, x, t1))
    p("vpsllq ${}, %ymm{}, %ymm{}".format(4, x, t0))
    p("vpermq ${}, %ymm{}, %ymm{}".format(int('10110001', 2), t1, t1))
    p("vpor %ymm{}, %ymm{}, %ymm{}".format(t0, t1, x_rol))

    for i in range(SOFIA_N // 8 // 2):
        poly_mul_GF4(xixj, x_rol, x_hihi, x_loxhi_lo, t0)
        store_xixj(56 + i, xixj)
        if i < SOFIA_N // 8 // 2 - 1:
            p("vpshufb %ymm{}, %ymm{}, %ymm{}".format(rol8_shufb, x_rol, x_rol))

    p("vpshufb shr64_shufb, %ymm{}, %ymm{}".format(x, x_rol))
    poly_mul_GF4(xixj, x_rol, x_hihi, x_loxhi_lo, t0)
    store_xixj(64, xixj)

    # now we have produced all quadratic terms

    # The general strategy here is to accumulate products of xixj * F, but keep
    # them separated into x_hilo * f_hilo and x_lohi * f_hilo, and only combine
    # them at the very end. This saves us from having to do a costly reduction
    # for every new x that gets multiplied into the accumulators.
    # The purpose here is to keep 'innerloop' as small as possible

    xixj = 0
    xixj_swapped = 1
    result = 2
    Fval = 3
    # We can actually free up register 1 and 3 as well, by storing the result
    # on the stack and mapping xixj and xixj_swapped to the same register,
    # but it turns out that this causes more harm than it does good.
    # An extra accumulation pair turns out to be of little extra value, while
    # re-using xixj appears to make pipelining less optimal.
    accums = [(4, 5), (6, 7), (8, 9), (10, 11), (12, 13), (14, 15)]
    n_accums = len(accums)

    # We will be adding and subtracting from this pointer; if it's the rsp,
    # it may result in weird behavior when the OS assumes we free stack space.
    # In any case, Valgrind complains otherwise.
    p("mov %rsp, %r10")
    p("mov $0, %r11")
    p("vpxor %ymm{}, %ymm{}, %ymm{}".format(result, result, result))

    p("outerloop:")
    for j in range(1 + (SOFIA_N % n_accums > 0)):  # Separate 'regular' and 'final' accumulation rounds.
        # last round needs less accumulators
        if j == 1:
            accums = accums[:SOFIA_N - (math.ceil(SOFIA_N / len(accums)) - 1) * len(accums)]

        # add the linear terms
        # the accumulators are still empty, so implicit mov instead of xor
        x = xixj
        x_swapped = xixj_swapped
        p("vmovdqa {}(%rsi), %ymm{}".format(0, x))
        p("vpermq ${}, %ymm{}, %ymm{}".format(int('01001110', 2), x, x_swapped))
        for k, acc in enumerate(accums):
            # the accumulators are still empty, so no need to xor
            p("vpand {}(%rdx), %ymm{}, %ymm{}".format(32 * 66 * k, x, acc[0]))
            p("vpand {}(%rdx), %ymm{}, %ymm{}".format(32 * 66 * k, x_swapped, acc[1]))
        p("add ${}, %rdx".format(32))

        # add the quadratic terms
        # use an explicit loop to prevent overflowing L1 instruction cache
        p("mov ${}, %r9".format(65))
        p("innerloop_{}:".format(j))

        p("vmovdqa {}(%r10), %ymm{}".format(0, xixj))
        p("add $32, %r10")
        p("vpermq ${}, %ymm{}, %ymm{}".format(int('01001110', 2), xixj, xixj_swapped))
        for k, acc in enumerate(accums):
            p("vpand {}(%rdx), %ymm{}, %ymm{}".format(32 * 66 * k, xixj, Fval))
            p("vpxor %ymm{}, %ymm{}, %ymm{}".format(Fval, acc[0], acc[0]))
        for k, acc in enumerate(accums):
            p("vpand {}(%rdx), %ymm{}, %ymm{}".format(32 * 66 * k, xixj_swapped, Fval))
            p("vpxor %ymm{}, %ymm{}, %ymm{}".format(Fval, acc[1], acc[1]))
        p("add ${}, %rdx".format(32))

        p("dec %r9")
        p("jnz innerloop_{}".format(j))
        p("sub ${}, %r10".format(32 * 65))
        p("add ${}, %rdx".format((len(accums) * 66 - 66) * 32))

        top_128bits = xixj
        t0 = xixj_swapped
        p("vmovdqa top_128bits, %ymm{}".format(top_128bits))
        for acc in accums:
            poly_F4_hilo_lohi_reduce(acc[0], acc[1], t0)

        ones = accums[0][1]
        p("vmovdqa ones_quad, %ymm{}".format(ones))

        if j == 0:
            shiftreg = accums[1][1]
            t1 = accums[2][1]
            p("vmovq %r11, %xmm{}".format(shiftreg))
            p("vpbroadcastq %xmm{}, %ymm{}".format(shiftreg, shiftreg))

        for i, acc in enumerate(acc[0] for acc in accums):  # already reduced into acc[0]
            # fold acc onto itself
            p("vpermq ${}, %ymm{}, %ymm{}".format(int('10110001', 2), acc, t0))
            p("vpxor %ymm{}, %ymm{}, %ymm{}".format(t0, acc, acc))

            p("vextracti128 $1, %ymm{}, %xmm{}".format(acc, t0))
            for accpart in [acc, t0]:
                p("vmovq %xmm{}, %r9".format(accpart))
                p("popcnt %r9, %r9")
                p("vmovq %r9, %xmm{}".format(accpart))
            p("vinserti128 $1, %xmm{}, %ymm{}, %ymm{}".format(t0, acc, acc))
            p("vpand %ymm{}, %ymm{}, %ymm{}".format(acc, ones, acc))

            # position the bit properly in the results bitstring
            if j == 0:
                if i == 0:
                    p("vpand mask_low7bit_in_q, %ymm{}, %ymm{}".format(shiftreg, t1))
                else:
                    p("vpand mask_low7bit_in_q, %ymm{}, %ymm{}".format(t1, t1))
                p("cmp ${}, %r11".format(64 - i))  # delay increments
                p("jl dontpermute_{}".format(i))
                p("vpermq ${}, %ymm{}, %ymm{}".format(int('10110001', 2), acc, acc))
                p("dontpermute_{}:".format(i))
                p("vpsllvq %ymm{}, %ymm{}, %ymm{}".format(t1, acc, acc))
                p("vpaddq %ymm{}, %ymm{}, %ymm{}".format(t1, ones, t1))
            else:
                p("vpermq ${}, %ymm{}, %ymm{}".format(int('10110001', 2), acc, acc))
                p("vpsllq ${}, %ymm{}, %ymm{}".format((i + n_accums * (SOFIA_N // n_accums)) % 64, acc, acc))
            p("vpxor %ymm{}, %ymm{}, %ymm{}".format(acc, result, result))

        if j == 0:
            p("add ${}, %r11".format(len(accums)))  # do all increments at once
            p("cmp ${}, %r11".format(n_accums * (SOFIA_N // n_accums)))
            p("jl outerloop")
    p("vmovdqa %ymm{}, 0(%rdi)".format(result))

    p("mov %r8, %rsp")
    p("ret")
