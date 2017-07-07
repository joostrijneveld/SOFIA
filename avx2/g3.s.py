import math

p = print
SOFIA_N = 128
INSTANCES = 3

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
    p(".global G3")
    p(".att_syntax prefix")

    p("G3:")

    p("mov %rsp, %r8")  # Use r8 to store the old stack pointer during execution.
    p("andq $-32, %rsp")  # Align rsp to the next 32-byte value, for vmovdqa.
    # allocate for 64 quadratics, and the half-row (32 for alignment)
    p("subq ${}, %rsp".format(3 * 32 * 65))

    # rdi: fx, rsi: x, rdx: y, rcx: F
    x = 0
    x_rol = 1
    xixj = 2
    x_hihi = 3
    x_loxhi_lo = 4

    y = 5
    y_rol = 6
    y_hihi = 7
    y_loxhi_lo = 8

    t0 = 15
    t1 = 14
    rol8_shufb = 13
    t2 = 12
    t3 = 11
    t4 = 10
    zeroes = 9
    p("vpxor %ymm{}, %ymm{}, %ymm{}".format(zeroes, zeroes, zeroes))
    p("vmovdqa rol8_shufb, %ymm{}".format(rol8_shufb))

    for inst in range(INSTANCES):
        p("vmovdqa {}(%rsi), %ymm{}".format(32 * inst, x))
        p("vmovdqa {}(%rdx), %ymm{}".format(32 * inst, y))

        p("vpermq ${}, %ymm{}, %ymm{}".format(int('11101110', 2), x, x_hihi))
        p("vpermq ${}, %ymm{}, %ymm{}".format(int('11101110', 2), y, y_hihi))

        p("vinserti128 $1, %xmm{}, %ymm{}, %ymm{}".format(x, zeroes, x_loxhi_lo))
        p("vinserti128 $1, %xmm{}, %ymm{}, %ymm{}".format(y, zeroes, y_loxhi_lo))
        p("vpxor %ymm{}, %ymm{}, %ymm{}".format(x_loxhi_lo, x, x_loxhi_lo))
        p("vpxor %ymm{}, %ymm{}, %ymm{}".format(y_loxhi_lo, y, y_loxhi_lo))

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
            p("vmovdqa %ymm{}, {}(%rsp)".format(xixj, (inst * 65 + i) * 32))
            # swapping after loading is actually faster than storing swapped values

        for i in range(SOFIA_N // 8 // 2):
            poly_mul_GF4(t1, x_rol if i > 0 else x, y_hihi, y_loxhi_lo, t0)
            poly_mul_GF4(t2, y_rol if i > 0 else y, x_hihi, x_loxhi_lo, t0)
            p("vpxor %ymm{}, %ymm{}, %ymm{}".format(t1, t2, xixj))
            store_xixj(i, xixj)
            if i < SOFIA_N // 8 // 2 - 1:
                p("vpshufb %ymm{}, %ymm{}, %ymm{}".format(rol8_shufb, x_rol if i > 0 else x, x_rol))
                p("vpshufb %ymm{}, %ymm{}, %ymm{}".format(rol8_shufb, y_rol if i > 0 else y, y_rol))

        for i in range(1, 4):
            p("vpsllq ${}, %ymm{}, %ymm{}".format(i, x, t0))
            p("vpsrlq ${}, %ymm{}, %ymm{}".format(8-i, x, t1))
            p("vpand no_{}lsb_in_byte, %ymm{}, %ymm{}".format(i, t0, t0))
            p("vpand only_{}lsb_in_byte, %ymm{}, %ymm{}".format(i, t1, t1))
            p("vpor %ymm{}, %ymm{}, %ymm{}".format(t0, t1, t1))

            p("vpsllq ${}, %ymm{}, %ymm{}".format(i, y, t0))
            p("vpsrlq ${}, %ymm{}, %ymm{}".format(8-i, y, t2))
            p("vpand no_{}lsb_in_byte, %ymm{}, %ymm{}".format(i, t0, t0))
            p("vpand only_{}lsb_in_byte, %ymm{}, %ymm{}".format(i, t2, t2))
            p("vpor %ymm{}, %ymm{}, %ymm{}".format(t0, t2, t2))

            for j in range(SOFIA_N // 8):  # for each byterotation
                poly_mul_GF4(t3, t1 if j == 0 else x_rol, y_hihi, y_loxhi_lo, t0)
                poly_mul_GF4(t4, t2 if j == 0 else y_rol, x_hihi, x_loxhi_lo, t0)
                p("vpxor %ymm{}, %ymm{}, %ymm{}".format(t3, t4, xixj))
                store_xixj((8+(i-1)*(SOFIA_N // 8)+j), xixj)
                if j < SOFIA_N // 8 - 1:
                    p("vpshufb %ymm{}, %ymm{}, %ymm{}".format(rol8_shufb, t1 if j == 0 else x_rol, x_rol))
                    p("vpshufb %ymm{}, %ymm{}, %ymm{}".format(rol8_shufb, t2 if j == 0 else y_rol, y_rol))

        p("vpsrlq ${}, %ymm{}, %ymm{}".format(60, x, t1))
        p("vpsllq ${}, %ymm{}, %ymm{}".format(4, x, t0))
        p("vpermq ${}, %ymm{}, %ymm{}".format(int('10110001', 2), t1, t1))
        p("vpor %ymm{}, %ymm{}, %ymm{}".format(t0, t1, x_rol))

        p("vpsrlq ${}, %ymm{}, %ymm{}".format(60, y, t1))
        p("vpsllq ${}, %ymm{}, %ymm{}".format(4, y, t0))
        p("vpermq ${}, %ymm{}, %ymm{}".format(int('10110001', 2), t1, t1))
        p("vpor %ymm{}, %ymm{}, %ymm{}".format(t0, t1, y_rol))

        for i in range(SOFIA_N // 8 // 2):
            poly_mul_GF4(t1, x_rol, y_hihi, y_loxhi_lo, t0)
            poly_mul_GF4(t2, y_rol, x_hihi, x_loxhi_lo, t0)
            p("vpxor %ymm{}, %ymm{}, %ymm{}".format(t1, t2, xixj))
            store_xixj(56 + i, xixj)
            if i < SOFIA_N // 8 // 2 - 1:
                p("vpshufb %ymm{}, %ymm{}, %ymm{}".format(rol8_shufb, x_rol, x_rol))
                p("vpshufb %ymm{}, %ymm{}, %ymm{}".format(rol8_shufb, y_rol, y_rol))

        p("vpshufb shr64_shufb, %ymm{}, %ymm{}".format(x, x_rol))
        p("vpshufb shr64_shufb, %ymm{}, %ymm{}".format(y, y_rol))
        poly_mul_GF4(t1, x_rol, y_hihi, y_loxhi_lo, t0)
        poly_mul_GF4(t2, y_rol, x_hihi, x_loxhi_lo, t0)
        p("vpxor %ymm{}, %ymm{}, %ymm{}".format(t1, t2, xixj))
        store_xixj(64, xixj)

    # now we have produced all quadratic terms

    # The general strategy here is to accumulate products of xixj * F, but keep
    # them separated into x_hilo * f_hilo and x_lohi * f_hilo, and only combine
    # them at the very end. This saves us from having to do a costly reduction
    # for every new x that gets multiplied into the accumulators.
    # The purpose here is to keep 'innerloop' as small as possible

    xixj = 0
    xixj_swapped = xixj
    Fload = [1, 2]
    Fval = 3
    # As opposed to regular MQ, register pressure is much higher now.
    # We store the intermediate results in the output buffer, and map xixj
    #  and xixj_swapped to the same register.
    accums = [[(4, 5), (6, 7)], [(8, 9), (10, 11)], [(12, 13), (14, 15)]]
    n_accums = len(accums[0])

    # We will be adding and subtracting from this pointer; if it's the rsp,
    # it may result in weird behavior when the OS assumes we free stack space.
    # In any case, Valgrind complains otherwise.
    p("mov %rsp, %r10")
    p("mov $0, %r11")

    for inst in range(INSTANCES):
        p("vmovdqa %ymm{}, {}(%rdi)".format(zeroes, inst * 32))

    # only 2 accums per MQ instance; no special rounds, since 2 divides 128
    p("outerloop:")

    # no need for linear terms in G
    p("add ${}, %rcx".format(32))

    # the first iteration doesn't xor; do that outside the loop
    for k in range(n_accums):
        p("vmovdqa {}(%rcx), %ymm{}".format(32 * 66 * k, Fload[k]))
    for inst in range(INSTANCES):
        p("vmovdqa {}(%r10), %ymm{}".format(inst * 65 * 32, xixj))
        for k, acc in enumerate(accums[inst]):
            p("vpand {}(%rcx), %ymm{}, %ymm{}".format(32 * 66 * k, xixj, acc[0]))
        p("vpermq ${}, %ymm{}, %ymm{}".format(int('01001110', 2), xixj, xixj_swapped))
        for k, acc in enumerate(accums[inst]):
            p("vpand {}(%rcx), %ymm{}, %ymm{}".format(32 * 66 * k, xixj_swapped, acc[1]))
    p("add ${}, %rcx".format(32))

    # add the rest of the quadratic terms
    # use an explicit loop to prevent overflowing L1 cache
    p("mov ${}, %r9".format(64))
    p("innerloop:")
    for k in range(n_accums):
        p("vmovdqa {}(%rcx), %ymm{}".format(32 * 66 * k, Fload[k]))
    for inst in range(INSTANCES):
        p("vmovdqa {}(%r10), %ymm{}".format(inst * 65 * 32 + 32, xixj))
        for k, acc in enumerate(accums[inst]):
            p("vpand %ymm{}, %ymm{}, %ymm{}".format(Fload[k], xixj, Fval))
            p("vpxor %ymm{}, %ymm{}, %ymm{}".format(Fval, acc[0], acc[0]))
        p("vpermq ${}, %ymm{}, %ymm{}".format(int('01001110', 2), xixj, xixj_swapped))
        for k, acc in enumerate(accums[inst]):
            p("vpand %ymm{}, %ymm{}, %ymm{}".format(Fload[k], xixj_swapped, Fval))
            p("vpxor %ymm{}, %ymm{}, %ymm{}".format(Fval, acc[1], acc[1]))
    p("add ${}, %rcx".format(32))

    p("add $32, %r10")
    p("dec %r9")
    p("jnz innerloop")
    p("sub ${}, %r10".format(32 * 64))
    p("add ${}, %rcx".format((len(accums[0]) * 66 - 66) * 32))

    top_128bits = xixj
    t0 = Fval
    p("vmovdqa top_128bits, %ymm{}".format(top_128bits))
    for inst in range(INSTANCES):
        for acc in accums[inst]:
            poly_F4_hilo_lohi_reduce(acc[0], acc[1], t0)

    ones = accums[0][0][1]
    p("vmovdqa ones_quad, %ymm{}".format(ones))

    shiftreg = accums[0][1][1]
    t1 = accums[1][0][1]
    result = accums[1][1][1]

    for inst in range(INSTANCES):
        p("vmovdqa {}(%rdi), %ymm{}".format(inst * 32, result))
        p("vmovq %r11, %xmm{}".format(shiftreg))
        p("vpbroadcastq %xmm{}, %ymm{}".format(shiftreg, shiftreg))
        for i, acc in enumerate(acc[0] for acc in accums[inst]):  # already reduced into acc[0]
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
            if i == 0:
                p("vpand mask_low7bit_in_q, %ymm{}, %ymm{}".format(shiftreg, t1))
            else:
                p("vpand mask_low7bit_in_q, %ymm{}, %ymm{}".format(t1, t1))
            p("cmp ${}, %r11".format(64 - i))  # delay increments
            p("jl dontpermute_{}_{}".format(i, inst))
            p("vpermq ${}, %ymm{}, %ymm{}".format(int('10110001', 2), acc, acc))
            p("dontpermute_{}_{}:".format(i, inst))
            p("vpsllvq %ymm{}, %ymm{}, %ymm{}".format(t1, acc, acc))
            p("vpaddq %ymm{}, %ymm{}, %ymm{}".format(t1, ones, t1))
            p("vpxor %ymm{}, %ymm{}, %ymm{}".format(acc, result, result))
        p("vmovdqa %ymm{}, {}(%rdi)".format(result, inst * 32))

    p("add ${}, %r11".format(len(accums[0])))  # do all increments at once
    p("cmp ${}, %r11".format(n_accums * (SOFIA_N // n_accums)))
    p("jl outerloop")

    p("mov %r8, %rsp")
    p("ret")
