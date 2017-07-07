#!/usr/bin/python

from math import log, ceil, factorial, floor

def nCr(n,r):
    return factorial(n) // factorial(r) // factorial(n-r)

security = 128
com_len = 256
hash_len = 256

q_N_M = [(4, 128),
         (4, 130),
         (5, 112),
         (7, 93),
         (7, 96),
         (8, 87),
         (11, 78),
         (13, 75),
         (16, 72),
         (17, 71),
         (19, 70),
         (23, 67),
         (29, 65),
         (31, 64),
         (32, 64),
         (37, 63),
         (39, 63),
         (41, 62),
         (43, 62),
         (47, 61),
         (51, 61),
         (53, 60),
         (57, 60),
         (59, 60),
         (61, 60),
         (64, 60),
         (127, 56),
         (251, 54)]

print("5-pass:")

results = []

for q, M in q_N_M:
    N = M
    rt_len = N * ceil(log(q) / log(2))
    e_len = M * ceil(log(q) / log(2))

    resp1_len = rt_len + e_len
    resp2_len = rt_len

    h1_len = resp1_len
    h2_len = resp2_len

    for t in range(2, q+1):
        error = 1/2.0 + 1.0 / (2.0 * t)
        r = ceil(-(2*security)* log(2) / log(error))
        size = com_len + (t - 1)*h1_len + (2 - 1)*h2_len
        size += resp1_len + resp2_len
        size *= r

        Fq_elements = r * ((t-1)*2*N + (2-1)*M + 2*N + M)
        idealized_compressed_size = r*com_len + ceil(Fq_elements * log(q) / log(2))

        Fq_per_quadword = floor(64 / (log(q) / log(2)))
        quadword_compressed_size = r*com_len + ceil(64 * Fq_elements / Fq_per_quadword)

        # add transcript hash
        size += hash_len
        idealized_compressed_size += hash_len
        quadword_compressed_size += hash_len

        results.append((N, q, t, r, size, idealized_compressed_size, quadword_compressed_size))

results.sort(key=lambda x: x[6])

for N, q, t, r, size, compressed_size, quadword_compressed_size in results:
        print("q = {:3}, N = M = {:3}, rounds = {:3}, t = {:3} => size: = {:.2f} KiB, compsize = {:.2f} KB, 64-bit compsize = {:.2f} KB"
              .format(q, N, r, t,
                      size / 8 / 1024,
                      compressed_size / 8 / 1024,
                      quadword_compressed_size / 8 / 1024))

print("3-pass:")

N = M = 256
q = 2

rt_len = N * ceil(log(q) / log(2))
e_len = M * ceil(log(q) / log(2))

resp_len = 2*rt_len + e_len
h_len = resp_len

t = 3
error = 2 / 3.0
r = ceil(-(2*security) * log(2) / log(error))
size = com_len + (t - 1)*h_len + resp_len
size *= r
size += hash_len

print("size: {} = {:.2f} KiB"
      .format(size, size / 8 / 1024))
