
from collections import Counter
from itertools import combinations

N = 128  # this is assumed implicitly all over the place
pairs = []


def poly_F4_mul_red(a, b):
    a_hi = (a >> 1) & 1
    a_lo = a & 1
    b_hi = (b >> 1) & 1
    b_lo = b & 1
    r_hi = (a_hi & b_lo) ^ (a_lo & b_hi) ^ (a_hi & b_hi)
    r_lo = (a_lo & b_lo) ^ (a_hi & b_hi)
    return (r_hi << 1) ^ r_lo


def ROL(l, d, n=N):
    d = n - d
    return l[d:n] + l[:d]


def SHR(l, d, n=N):
    return l[d:n] + [None] * 4


def apply_pairs(pairs, xval):
    xval = [int(xval[i:i+2], 16) for i in range(0, len(xval), 2)]
    # convert bytes to bits
    xbits = []
    for val in xval:
        for j in range(8):
            xbits.append((val >> j) & 1)
    xlow = xbits[:N]
    xhigh = xbits[N:]
    x = [0] * N
    print(xlow)
    print(xhigh)
    for i in range(N):
        x[i] = xlow[i] + (xhigh[i] << 1)
    print(x)
    result = 0
    for a, b in pairs:
        if -1 in (a, b):
            continue
        print(x[a], x[b])
        print("r: ", poly_F4_mul_red(x[a], x[b]))
        result ^= poly_F4_mul_red(x[a], x[b])
    return ('lo', result & 1, 'hi', (result >> 1) & 1)


x = list(range(N))
x_rol = list(range(N))


def apply_shufmask(x_rol, mask):
    x = [None] * N
    for i, idx in enumerate(mask):
        x[i*8:(i+1)*8] = x_rol[idx*8:8*(idx+1)]
    return x

for j in range(N // 8 // 2):  # for each byterotation
    for k in range(N):  # for each element
        pairs.append(tuple(sorted([x[k], x_rol[k]])))
    x_rol = ROL(x_rol, 8)

x_rol = apply_shufmask(x, 2*[i for i in range(N // 8)][:N // 16])
# mask out bottom half; already had those
x_rol = ([-1] * (N // 2)) + x_rol[N // 2:]
for k in range(N):  # for each element
    pairs.append(tuple(sorted([x[k], x_rol[k]])))

for i in range(1, 4):  # for each bitrotation
    # rotate each byte by 1 bit
    x_rol = sum([ROL(x[j:j+8], i, 8) for j in range(0, len(x), 8)], [])
    for j in range(N // 8):  # for each byterotation
        for k in range(N):  # for each element
            pairs.append(tuple(sorted([x[k], x_rol[k]])))
        x_rol = ROL(x_rol, 8)

# rotate the entire thing by 4 bits
x_rol = ROL(x, 4)
for i in range(N // 8 // 2):
    for k in range(N):  # for each element
        pairs.append(tuple(sorted([x_rol[k], x[k]])))
    if i < N // 8 // 2 - 1:
        x_rol = ROL(x_rol, 8)

# xstring = "00CB96612CF7C28D5823EEB9844F1AE5B07B4611DCA7723D08D39E6934FFCA95"
# r = apply_pairs(pairs, xstring)
# print(r)
# import sys; sys.exit(0)

# We should not have any duplicates.
print(Counter(pairs))
print(len(Counter(pairs)))
dups = [k for k, v in Counter(pairs).items() if v > 1]
print('dups:', dups)
print(len(dups))

# All pairs should be generated.
allpairs = set([tuple(sorted((x, y))) for x in range(N) for y in range(N)])
pairs = set(pairs)
print('diff:', sorted(list(allpairs - pairs)))
print('diff:', [(a, b) for a, b in sorted(list(allpairs - pairs)) if (a - b) % 8 == 0])
print('diff:', [(a, b) for a, b in sorted(list(allpairs - pairs)) if (a - b) % 8 != 0])
print(len(allpairs - pairs))
