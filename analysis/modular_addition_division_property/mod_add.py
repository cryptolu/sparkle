'''
Script to verify recursive division property propagation (MILP or SMT)
of modular addition (proposed in [1]).
Compares satisfying trails to all valid trails computed
from the lookup-table of modular addition (up to 10-bit arguments).
Note that redundant trails do not affect results of the division property search
(though might affect the performance).

[1] Beierle et al. Alzette: a 64-bit ARX-box (feat. CRAX and TRAX). CRYPTO 2020, to appear.

Usage:

$ python3 -V
Python 3.8.5

$ python3 -m pip install -U bint

$ python3 mod_add.py 7
7 x 7 -> 7 bits
23487 points in size-reduced DPPT of modular addition
42627 points before size-reducing DPPT using recursive method
23487 points in size-reduced DPPT using recursive method
reduced sets match!

$ python3 mod_add.py 8
8 x 8 -> 8 bits
108388 points in size-reduced DPPT of modular addition
225948 points before size-reducing DPPT using recursive method
108388 points in size-reduced DPPT using recursive method
reduced sets match!

$ python3 mod_add.py 9
9 x 9 -> 9 bits
maximum output division weight 9
500361 points in size-reduced DPPT of modular addition
1197993 points before size-reducing DPPT using recursive method
500361 points in size-reduced DPPT using recursive method
reduced sets match!

$ python3 mod_add.py 10
10 x 10 -> 10 bits
maximum output division weight 10
2310261 points in size-reduced DPPT of modular addition
6352419 points before size-reducing DPPT using recursive method
2310261 points in size-reduced DPPT using recursive method
reduced sets match!
'''

import sys
from bint import Bin
from collections import defaultdict


def anf_in_place(bf, n):
    for k in range(n):
        halfstep = 1 << k
        step = 2 << k
        for i in range(0, len(bf), step):
            for j in range(0, halfstep):
                bf[i + j + halfstep] ^= bf[i + j]
    return bf

def size_reduce(kset):
    result = set()
    for x in kset:
        for y in kset:
            # x >= y (bitwise)
            if x != y and (x & y == y):
                break
        else:
            result.add(x)
    return result

def sbox_division(sbox, n, m):
    """
    Compute the reduced DPPT of n x m bit S-box
    Optimized a bit
    """
    assert 1 << n == len(sbox)
    assert max(sbox) < 1 << m

    by_k = defaultdict(set)
    # iterate over all products of coordinates
    for u in range(2**m):
        bf = [int(sbox[x] & u == u) for x in range(2**n)]
        anf = anf_in_place(bf, n)
        # save the product mask per each monomial it generates
        for k, val in enumerate(anf):
            if val:
                by_k[k].add(u)

    for k in by_k:
        by_k[k] = size_reduce(by_k[k])

    by_hw = defaultdict(list)
    for x in range(2**n):
        by_hw[Bin(x).hw()].append(x)

    # propagate info to "lower" monomials (at the input)
    # do in levels by HW
    for mask_hw in reversed(range(n + 1)):
        for mask in by_hw[mask_hw]:
            for bit in range(n):
                if mask & (1 << bit) == 0:
                    continue
                by_k[mask ^ (1 << bit)] |= by_k[mask]

    for k in by_k:
        by_k[k] = size_reduce(by_k[k])

    return tuple(sorted(by_k[k]) for k in range(2**n))


# Part 1
# modular addition sbox
# n x n -> n bits
n = int(sys.argv[1])
print(f"{n} x {n} -> {n} bits")
mod = 2**n
sbox = [(x + y) % mod for x in range(mod) for y in range(mod)]

div = sbox_division(sbox, 2*n, n)
points = set()
max_out_weight = 0
for k, div in enumerate(div):
    for out in div:
        points.add(Bin(k, 2*n).tuple + Bin(out, n).tuple)
        max_out_weight = max(max_out_weight, Bin(out).hw())

print("maximum output division weight", max_out_weight)
print(len(points), "points in size-reduced DPPT of modular addition")


# Part 2
# now compute points satisfied by using our recursive equations

# MILP inequalities
# sum(x, y, c, maj(x,y,c), xor(x,y,c), const) >= 0
maj_eq = [
    (-1, -1, -1, 2, 1, 0),
    (1, 1, 1, -2, -2, 1),
]

def satisfy(pt, eq):
    assert len(pt) + 1 == len(eq)
    res = sum(a * b for a, b in zip(pt, eq)) + eq[-1]
    return res >= 0

# SMT conditions
def SMT_satisfy(pt):
    a, b, c, cc, y = pt
    na = not a
    nb = not b
    nc = not c
    ncc = not cc
    ny = not y

    if (cc & y) and not (a & b & c):
        return False
    if (ncc & ny) and not (na & nb & nc):
        return False
    if (ncc & y) and not ((a ^ b ^ c) & (na | nc)):
        return False
    if (cc & ny) and not ((a | b | c) & (na | nb | nc)):
        return False
    return True

def recursive_collect(i, x, y, z, cin):
    """
    Find recursively all activity patterns
    that satisfy the two majority equations.
    """
    if i == n:
        if cin == 1:
            # if carry bit is active,
            # we stop the propagation
            # (since it's can not be used anywhere)
            return
        out.append(z)
        return
    # zbit: activity bit of the output bit
    # cout: activity bit of the carry bit
    for zbit, cout in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        pt = x[i], y[i], cin, cout, zbit
        check_eq = all(satisfy(pt, eq) for eq in maj_eq)
        check_smt = SMT_satisfy(pt)
        # verify that MILP and SMT models are equal
        assert check_smt == check_eq
        if check_eq:
            recursive_collect(i + 1, x, y, z + (zbit,), cin=cout)


points2 = set()
points2red = set()
for k in range(2**(2*n)):
    # k: input division property
    k = Bin(k, 2*n).tuple
    x = k[0*n:1*n][::-1]
    y = k[1*n:2*n][::-1]

    out = []
    recursive_collect(i=0, x=x, y=y, z=(), cin=0)
    kset = set()
    for z in out:
        p = x[::-1] + y[::-1] + z[::-1]
        points2.add(p)
        kset.add(Bin(z[::-1]).int)
    kset = size_reduce(kset)
    points2red.update({
        x[::-1] + y[::-1] + Bin(z, n).tuple for z in kset
    })

print(len(points2), "points before size-reducing DPPT using recursive method")
print(len(points2red), "points in size-reduced DPPT using recursive method")

assert points == points2red
print("reduced sets match!")
