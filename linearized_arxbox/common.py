#-*- coding:utf-8 -*-

BITS = 32
MASK = (1 << BITS) - 1

CS = 0xB7E15162, 0xBF715880, 0x38B4DA56, 0x324E7738, 0xBB1185EB, 0x4F7C7B57, 0xCFBFA1C8, 0xC2B3293D, 0, 0xFFffFFff

ROTS = (
(31, 24),
(17, 17),
( 0, 31),
(24, 16),
)

CS = [c & MASK for c in CS]
ROTS = [(a % BITS, b % BITS) for a, b in ROTS]


def rol(x, n, bits=BITS):
    mask = (1 << bits) - 1
    n %= bits
    x &= mask
    return ((x << n) | (x >> (bits - n))) & mask

def ror(x, n, bits=BITS):
    return rol(x, bits - n, bits)

def xor(a, b):
    return a ^ b

def tobin(x, n):
    """MSB to LSB binary form"""
    assert 0 <= x < 1<<n
    return tuple(map(int, bin(x).lstrip("0b").rjust(n, "0")))

def frombin(v):
    """MSB to LSB binary form"""
    return int("".join(map(str, v)), 2 )

def round(x, y, rots, c):
    x += ror(y, rots[0])
    x &= MASK
    y ^= ror(x, rots[1])
    x ^= c
    return x, y

def xround(x, y, rots, c, carries=0):
    x ^= ror(y, rots[0]) ^ carries
    y ^= ror(x, rots[1])
    x ^= c
    return x, y

def ARXbox(x, y, c, rounds=4):
    for rno in xrange(rounds):
        x, y = round(x, y, ROTS[rno], c)
    return x, y

def XRXbox(x, y, c, rounds=4):
    for rno in xrange(rounds):
        x, y = xround(x, y, ROTS[rno], c)
    return x, y


from symbolic import Bit
from cryptools.py.containers import Vector

varbits = Vector(
    [Bit("x%d" % i) for i in xrange(1, BITS+1)] +
    [Bit("y%d" % i) for i in xrange(1, BITS+1)]
)

var2id = {v: i for i, v in enumerate(varbits)}

