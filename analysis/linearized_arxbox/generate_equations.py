from sage.all import *
from common import *

cid = int(sys.argv[1])
C = CS[cid]

print "CONST C = %08x (index %d)" % (C, list(CS).index(C))

mid = varbits

c = Vector(map(Bit, tobin(C, BITS)))

PRINT_EQ = 0

if PRINT_EQ:
    c = Vector([Bit("c%d" % i) for i in xrange(1, BITS+1)])

ONES = Vector(Bit(1) for _ in xrange(BITS))

def latex_str(eq):
    import re
    for i in xrange(8):
        s1 = str(eq[-BITS+1+i][0]).replace("^", r"\oplus")
        s1 = re.sub(r"([cxy])(\d+)", r"\1_{\2+i}", s1)
        s2 = str(eq[-BITS+1+i][1]).replace("^", r"\oplus")
        s2 = re.sub(r"([cxy])(\d+)", r"\1_{\2+i}", s2)
        print
        print r"(%s) \cdot (%s) = 0, \\" % (s1, s2)
    return ""

def make_eqs(x, y, e, sub=True):
    res = []
    for i in xrange(BITS-1):
        if e[i] == e[i+1]:
            # quadratic x1
            res.append((x[i+1] + e[i+1], y[i+1] + e[i+1]))
        else:
            # linear x2
            res.append((x[i+1] + e[i+1] + 1, Bit.ONE))
            res.append((y[i+1] + e[i+1] + 1, Bit.ONE))
    return res

# carry pattern
ES_INT = 0, 0, 0, 0
# ES_INT = 0xfffffffe, 0, 0, 0
# ES_INT = 0xfefefefe, 0xfefefefe, 0xfefefefe, 0xfefefefe

ES_STR = ".".join("%08x" % e for e in ES_INT)
ES = [Vector(tobin(e, BITS)) for e in ES_INT]

print "CARRI-ES:", ES_STR

eqs = []

# forwards (addition)
l, r = mid.split()
l = l ^ c

DIV = "\n     "

eqs += make_eqs(l, r.ror(ROTS[2][0]), ES[2])
if PRINT_EQ: print "r2add", latex_str(eqs)
l = l ^ r.ror(ROTS[2][0]) ^ ES[2]
r = r ^ l.ror(ROTS[2][1])
l = l ^ c

eqs += make_eqs(l, r.ror(ROTS[3][0]), ES[3])
if PRINT_EQ: print "r3add", latex_str(eqs)
l = l ^ r.ror(ROTS[3][0]) ^ ES[3]
r = r ^ l.ror(ROTS[3][1])
l = l ^ c

# backwards (subtraction)
l, r = mid.split()

r = r ^ l.ror(ROTS[1][1])
eqs += make_eqs(l ^ ONES, r.ror(ROTS[1][0]), ES[1])
if PRINT_EQ: print "r1sub", latex_str(eqs)
l = l ^ r.ror(ROTS[1][0]) ^ ES[2]

l = l ^ c
r = r ^ l.ror(ROTS[0][1])
eqs += make_eqs(l ^ ONES, r.ror(ROTS[0][0]), ES[0])
if PRINT_EQ: print "r0sub", latex_str(eqs)
l = l ^ r.ror(ROTS[0][0])

if PRINT_EQ: quit()

monos = []
monos += map(tuple, Combinations(range(BITS * 2), 2))
monos += map(tuple, Combinations(range(BITS * 2), 1))
monos += [()]

polys = [a * b for a, b in eqs]

bit_to_row_cache = {}

def bit_to_row(b):
    if b not in bit_to_row_cache:
        row = vector(GF(2), len(monos))
        for mono in poly.anf:
            if mono == ():
                row[mono2i[()]] = 1
            elif len(mono) == 1:
                v1 = var2id[mono[0]]
                row[mono2i[(v1,)]] = 1
            elif len(mono) == 2:
                v1 = var2id[mono[0]]
                v2 = var2id[mono[1]]
                if v1 > v2: v1, v2 = v2, v1
                row[mono2i[(v1, v2)]] = 1
            else:
                assert 0, mono
        bit_to_row_cache[b] = row
    return bit_to_row_cache[b]

def row_to_bit(row):
    res = Bit(0)
    for mono, take in zip(monos, row):
        if not take:
            continue
        if len(mono) == 2:
            res += mid[mono[0]] * mid[mono[1]]
        elif len(mono) == 1:
            res += mid[mono[0]]
        else:
            res += Bit.ONE
    return res

mono2i = {mono: i for i, mono in enumerate(monos)}
m = matrix(GF(2), len(eqs), len(monos))
mpolys = [None] * len(eqs)
for y, poly in enumerate(polys):
    m[y] = bit_to_row(poly)
    mpolys[y] = poly

bin2 = binomial(2*BITS, 2)

print "bin2", bin2
print "rank", m[:,:bin2].rank()
print m.nrows(), "x", m.ncols()
print

eqs = list(eqs)

def matrix_row_basis_indices(m):
    """For efficient adding a bunch of vectors to a vector"""
    good = set(range(m.nrows()))
    target_rank = m.rank()
    ker = m.left_kernel().matrix()
    ker = ker[:,::-1].echelon_form()[:,::-1]
    for rel in ker:
        for i, take in reversed(list(enumerate(rel))):
            if take:
                assert i in good
                good.remove(i)
                break
    assert len(good) == target_rank
    return good

for ieq, (a, b) in enumerate(eqs[::]):
    m2 = matrix(GF(2), 2*BITS, len(monos))
    print "============"
    print "Equation #%d" % ieq
    print "eqa", a, "eqb", b
    print

    aspace = set()
    bspace = set()
    for eq, space in ((a, aspace), (b, bspace)):
        m2polys = []
        for y, v in enumerate(mid):
            poly = v * eq
            m2.set_row(y, bit_to_row(poly))
            m2polys.append(poly)
        test = m.stack(m2)
        testpolys = mpolys + m2polys
        for sol in test[:,:bin2].left_kernel().basis():
            res = Bit(0)
            for i, take in enumerate(sol):
                if take:
                    res += testpolys[i]
            space.add(res)

    aspace = {p for p in aspace if p}
    bspace = {p for p in bspace if p}

    m2 = []
    m2polys = []
    res = {}
    for ca in aspace:
        for cb in bspace:
            poly = ca * cb
            if poly not in res:
                row = bit_to_row(poly)
                print "new poly", poly
                print "a:", a, ", b:", b
                print "ca:", ca, ", cb:", cb
                print
                m2.append(row)
                m2polys.append(poly)
                res[poly] = (ca, cb)

    good = matrix_row_basis_indices(matrix(GF(2), m2))

    m2      = [item for i, item in enumerate(m2)      if i in good]
    m2polys = [item for i, item in enumerate(m2polys) if i in good]

    if not m2: continue
    old_mpolys = mpolys[::]

    m = m.stack(matrix(GF(2), m2))
    mpolys += m2polys
    good = matrix_row_basis_indices(m)

    m = matrix(GF(2), [item for i, item in enumerate(m) if i in good])
    mpolys = [item for i, item in enumerate(mpolys) if i in good]

    mpolys[:len(old_mpolys)] == old_mpolys
    assert mpolys[:len(old_mpolys)] == old_mpolys

    for poly in mpolys[len(old_mpolys):]:
        row = bit_to_row(poly)
        print "check", poly

        ca, cb = res[poly]
        eqs.append((ca, cb))
        print "new eq", ca, " | ", cb

    print len(eqs), ":", m.nrows(), "x", m.ncols()#, ":", m.rank(), m[:,:bin2].rank()
    print

f = open("eqs/%d_%08x_%s.txt" % (BITS, C, ES_STR), "w")
for eq in eqs:
    print eq[0], "|", eq[1]
    f.write("%s | %s\n" % (eq[0], eq[1]))
f.close()

print "total", len(eqs), "equations"
