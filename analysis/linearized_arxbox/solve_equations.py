from sage.all import *

from common import *

try:
    # zero carry pattern
    cid = int(sys.argv[1])
    C = CS[cid]
    ES_INT = 0, 0, 0, 0
    ES_STR = ".".join("%08x" % e for e in ES_INT)
    fname = "eqs/%d_%08x_%s.txt" % (BITS, C, ES_STR)
except ValueError:
    # arbitrary carry pattern
    fname = sys.argv[1]
    info = fname
    if "/" in info:
        info = info[info.rfind("/")+1:]
    assert info.endswith(".txt")
    info = info[:-4]

    bits, const, es = info.split("_")
    assert BITS == int(bits)
    C = int(const, 16)
    ES_INT = tuple(int(e, 16) for e in es.split("."))

ES_STR = ".".join("%08x" % e for e in ES_INT)
ES = [Vector(tobin(e, BITS)) for e in ES_INT]

print "CONST C = %08x (index %d)" % (C, list(CS).index(C))
print "CARRI-ES:", ES_STR

eqs = []

def parse(s):
    res = Bit(0)
    for c in s.strip().split(" ^ "):
        if c == "1":
            res ^= Bit.ONE
        else:
            assert c[0] in "xy"
            res ^= varbits[var2id[c]]
    return res

eqs = []
for line in open(fname):
    if line.strip() == "": continue
    l, r = line.split(" | ")
    eqs.append((parse(l), parse(r)))
    print eqs[-1][0], "  *  ", eqs[-1][1], "  = 0"
print "%d total equations" % len(eqs)

order = []
sols = []

SOL_FNAME = "sol/%d_%08x_%s.txt" % (BITS, C, ES_STR)
open(SOL_FNAME, "w").close()

import time
T0 = TT = time.time()
NCHECKED = 0

def recurse(known, eqs):
    global order
    # on 2*BITS-1 the system will be linear
    assert len(known) < 2*BITS

    # 1
    if len(known) >= len(order):
        assert len(order) == len(known)
        clauses = []
        for l, r in eqs:
            lv = l.variables()
            rv = r.variables()
            if len(lv) >= 1 and len(rv) >= 1:
                clauses.append(lv)
                clauses.append(rv)
        cnt = {}
        for clause in clauses:
            for v in clause:
                if v not in cnt:
                    cnt[v] = [0] * 20
                cnt[v][len(clause)] += 1
        best = (99999, 99999), None
        for v, vcnt in cnt.items():
            score = [(i, -c) for i, c in enumerate(vcnt) if c]
            best = min(best, (score, v))
        var = best[1]
        varid = var2id[var]
        print "ORDER[%d] = %d" % (len(order), varid)
        order.append(varid)
    else:
        varid = order[len(known)]
        var = varbits[varid]

    for value in (1, 0):
        known[varid] = value
        # print "try", varid, value

        eqs2 = []
        for l, r in eqs:
            l = l.subs(var, Bit(value))
            r = r.subs(var, Bit(value))
            eqs2.append((l, r))
        res = consistent(known, eqs2)

        if res == True:
            recurse(known, eqs2)
        else:
            global NCHECKED
            NCHECKED += 2**(2*BITS - len(known))

        del known[varid]

    global TT
    if time.time() - TT > 5 and NCHECKED:
        took = time.time() - T0
        print "TIME %5.2fs" % took,
        print "CHECKED", "2^%.2f" % math.log(NCHECKED, 2),
        print "ESTIMATED LEFT", "%5.2fs" % (took * (2.0**(2*BITS) - NCHECKED) / NCHECKED)
        # order = order[:len(known)+1] # experimental, refresh order sometimes

        TT = time.time()

def consistent(known, eqscur):
    only_linear = True
    lineq = []

    # print sorted(known.items())
    for ieq, (l, r) in enumerate(eqscur):
        if l.varcount() > r.varcount():
            l, r = r, l

        if l.varcount() >= 1 and r.varcount() >= 1:
            only_linear = False
            continue
        if l.varcount() == 0 and r.varcount() == 0:
            if l.const() == r.const() == 1:
                return False
            continue
        if l.varcount() == 0 and r.varcount() >= 1:
            if l.const() == 1:
                lineq.append(r)
            continue

    mat = matrix(GF(2), len(lineq) + len(known), 2*BITS)
    target = vector(GF(2), len(lineq) + len(known))
    for y, eq in enumerate(lineq):
        for v in eq.variables():
            mat[y, var2id[v]] = 1
        target[y] = eq.const()

    y = len(lineq)
    for varid in known:
        vec = vector(GF(2), 2*BITS)
        mat[y,varid] = 1
        y += 1

    try:
        sol = mat.solve_right(target)
    except ValueError as err:
        return False
    except ZeroDivisionError as err:
        return False

    if only_linear:
        print "matrix", mat.nrows(), "x", mat.ncols(), "rank", mat.rank()
        print "solutions", "2^%d" % mat.right_kernel().dimension()
        for z in mat.right_kernel():
            solz = sol + z
            for i, v in known.items():
                solz[i] = v

            sols.append(tuple(solz))

            l = frombin(solz[:BITS])
            r = frombin(solz[BITS:])
            r = r ^ ror(l, ROTS[1][1])
            l = l ^ ror(r, ROTS[1][0]) ^ ES_INT[1]

            l = l ^ C
            r = r ^ ror(l, ROTS[0][1])
            l = l ^ ror(r, ROTS[0][0]) ^ ES_INT[0]
            print "solution: %08x:%08x" % (l, r)
            with open(SOL_FNAME, "a") as f:
                f.write("(0x%08x, 0x%08x)\n" % (l, r))
        return "linear"
    return True

recurse(known={}, eqs=eqs)

print "======================="
print "total solutions:", len(sols)
