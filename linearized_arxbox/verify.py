#-*- coding:utf-8 -*-

import ast, glob, sys
from common import *

if len(sys.argv) == 1:
    fnames = glob.glob("sol/*.txt")
else:
    fnames = sys.argv[1:]

for fname in fnames:
    info = fname
    if "/" in info:
        info = info[info.rfind("/")+1:]
    assert info.endswith(".txt")
    info = info[:-4]

    try:
        bits, const, es = info.split("_")
    except ValueError:
        bits, const = info.split("_")
        es = "0:0:0:0"

    if BITS != int(bits):
        continue

    C = int(const, 16)
    ES_INT = tuple(int(e, 16) for e in es.split(":"))
    ES_STR = ":".join("%08x" % e for e in ES_INT)
    ES = [Vector(tobin(e, BITS)) for e in ES_INT]

    try:
        sols = open(fname).read().strip().splitlines()
        sols = map(ast.literal_eval, sols)
    except IOError:
        continue

    print "CONST C = %08x (index %d)" % (C, list(CS).index(C))
    print "CARRI-ES:", ES_STR

    print "%d solutions, %d unique" % (len(sols), len(set(sols)))
    for l, r in sorted(set(sols)):
        print "INPUT %08x %08x ->" % (l, r),
        x, y = l, r
        for i in xrange(4):
            res1 = round(x, y, ROTS[i], C)
            res2 = xround(x, y, ROTS[i], C, carries=ES_INT[i])
            # print
            # print i, "%08x:%08x vs %08x:%08x" % (res1 + res2)
            assert res1 == res2
            x, y = res1
            if i < 3:
                print "%08x %08x ->" % (x, y),
            else:
                print "%08x %08x OUTPUT" % (x, y),
        print
    print
