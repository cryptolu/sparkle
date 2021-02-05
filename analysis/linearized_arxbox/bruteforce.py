import sys
from common import *

def iter_solutions():
    for l in xrange(2**BITS):
        for r in xrange(2**BITS):
            inp = x, y = l, r
            for i in xrange(4):
                test1 = xround(x, y, ROTS[i], C)
                test2 = round(x, y, ROTS[i], C)
                if test1 != test2:
                    break
                x, y = test1
            else:
                yield l, r

if __name__ == '__main__':
    cid = int(sys.argv[1])
    C = CS[cid]
    print "CONST C = %08x (index %d)" % (C, list(CS).index(C))
    num = 0
    for l, r in iter_solutions():
        print "Input %08x:%08x" % (l, r)
        num += 1
    print "Total:", num
