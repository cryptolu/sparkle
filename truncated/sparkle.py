'''
time sage sparkle.py INV:matrices/L8.txt >results/invL8.txt

real    10m43.031s
'''

from sage.all import *

from ttt import *

# argument - path to file with matrix,
# prepend with INV: for inverse of the matrix

MATRIX_FILE = sys.argv[1]
INVERSE = False
if MATRIX_FILE.startswith("INV:"):
    INVERSE = True
    MATRIX_FILE = MATRIX_FILE[4:]

m = open(MATRIX_FILE).read()

# for SPARKLE, matrices are 16x parallel
n = len(m.split()[0])
m = "".join(m.split())
m = matrix(GF(2), n, n, map(int, m))
MATRIX = matrix_kronecker(m, 16)
M = 64
T = MATRIX.nrows() // 64

assert MATRIX.nrows() == MATRIX.ncols() == T * M

if INVERSE:
    print "doing inverse..."
    MATRIX = ~MATRIX
    print "done"
    MATRIX_FILE += ".inverse"

try:
    import ast
    with open(MATRIX_FILE + ".transitions") as fd:
        transitions = ast.literal_eval(fd.read())
    assert transitions
    print "Loaded transitions from file"
except Exception as err:
    print "Error loading cached transitions:", err

    transitions = compute_ttt(MATRIX, T=T, M=M, debug=True)
    with open(MATRIX_FILE + ".transitions", "w") as fd:
        fd.write(`transitions`)

# in exponent
epsilon = 0.01
print "epsilon =", epsilon

# best trails

print "# Best truncated trails:"
format = "frac"
transitions = compute_ttt(MATRIX, T=T, M=M, debug=False, output_format=format)
prevrno = None
for rno, prob, generic, trail in best_truncated_trails(transitions, branch_size=M, transitions_format=format):
    if rno != prevrno:
        print "%d rounds:" % rno
        prevrno = rno

    # insignificant probabilities
    if math.log(prob, 2) < -64*sum(1-v for v in trail[-1]) + epsilon:
        continue

    diff = prob - generic
    if diff:
        diff = "2^%.3f" % math.log(diff, 2)
    if prob:
        prob = "2^%.3f" % math.log(prob, 2)
    if generic:
        generic = "2^%.3f" %  math.log(generic, 2)
    print "   ",
    print " -> ".join(map(str, trail)), ":",
    print prob, "vs", generic, "diff", diff
print "============================"

# best hulls
print "# Best truncated hulls:"
format = "frac"
transitions = compute_ttt(MATRIX, T=T, M=M, debug=False, output_format=format)
prevrno = None

# for one specific input vec?
# vec = [0] * T
# vec[0] = 1
# for rno, prob, generic, trail in best_truncated_hulls(vec, transitions, branch_size=M, transitions_format=format):

# for all input vecs
for rno, prob, generic, trail in best_truncated_hulls_all(transitions, branch_size=M, transitions_format=format):
    if rno >= 5:
        break
    if rno != prevrno:
        print "%d rounds:" % rno
        prevrno = rno

    if math.log(prob, 2) < -64*sum(1-v for v in trail[-1]) + epsilon:
        continue

    diff = prob - generic
    if diff:
        diff = "2^%.3f" % math.log(diff, 2)
    if prob:
        prob = "2^%.3f" % math.log(prob, 2)
    if generic:
        generic = "2^%.3f" %  math.log(generic, 2)
    print "   ",
    print trail, ":",
    print prob, "vs", generic, "diff", diff
