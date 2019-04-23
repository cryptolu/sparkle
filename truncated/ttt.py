'''
Truncated Differentials

see Section 11.4 of the thesis
A. Udovenko, Design and Cryptanalysis of Symmetric-Key Algorithms in Black and White-box Models. Doctoral Thesis. University of Luxembourg, 2019
'''

from sage.all import *
from itertools import product
from collections import defaultdict

def frombin(v):
    '''MSB to LSB binary form'''
    return int("".join(map(str, v)), 2 )

def hw(v):
    return Integer(v).popcount()

def submatrix(mat, output_mask, input_mask):
    MX = mat.ncols() // len(input_mask)
    MY = mat.ncols() // len(output_mask)
    if min(sum(input_mask), sum(output_mask)) == 0:
        return matrix(GF(2), sum(output_mask)*MY, sum(input_mask)*MX)
    res = None
    for x, take in enumerate(input_mask):
        if not take: continue
        sub = mat[:,x*MX:(x+1)*MX]
        if res is None:
            res = sub
        else:
            res = res.augment(sub)
    mat = res
    res = None
    for y, take in enumerate(output_mask):
        if not take: continue
        sub = mat[y*MY:(y+1)*MY]
        if res is None:
            res = sub
        else:
            res = res.stack(sub)
    return res

def compute_ttt(mat, T=None, M=None, debug=False, output_format="log"):
    '''T x T matrix of M x M submatrices
    output_format: log, frac, pair
        pair is (number of inputs fitting the transition, , number of input fitting the input mask)
    '''
    if T is None: T = mat.nrows() // M
    if M is None: M = mat.nrows() // T
    assert mat.nrows() == mat.ncols() == T * M
    assert output_format in ("log", "frac", "pair")

    transitions = {}
    if debug: print >>sys.stderr, "Computing transitions of %d x %d matrix (%d x %d subdivision)" % (mat.nrows(), mat.ncols(), T, T)
    if debug: print >>sys.stderr, "I. Ranks of Submatrices & Loose Transitions"

    masks = list(product(range(2), repeat=T))

    tab_loose = []
    for input_mask in masks:
        for output_mask in masks:
            sub = submatrix(mat, [1-i for i in output_mask], input_mask)
            r =  sub.rank()
            tab_loose.append(2**(sum(input_mask)*M - r))

    if debug: print >>sys.stderr, "II. Sum over Sub-transitions"

    # undo submask sum: loose -> exact
    tab_exact = undo_submask_sum(tab_loose)

    # undo submask sum: exact dimensions by input mask
    # can be computed faster since only weight of mask matters
    mask_card = compute_mask_cards(T, branch_size=M)

    transitions = {}
    index = 0
    for input_mask in masks:
        input_mask_card = mask_card[sum(input_mask)]

        outputs = {}
        for output_mask in masks:
            exact_card = tab_exact[index]
            index += 1

            if exact_card == 0:
                continue

            assert 0 <= exact_card <= input_mask_card
            if output_format == "pair":
                outputs[output_mask] = exact_card, input_mask_card
                continue

            prob = QQ(exact_card) / input_mask_card
            assert 0 <= prob <= 1
            if output_format == "frac":
                outputs[output_mask] = prob
                continue

            sprob = "2^{%.5f} = %s" % (math.log(prob, 2), prob)

            if debug: print >>sys.stderr, "prob[0x%02x, 0x%02x] = %s" % (frombin(input_mask), frombin(output_mask), sprob)
            problog = math.log(prob, 2)
            if output_format == "log":
                outputs[output_mask] = problog
                continue

            assert False
        transitions[input_mask] = outputs
    return transitions

def undo_submask_sum(tab):
    '''Inverse of Movius-like transform over integers'''
    if len(tab) == 1:
        return tab[::]
    h = len(tab) / 2
    res0 = undo_submask_sum(tab[:h])
    subtab = [a - b for a, b in zip(tab[h:], tab[:h])]
    res1 = undo_submask_sum(subtab)
    return res0 + res1

def matrix_kronecker(mat, k):
    '''Kronecker power (expand to k parallel versions)'''
    n = mat.nrows()
    assert mat.ncols() == n

    ring = mat.base_ring()
    I = identity_matrix(ring, k)
    Z = zero_matrix(ring, k)

    matbig = matrix(ring, k*n, k*n)
    for i in xrange(n):
        for j in xrange(n):
            matbig[i*k:i*k+k, j*k:j*k+k] = I if mat[i,j] else Z
    return matbig

def best_truncated_trails(transitions, branch_size, transitions_format="log"):
    '''
    transitions: {input_mask: {output_mask: log(prob)}}
    branch_size: for distinguishing effective trails
    yields (number of rounds, prob frac or log, generic frac or log, trail)
    '''
    assert transitions_format in ("log", "frac")
    PROB1 = 1 if transitions_format == "frac" else 0
    PROB0 = 0 if transitions_format == "frac" else -Infinity
    d = {mask: PROB1 for mask in transitions if sum(mask) > 0}
    prev = {mask: (mask,) for mask in d}
    rno = 0

    presice_generic = True
    if presice_generic:
        for mask in transitions:
            T = len(mask)
            break
        mask_card = compute_mask_cards(T, branch_size=branch_size)
        precise = {}
        for wt in xrange(1, T+1):
            generic = QQ(mask_card[wt]) / (2**(T*branch_size)-1)
            if transitions_format == "log":
                generic = math.log(generic, 2)
            precise[wt] = generic
            # print >>sys.stderr, "generic: mask weight %d, prob %s" % (wt, generic)

    while True:
        rno += 1

        d2 = defaultdict(lambda: PROB0)
        prev2 = {}
        good = 0
        for input_mask in d:
            for output_mask, prob in transitions[input_mask].items():
                numzero = sum(1-v for v in output_mask)

                # warning: generic is not really precise here
                if transitions_format == "log":
                    prob2 = d[input_mask] + prob
                    generic = -numzero * branch_size
                else:
                    prob2 = d[input_mask] * prob
                    generic = QQ(2)**(-numzero * branch_size)
                if presice_generic: generic = precise[sum(output_mask)]

                if prob2 > generic:
                    good = 1
                    if prob2 > d2[output_mask]:
                        d2[output_mask] = prob2
                        prev2[output_mask] = prev[input_mask] + (output_mask,)

        for mask in sorted(d2, key=d2.__getitem__):
            if d2[mask] > PROB0:
                if transitions_format == "log":
                    generic = -numzero * branch_size
                else:
                    generic = QQ(2)**(-numzero * branch_size)
                if presice_generic: generic = precise[sum(mask)]

                yield (rno, d2[mask], generic, prev2[mask])

        if not good:
            break
        d, prev = d2, prev2

def best_truncated_hulls(input_mask, transitions, branch_size, transitions_format="log"):
    '''
    transitions: {input_mask: {output_mask: log(prob)}}
    branch_size: for distinguishing effective trails
    '''
    assert transitions_format == "frac"
    PROB1 = 1
    PROB0 = 0

    T = len(input_mask)
    mask_card = compute_mask_cards(T, branch_size=branch_size)

    precise = {0: 0}
    for wt in xrange(1, T+1):
        precise[wt] = QQ(mask_card[wt]) / (2**(T*branch_size)-1)

    m = transitions_matrix(T, transitions)
    vec = vector(QQ, 2**T)
    vec[frombin(input_mask)] = 1

    masks = list(product(range(2), repeat=T))

    rno = 0
    while True:
        rno += 1
        vec = m * vec
        assert sum(vec) == 1

        for mask_int, prob in sorted(enumerate(vec), key=lambda (i, val): val):
            mask = masks[mask_int]
            generic = precise[sum(mask)]
            if prob > generic:
                yield (rno, prob, generic, (input_mask, mask,))

def best_truncated_hulls_all(transitions, branch_size, transitions_format="log"):
    '''
    transitions: {input_mask: {output_mask: log(prob)}}
    branch_size: for distinguishing effective trails
    '''
    assert transitions_format == "frac"
    PROB1 = 1
    PROB0 = 0

    for input_mask in transitions:
        T = len(input_mask)
        break
    mask_card = compute_mask_cards(T, branch_size=branch_size)

    precise = {0: 0}
    for wt in xrange(1, T+1):
        precise[wt] = QQ(mask_card[wt]) / (2**(T*branch_size)-1)

    m1 = transitions_matrix(T, transitions)
    masks = list(product(range(2), repeat=T))

    m = identity_matrix(QQ, m1.nrows())

    rno = 0
    while True:
        rno += 1
        m = m * m1

        for input_mask_int, input_mask in enumerate(masks):
            if not input_mask_int: continue
            for mask_int, prob in sorted(enumerate(m.column(input_mask_int)), key=lambda (i, val): val):
                mask = masks[mask_int]
                generic = precise[sum(mask)]
                if prob > generic:
                    yield (rno, prob, generic, (input_mask, mask,))

def transitions_matrix(T, transitions, transitions_format="frac"):
    assert transitions_format == "frac"
    m = matrix(QQ, 2**T, 2**T)
    for input_mask in transitions:
        for output_mask in transitions[input_mask]:
            m[frombin(output_mask), frombin(input_mask)] = transitions[input_mask][output_mask]
    return m

def compute_mask_cards(branches, branch_size):
    res = [2**(i*branch_size) for i in xrange(branches+1)]
    for i in xrange(branches + 1):
        for j in xrange(i + 1, branches + 1):
            res[j] -= res[i] * binomial(j, i)
    return res
