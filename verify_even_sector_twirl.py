#!/usr/bin/env python3
"""Exact finite-dimensional check of the manuscript's even-sector twirling corollary.

This checks the dephasing and GL(3,2) permutation algebra on all 64
matrix-unit basis elements. It does *not* independently construct the
underlying Majorana braid circuit; that is supplied by the analytic
Gaussian-normalizer argument in the manuscript.
Python >=3.10, standard library only. Assertions are deliberately avoided.
"""
from fractions import Fraction as Q
from itertools import product


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def dot(x, y):
    return (x & y).bit_count() & 1


def matvec(bits, x):
    xx = ((x >> 2) & 1, (x >> 1) & 1, x & 1)
    return sum((sum(bits[3*i+j]*xx[j] for j in range(3)) & 1) << (2-i)
               for i in range(3))


def gl3():
    return [b for b in product((0, 1), repeat=9)
            if len({matvec(b, x) for x in range(8)}) == 8]


def dephased_matrix_unit(r, s):
    """Dephasing coefficient of |r><s| under all eight stabilizer products."""
    return Q(sum((-1)**(dot(v, r) ^ dot(v, s)) for v in range(8)), 8)


def main():
    GL = gl3()
    require(len(GL) == 168, "GL(3,2) has wrong cardinality")
    # Matrix units E_rs form a basis of all 8x8 complex matrices. Test the
    # combined map on each basis unit, not just on one convenient state.
    for r in range(8):
        for s in range(8):
            coef = dephased_matrix_unit(r, s)
            require(coef == Q(int(r==s)), f"dephasing failed E[{r},{s}]")
            output = [Q(0)]*8
            if coef:
                for M in GL:
                    t = matvec(M, r)
                    output[t] += coef / len(GL)
            target = ([Q(1)]+[Q(0)]*7) if (r==s==0) else (
                [Q(0)]+[Q(1, 7)]*7 if r==s else [Q(0)]*8)
            require(output == target, f"twirl failed for E[{r},{s}]")
    # A concrete non-isotropic diagonal distribution.
    # Linearity and the matrix-unit test imply this identity for every
    # even-sector density matrix, including complex-coherent inputs.
    diagonal = [Q(31, 100), Q(12,100), Q(10,100), Q(9,100),
                Q(8,100), Q(13,100), Q(11,100), Q(6,100)]
    require(sum(diagonal) == 1, "invalid illustrative distribution")
    output = [Q(0)]*8
    for r, value in enumerate(diagonal):
        for M in GL:
            output[matvec(M, r)] += value / len(GL)
    require(output == [diagonal[0]]+[Q(1-diagonal[0],7)]*7,
            "concrete twirl does not reproduce isotropic family")
    # Also test a bona fide pure state with nonzero off-diagonal elements.
    # Its dephased/twirled target population is |<0|psi>|^2.
    amplitudes = [Q(3), Q(1), Q(-2), Q(0), Q(1), Q(1), Q(0), Q(-1)]
    norm = sum(a*a for a in amplitudes)
    pure_out = [Q(0)]*8
    for r in range(8):
        for s in range(8):
            coeff = amplitudes[r]*amplitudes[s]/norm
            d = dephased_matrix_unit(r, s)
            if d:
                for M in GL:
                    pure_out[matvec(M, r)] += coeff*d/len(GL)
    F = amplitudes[0]*amplitudes[0]/norm
    require(pure_out == [F]+[Q(1-F,7)]*7,
            "coherent even-sector pure-state twirl failed")
    # From the finite-dimensional algebra to distillability also uses the
    # physical Gaussian implementation of the symmetry and H12 certificate.
    print("[ok] 64 matrix units; 8-element syndrome dephasing; 168 GL frame maps")
    print("[ok] target population unchanged; nonzero syndromes averaged uniformly")
    print("ALL EVEN-SECTOR TWIRL CHECKS PASSED")


if __name__ == '__main__':
    main()
