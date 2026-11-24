#!/usr/bin/env python3
"""Exact independent integer-weight comparison of H12 and its all-S control.

The original H12 differs from H12_BAL only at stage four, where R_7 uses
an independent raw source and S_7 uses a second processed block.

Python >= 3.10; no third-party dependencies.
"""
from fractions import Fraction
from decimal import Decimal, localcontext


class VerificationError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise VerificationError(message)

E = (0, 1, 2, 1, 4, 1, 2, 1)
WORD = (1, 2, 4, 7, 1, 3, 7, 6, 3, 4, 1, 6)
FSTAR_UPPER = Fraction(664887507913, 1099511627776)


def parity(u, v):
    return (u & v).bit_count() & 1


def combine(x, y, ell):
    z = [0] * 8
    for r in range(8):
        for s in range(8):
            if parity(ell, r) == parity(ell, s):
                out = r ^ s ^ (E[ell] if parity(ell, s) else 0)
                z[out] += x[r] * y[s]
    return z


def run(F, fresh):
    F = Fraction(F)
    # Integer weights for F at the target label and uniform other labels.
    raw = [7 * F.numerator] + [F.denominator - F.numerator] * 7
    current = raw[:]
    cost = Fraction(1)
    for j, ell in enumerate(WORD, 1):
        mode_is_R = (j == 4 and fresh)
        right = raw if mode_is_R else current
        a = Fraction(sum(combine(current, right, ell)),
                     sum(current) * sum(right))
        if a <= 0:
            raise ArithmeticError(f"Zero acceptance at stage {j}")
        cost = (cost + 1 if mode_is_R else 2 * cost) / a
        current = combine(current, right, ell)
    require(sum(current) > 0, "zero accepted weight")
    return Fraction(current[0], sum(current)), cost


def dec(q, digits=15):
    with localcontext() as ctx:
        ctx.prec = digits + 25
        return f"{Decimal(q.numerator) / Decimal(q.denominator):.{digits}f}"


def main():
    for F in (Fraction(5995, 10000), Fraction(5997, 10000), Fraction(3, 5)):
        fh, ch = run(F, fresh=True)
        fs, cs = run(F, fresh=False)
        print(f"raw F={dec(F,4)}")
        print(f"  H12:       fidelity={dec(fh)}, expected raw cost={dec(ch,6)}")
        print(f"  all-S:     fidelity={dec(fs)}, expected raw cost={dec(cs,6)}")
        if F == Fraction(5995, 10000):
            require(fh > FSTAR_UPPER > fs, "H12 basin/control ordering failed")
            require(Fraction(551387, 1000000) < fs < Fraction(551388, 1000000), "F=0.5995 control bound failed")
            require(Fraction(605026, 1000000) < fh < Fraction(605027, 1000000), "F=0.5995 H12 bound failed")
        if F == Fraction(5997, 10000):
            require(fh > FSTAR_UPPER > fs, "H12 basin/control ordering failed")
            require(Fraction(562269, 1000000) < fs < Fraction(562270, 1000000), "F=0.5997 control bound failed")
        # Verify the exact small-cost ordering, not just the displayed rounding.
        require(ch < cs, "fresh-input stage did not lower cost for this word")
        require(0 < fh <= 1 and 0 < fs <= 1, "invalid output fidelity")
        print()
    print("ALL CONTROL COMPARISONS PASSED")


if __name__ == "__main__":
    main()
