#!/usr/bin/env python3
"""
Exact verifier for the heterogeneous-input Majorana-ancilla bootstrap
introduced in the September 2026 revision of the Phase-1 notes.

The theorem-level checks use only Python's standard library and exact
fractions.Fraction arithmetic.  The optional operation-noise continuation
at the end uses mpmath and is explicitly numerical, not part of the exact
certificate.
"""

from __future__ import annotations
from fractions import Fraction
from typing import Sequence

try:
    import mpmath as mp
except ImportError:
    mp = None

Vec = list[Fraction]
E_CANONICAL = {1: 1, 2: 2, 3: 1, 4: 4, 5: 1, 6: 2, 7: 1}

# S = self-pair the current stream; R = pair current stream with a fresh raw input.
H12: tuple[tuple[str, int], ...] = (
    ("S", 1), ("S", 2), ("S", 4), ("R", 7),
    ("S", 1), ("S", 3), ("S", 7), ("S", 6),
    ("S", 3), ("S", 4), ("S", 1), ("S", 6),
)

FSTAR_UPPER = Fraction(664887507913, 1099511627776)


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def dot(a: int, b: int) -> int:
    return (a & b).bit_count() & 1


def isotropic_fourier(F: Fraction) -> Vec:
    require(Fraction(0) <= F <= Fraction(1), "F outside [0,1]")
    x = (8 * F - 1) / 7
    return [Fraction(1)] + [x] * 7


def isotropic_probability(F: Fraction) -> Vec:
    require(Fraction(0) <= F <= Fraction(1), "F outside [0,1]")
    q = (1 - F) / 7
    return [F] + [q] * 7


def fourier_of_distribution(P: Sequence[Fraction]) -> Vec:
    require(len(P) == 8, "distribution must have length 8")
    require(sum(P) == 1, "distribution not normalized")
    require(all(p >= 0 for p in P), "negative probability")
    return [
        sum(((-1) ** dot(v, s)) * P[s] for s in range(8))
        for v in range(8)
    ]


def check_fourier(
    x1: Sequence[Fraction],
    x2: Sequence[Fraction],
    ell: int,
) -> tuple[Vec, Fraction]:
    e = E_CANONICAL[ell]
    den = 1 + x1[ell] * x2[ell]
    require(den > 0, "nonpositive conditional denominator")
    y = [Fraction(0)] * 8
    y[0] = Fraction(1)
    for v in range(1, 8):
        if dot(v, e) == 0:
            y[v] = (
                x1[v] * x2[v] + x1[v ^ ell] * x2[v ^ ell]
            ) / den
        else:
            y[v] = (
                x1[v] * x2[v ^ ell] + x1[v ^ ell] * x2[v]
            ) / den
    return y, den / 2


def check_probability(
    P1: Sequence[Fraction],
    P2: Sequence[Fraction],
    ell: int,
) -> tuple[Vec, Fraction]:
    e = E_CANONICAL[ell]
    out = [Fraction(0)] * 8
    acc = Fraction(0)
    for r in range(8):
        for s in range(8):
            w = P1[r] * P2[s]
            if dot(ell, r) == dot(ell, s):
                acc += w
                t = r ^ s ^ (e if dot(ell, s) else 0)
                out[t] += w
    require(acc > 0, "zero acceptance probability")
    return [z / acc for z in out], acc


def run_fourier(
    F: Fraction,
    schedule: Sequence[tuple[str, int]] = H12,
) -> tuple[Fraction, Fraction, int, list[Fraction], Vec]:
    raw = isotropic_fourier(F)
    x = list(raw)
    cost = Fraction(1)
    leaves = 1
    acceptances: list[Fraction] = []
    for mode, ell in schedule:
        if mode == "S":
            xnew, A = check_fourier(x, x, ell)
            cost = 2 * cost / A
            leaves *= 2
        elif mode == "R":
            xnew, A = check_fourier(x, raw, ell)
            cost = (cost + 1) / A
            leaves += 1
        else:
            raise VerificationError(f"unknown schedule mode {mode!r}")
        x = xnew
        acceptances.append(A)
    Fout = sum(x) / 8
    return Fout, cost, leaves, acceptances, x


def run_probability(
    F: Fraction,
    schedule: Sequence[tuple[str, int]] = H12,
) -> tuple[Fraction, Fraction, int, list[Fraction], Vec]:
    raw = isotropic_probability(F)
    P = list(raw)
    cost = Fraction(1)
    leaves = 1
    acceptances: list[Fraction] = []
    for mode, ell in schedule:
        if mode == "S":
            Pnew, A = check_probability(P, P, ell)
            cost = 2 * cost / A
            leaves *= 2
        elif mode == "R":
            Pnew, A = check_probability(P, raw, ell)
            cost = (cost + 1) / A
            leaves += 1
        else:
            raise VerificationError(f"unknown schedule mode {mode!r}")
        P = Pnew
        acceptances.append(A)
    return P[0], cost, leaves, acceptances, P


def decimal_enclosure(x: Fraction, digits: int) -> tuple[str, str]:
    scale = 10 ** digits
    q = (x.numerator * scale) // x.denominator
    lo = Fraction(q, scale)
    hi = lo if lo == x else Fraction(q + 1, scale)

    def render(v: Fraction) -> str:
        n = v.numerator * scale // v.denominator
        s = str(n).zfill(digits + 1)
        return s[:-digits] + "." + s[-digits:] if digits else s

    return render(lo), render(hi)


def exact_checks() -> None:
    F = Fraction(5995, 10000)
    FF, CF, NF, AF, x = run_fourier(F)
    FP, CP, NP, AP, P = run_probability(F)

    require(FF == FP, "Fourier and probability-space fidelities differ")
    require(CF == CP, "Fourier and probability-space cost recurrences differ")
    require(NF == NP == 2304, f"unexpected raw-leaf count: {NF}, {NP}")
    require(AF == AP, "acceptance probabilities differ between implementations")
    require(FF > FSTAR_UPPER, "heterogeneous endpoint does not enter R4 basin")

    # Check the Fourier/probability representations of the final state exactly.
    require(list(x) == fourier_of_distribution(P),
            "final Fourier vector does not match direct probability calculation")

    flo, fhi = decimal_enclosure(FF, 30)
    ylo, yhi = decimal_enclosure(1 / CF, 20)
    clo, chi = decimal_enclosure(CF, 6)

    print("H12 exact theorem at F = 5995/10000")
    print("schedule =", H12)
    print("raw leaves in one fully accepted schedule =", NF)
    print("Fout in [", flo, ",", fhi, "]")
    print("Fstar upper =", FSTAR_UPPER)
    print("strict basin entry =", FF > FSTAR_UPPER)
    print("expected raw-copy cost in [", clo, ",", chi, "]")
    print("pooled yield in [", ylo, ",", yhi, "]")
    print("reduced Fout numerator bits =", FF.numerator.bit_length())
    print("reduced Fout denominator bits =", FF.denominator.bit_length())

    # Same schedule at the old W25 endpoint, useful for the robustness comparison.
    Fold = Fraction(5997, 10000)
    Fo, Co, No, _, _ = run_fourier(Fold)
    folo, fohi = decimal_enclosure(Fo, 30)
    yolo, yohi = decimal_enclosure(1 / Co, 20)
    print()
    print("H12 at F = 5997/10000")
    print("Fout in [", folo, ",", fohi, "]")
    print("pooled yield in [", yolo, ",", yohi, "]")
    require(Fo > FSTAR_UPPER, "H12 fails at 5997/10000")


# ---------------------- numerical operation-noise section ----------------------

def _iso_mp(F):
    x = (8 * F - 1) / 7
    return [mp.mpf(1)] + [x] * 7


def _check_mp(x1, x2, ell, lam):
    e = E_CANONICAL[ell]
    den = 1 + x1[ell] * x2[ell]
    kappa = 1 - mp.mpf(8) * lam / 7
    y = [mp.mpf(1)] + [mp.mpf(0)] * 7
    for v in range(1, 8):
        if dot(v, e) == 0:
            ideal = (
                x1[v] * x2[v] + x1[v ^ ell] * x2[v ^ ell]
            ) / den
        else:
            ideal = (
                x1[v] * x2[v ^ ell] + x1[v ^ ell] * x2[v]
            ) / den
        y[v] = kappa * ideal
    return y, den / 2


def _protocol_mp(F, lam, schedule=H12):
    raw = _iso_mp(F)
    x = list(raw)
    for mode, ell in schedule:
        if mode == "S":
            x, _ = _check_mp(x, x, ell, lam)
        else:
            x, _ = _check_mp(x, raw, ell, lam)
    return mp.fsum(x) / 8


def _r4_noisy_mp(F, lam):
    x = _iso_mp(F)
    for ell in (1, 2, 4, 7):
        x, _ = _check_mp(x, x, ell, lam)
    return mp.fsum(x) / 8


def _lower_root(lam):
    f = lambda z: _r4_noisy_mp(z, lam) - z
    guess = mp.mpf("0.6047116657215318756") + mp.mpf("2") * lam
    return mp.findroot(
        f,
        (guess - mp.mpf("0.01"), guess + mp.mpf("0.01")),
    )


def numerical_checks() -> None:
    if mp is None:
        print("\nmpmath unavailable: numerical noise checks skipped")
        return
    mp.mp.dps = 70

    Fstar = _lower_root(mp.mpf("0"))
    hF = lambda F: _protocol_mp(F, mp.mpf("0")) - Fstar
    Fcross = mp.findroot(
        hF, (mp.mpf("0.59945"), mp.mpf("0.59955"))
    )

    def noise_cross(Fstr: str, a: str, b: str):
        F = mp.mpf(Fstr)
        h = lambda lam: _protocol_mp(F, lam) - _lower_root(lam)
        return mp.findroot(h, (mp.mpf(a), mp.mpf(b)))

    l5997 = noise_cross("0.5997", "1e-4", "2e-4")
    l5995 = noise_cross("0.5995", "1e-6", "1e-5")

    print()
    print("NUMERICAL ONLY (not part of the exact theorem)")
    print("H12 noiseless basin crossing F =", mp.nstr(Fcross, 40))
    print("H12 operation-noise crossing at F=0.5997 =", mp.nstr(l5997, 40))
    print("H12 operation-noise crossing at F=0.5995 =", mp.nstr(l5995, 40))


if __name__ == "__main__":
    exact_checks()
    numerical_checks()
    print("\nALL HETEROGENEOUS-BOOTSTRAP CHECKS PASSED")
