#!/usr/bin/env python3
"""
Extended independent verifier for the research notes
"Purifying a noisy non-Gaussian ancilla for Majorana quantum computation".

This revision preserves the mathematical content of the supplied verifier and
strengthens its traceability, test coverage, and failure behaviour.

Requirements
------------
    Python >= 3.10
    SymPy  >= 1.12
    mpmath (tested with 1.3.0; no minimum version is asserted here)

Proof-status convention
-----------------------
* Exact/certificate checks use integers, fractions.Fraction, or exact SymPy
  expressions.
* The W25 certificate uses rational interval arithmetic with explicit outward
  rounding to the decimal grid 10^(-80). No binary floating-point number is
  used in that certificate.
* Results for the effective post-check syndrome-noise channel are numerical
  high-precision checks performed with mpmath. They are not promoted here to
  interval-certified theorems.

Important implementation note
-----------------------------
The verifier deliberately does NOT use Python ``assert`` statements for
scientific checks. Assertions can be removed by ``python -O``. All verification
conditions are enforced by explicit calls to ``require`` and therefore remain
active in optimized mode.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import platform
import re
import sys
from typing import Callable, Sequence, cast

import mpmath as mp
import sympy as sp


__version__ = "1.2"


class VerificationError(RuntimeError):
    """Raised when a scientific verification condition is not satisfied."""


def require(condition: bool, message: str) -> None:
    """Enforce a verification condition even when Python is run with ``-O``."""
    if not condition:
        raise VerificationError(message)


def _numeric_version(version: str) -> tuple[int, ...]:
    """Extract the leading numeric components of a package version string."""
    parts = re.findall(r"\d+", version)
    return tuple(int(p) for p in parts[:3])


def verify_environment() -> None:
    """Check requirements that are actually documented as minima.

    The supplied verifier documents Python >= 3.10 and SymPy >= 1.12.
    mpmath is an explicit dependency of the full numerical run, but this
    extended verifier does not claim an independently established minimum
    mpmath version; the reference run records the tested version instead.
    """
    require(sys.version_info >= (3, 10), "Python >= 3.10 is required")
    require(_numeric_version(sp.__version__) >= (1, 12), "SymPy >= 1.12 is required")


def environment_string() -> str:
    """Return a compact, reproducibility-oriented environment description."""
    return (
        f"Python {platform.python_version()}; "
        f"SymPy {sp.__version__}; mpmath {mp.__version__}; "
        f"{platform.platform()}"
    )


# ---------------------------------------------------------------------------
# F_2^3 utilities and the Majorana label table
# ---------------------------------------------------------------------------

Vec3 = tuple[int, int, int]
Matrix = tuple[Vec3, Vec3, Vec3]
I3: Matrix = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


def dot(a: int, b: int) -> int:
    """Binary scalar product on F_2^3 for integers encoded as three bits."""
    return (a & b).bit_count() & 1


def int_to_vec(x: int) -> Vec3:
    require(0 <= x < 8, f"F_2^3 label outside 0..7: {x}")
    return ((x >> 2) & 1, (x >> 1) & 1, x & 1)


def vec_to_int(v: Vec3) -> int:
    require(len(v) == 3 and all(bit in (0, 1) for bit in v), f"not a binary 3-vector: {v}")
    return (v[0] << 2) | (v[1] << 1) | v[2]


def mat_vec(M: Matrix, x: int) -> int:
    v = int_to_vec(x)
    out = tuple(
        sum(M[i][j] * v[j] for j in range(3)) & 1
        for i in range(3)
    )
    return vec_to_int(cast(Vec3, out))


def mat_mul(A: Matrix, B: Matrix) -> Matrix:
    out = tuple(
        tuple(
            sum(A[i][k] * B[k][j] for k in range(3)) & 1
            for j in range(3)
        )
        for i in range(3)
    )
    return cast(Matrix, out)


def transpose(M: Matrix) -> Matrix:
    return cast(Matrix, tuple(tuple(M[j][i] for j in range(3)) for i in range(3)))


def rank_mod2(M: Matrix) -> int:
    """Gaussian elimination over F_2."""
    rows = [list(row) for row in M]
    rank = 0
    for col in range(3):
        pivot = next((r for r in range(rank, 3) if rows[r][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for r in range(3):
            if r != rank and rows[r][col]:
                rows[r] = [a ^ b for a, b in zip(rows[r], rows[rank])]
        rank += 1
    return rank


def all_gl3() -> list[Matrix]:
    out: list[Matrix] = []
    for bits in product((0, 1), repeat=9):
        M = cast(Matrix, (bits[0:3], bits[3:6], bits[6:9]))
        if rank_mod2(M) == 3:
            out.append(M)
    return out


GL3 = all_gl3()
require(len(GL3) == 168, f"expected |GL(3,2)|=168, obtained {len(GL3)}")


@lru_cache(maxsize=None)
def inverse(M: Matrix) -> Matrix:
    """Return the inverse of an element of GL(3,2)."""
    for N in GL3:
        if mat_mul(M, N) == I3 and mat_mul(N, M) == I3:
            return N
    raise VerificationError("matrix declared invertible has no inverse")


# Supports from the notes, Eqs. (4)-(6):
# S1=-c1 c2 c5 c6, S2=-c2 c3 c6 c7, S3=-c1 c2 c3 c4.
STABILIZER_SUPPORTS = (
    frozenset((1, 2, 5, 6)),
    frozenset((2, 3, 6, 7)),
    frozenset((1, 2, 3, 4)),
)

# Labels in Majorana-index order p=1,...,8, as in Eq. (16) of the notes.
LABELS = [0b101, 0b111, 0b011, 0b001, 0b100, 0b110, 0b010, 0b000]


def majorana_commutes_with_monomial(p: int, support: frozenset[int]) -> bool:
    """Return whether c_p commutes with c_A, using the Majorana sign rule."""
    require(1 <= p <= 8, "Majorana index must be in 1..8")
    # c_p c_A = (-1)^(|A|-1_{p in A}) c_A c_p for |{p}|=1.
    exponent = len(support) - (1 if p in support else 0)
    return exponent % 2 == 0


def derive_labels_from_stabilizers() -> list[int]:
    """Derive lambda(p) from commutation/anticommutation with S1,S2,S3."""
    labels: list[int] = []
    for p in range(1, 9):
        bits = tuple(
            0 if majorana_commutes_with_monomial(p, support) else 1
            for support in STABILIZER_SUPPORTS
        )
        labels.append(vec_to_int(cast(Vec3, bits)))
    return labels


def permutation_parity(perm: Sequence[int]) -> int:
    inversions = sum(
        perm[i] > perm[j]
        for i in range(len(perm))
        for j in range(i + 1, len(perm))
    )
    return inversions & 1


def verify_basic_f2_and_labels() -> None:
    """Verify encoding identities, dot products, and the Majorana label table."""
    for x in range(8):
        require(vec_to_int(int_to_vec(x)) == x, f"int/vector round trip failed for {x}")
    for a in range(8):
        va = int_to_vec(a)
        for b in range(8):
            vb = int_to_vec(b)
            reference = sum(x * y for x, y in zip(va, vb)) % 2
            require(dot(a, b) == reference, f"dot product mismatch for {a},{b}")
    derived = derive_labels_from_stabilizers()
    require(derived == LABELS, f"Majorana label table mismatch: derived {derived}, stored {LABELS}")
    require(len(set(LABELS)) == 8, "Majorana labels are not a bijection onto F_2^3")


def verify_affine_structure() -> None:
    """Exhaust the finite affine geometry and the 28 conjugate checks."""
    label_to_index = {x: i for i, x in enumerate(LABELS)}
    require(len(label_to_index) == 8, "label_to_index is not bijective")

    for M in GL3:
        perm = [label_to_index[mat_vec(M, x)] for x in LABELS]
        require(permutation_parity(perm) == 0, f"odd GL(3,2) label permutation for M={M}")

    planes = {
        frozenset(x for x in range(8) if dot(v, x) == b)
        for v in range(1, 8)
        for b in (0, 1)
    }
    require(len(planes) == 14, f"expected 14 affine planes, obtained {len(planes)}")

    affine_maps = set()
    for M in GL3:
        for t in range(8):
            image_tuple = tuple(mat_vec(M, x) ^ t for x in range(8))
            affine_maps.add(image_tuple)
            for P in planes:
                image = frozenset(image_tuple[x] for x in P)
                require(image in planes, f"affine map failed to preserve planes: M={M}, t={t}, P={P}")
    require(len(affine_maps) == 1344, f"expected |AGL(3,2)|=1344, obtained {len(affine_maps)}")

    # Converse normalizer: exhaustive over all 8! point permutations.
    normalizer = set()
    for perm in permutations(range(8)):
        if all(frozenset(perm[x] for x in P) in planes for P in planes):
            normalizer.add(perm)
    require(normalizer == affine_maps, "plane-preserving permutation normalizer is not AGL(3,2)")

    # Every admissible (ell,e) is a conjugate of the base pair a=a=(001).
    a = 0b001
    realized = set()
    for M in GL3:
        ell = mat_vec(transpose(M), a)
        e = mat_vec(inverse(M), a)
        require(dot(ell, e) == 1, f"inadmissible conjugated pair ell={ell}, e={e}")
        realized.add((ell, e))
    admissible = {
        (ell, e)
        for ell in range(1, 8)
        for e in range(8)
        if dot(ell, e) == 1
    }
    require(realized == admissible, "conjugation does not realize exactly all admissible checks")
    require(len(realized) == 28, f"expected 28 checks, obtained {len(realized)}")

    # Check the conjugated syndrome map on all M and all 64 input pairs.
    for M in GL3:
        Minv = inverse(M)
        ell = mat_vec(transpose(M), a)
        e = mat_vec(Minv, a)
        for r in range(8):
            for s in range(8):
                r1, s1 = mat_vec(M, r), mat_vec(M, s)
                acc_base = dot(a, r1) == dot(a, s1)
                acc_general = dot(ell, r) == dot(ell, s)
                require(acc_base == acc_general, "base/general acceptance mismatch")
                if acc_base:
                    u1 = r1 ^ s1 ^ (a if dot(a, s1) else 0)
                    u_conj = mat_vec(Minv, u1)
                    u_general = r ^ s ^ (e if dot(ell, s) else 0)
                    require(u_conj == u_general, "base/general output-syndrome mismatch")


# ---------------------------------------------------------------------------
# Fourier check dynamics
# ---------------------------------------------------------------------------

E_CANONICAL = {1: 1, 2: 2, 3: 1, 4: 4, 5: 1, 6: 2, 7: 1}


def _validate_check(ell: int, e: int) -> None:
    require(1 <= ell <= 7, f"ell must be nonzero in F_2^3, obtained {ell}")
    require(0 <= e <= 7, f"e must lie in F_2^3, obtained {e}")
    require(dot(ell, e) == 1, f"check requires <ell,e>=1, obtained ell={ell}, e={e}")


def check_fraction(
    x1: list[Fraction],
    x2: list[Fraction],
    ell: int,
    e: int | None = None,
) -> tuple[list[Fraction], Fraction]:
    if e is None:
        e = E_CANONICAL[ell]
    _validate_check(ell, e)
    require(len(x1) == 8 and len(x2) == 8, "Fourier vectors must have length 8")
    den = Fraction(1) + x1[ell] * x2[ell]
    if den <= 0:
        raise VerificationError("nonpositive conditional denominator")
    y = [Fraction(0)] * 8
    y[0] = Fraction(1)
    for v in range(1, 8):
        if dot(v, e) == 0:
            y[v] = (x1[v] * x2[v] + x1[v ^ ell] * x2[v ^ ell]) / den
        else:
            y[v] = (x1[v] * x2[v ^ ell] + x1[v ^ ell] * x2[v]) / den
    return y, den / 2


def fourier_of_distribution(P: list[Fraction]) -> list[Fraction]:
    require(len(P) == 8, f"distribution length must be 8, obtained {len(P)}")
    require(sum(P) == 1, f"distribution is not normalized: sum={sum(P)}")
    require(all(p >= 0 for p in P), "distribution contains a negative probability")
    return [
        sum(((-1) ** dot(v, s)) * P[s] for s in range(8))
        for v in range(8)
    ]


def direct_check_distribution(
    P1: list[Fraction],
    P2: list[Fraction],
    ell: int,
    e: int,
) -> tuple[list[Fraction], Fraction]:
    _validate_check(ell, e)
    require(len(P1) == 8 and len(P2) == 8, "input distributions must have length 8")
    require(sum(P1) == 1 and all(p >= 0 for p in P1), "P1 is not a probability distribution")
    require(sum(P2) == 1 and all(p >= 0 for p in P2), "P2 is not a probability distribution")
    out = [Fraction(0)] * 8
    acc = Fraction(0)
    for r in range(8):
        for s in range(8):
            w = P1[r] * P2[s]
            if dot(ell, r) == dot(ell, s):
                acc += w
                t = r ^ s ^ (e if dot(ell, s) else 0)
                out[t] += w
    if acc <= 0:
        raise VerificationError("conditional check has zero acceptance probability")
    return [z / acc for z in out], acc


def verify_fourier_formula() -> None:
    """Verify Theorem 6.7 both by regression and exhaustively on basis states.

    The unnormalized acceptance/Fourier numerators are bilinear in the two input
    distributions. Therefore checking all 8x8 pairs of point distributions for
    every admissible (ell,e) verifies those bilinear identities on a basis; by
    linearity they then hold for arbitrary input distributions. Conditional
    formulas follow whenever acceptance is nonzero.
    """
    # A generic exact-rational regression pair.
    P1raw = [Fraction(j + 1, 36) for j in range(8)]
    P2raw = [Fraction(8 - j, 36) for j in range(8)]
    x1 = fourier_of_distribution(P1raw)
    x2 = fourier_of_distribution(P2raw)
    for ell in range(1, 8):
        for e in range(8):
            if dot(ell, e) != 1:
                continue
            Pout, acc_direct = direct_check_distribution(P1raw, P2raw, ell, e)
            x_direct = fourier_of_distribution(Pout)
            x_formula, acc_formula = check_fraction(x1, x2, ell, e)
            require(acc_direct == acc_formula, f"acceptance formula mismatch for ell={ell}, e={e}")
            require(x_direct == x_formula, f"Fourier update mismatch for ell={ell}, e={e}")

    # Exhaustive basis verification of the underlying bilinear identities.
    for ell in range(1, 8):
        for e in range(8):
            if dot(ell, e) != 1:
                continue
            for r in range(8):
                x_r = [Fraction((-1) ** dot(v, r)) for v in range(8)]
                for s in range(8):
                    x_s = [Fraction((-1) ** dot(v, s)) for v in range(8)]
                    accepted = dot(ell, r) == dot(ell, s)
                    acc_direct = Fraction(1 if accepted else 0)
                    acc_formula = Fraction(1 + x_r[ell] * x_s[ell], 2)
                    require(acc_direct == acc_formula, "basis acceptance identity failed")
                    t = r ^ s ^ (e if dot(ell, s) else 0) if accepted else 0
                    for v in range(8):
                        if dot(v, e) == 0:
                            numerator = Fraction(
                                x_r[v] * x_s[v] + x_r[v ^ ell] * x_s[v ^ ell], 2
                            )
                        else:
                            numerator = Fraction(
                                x_r[v] * x_s[v ^ ell] + x_r[v ^ ell] * x_s[v], 2
                            )
                        direct_numerator = (
                            Fraction((-1) ** dot(v, t)) if accepted else Fraction(0)
                        )
                        require(
                            numerator == direct_numerator,
                            f"basis Fourier numerator failed for ell={ell}, e={e}, r={r}, s={s}, v={v}",
                        )


def isotropic_fraction(F: Fraction) -> list[Fraction]:
    require(Fraction(0) <= F <= Fraction(1), f"fidelity outside [0,1]: {F}")
    x = (8 * F - 1) / 7
    return [Fraction(1)] + [x] * 7


def _validate_check_word(word: Sequence[int], *, allow_empty: bool = True) -> None:
    """Reject malformed check sequences before a calculation starts."""
    require(not isinstance(word, (str, bytes)), "check word must be a sequence of integer labels")
    if not allow_empty:
        require(len(word) > 0, "check word must not be empty")
    require(
        all(type(ell) is int and ell in E_CANONICAL for ell in word),
        f"check word contains a noncanonical label: {list(word)}",
    )


def apply_word_fraction(
    F: Fraction,
    word: list[int],
) -> tuple[Fraction, list[Fraction], list[Fraction]]:
    _validate_check_word(word)
    x = isotropic_fraction(F)
    acceptances: list[Fraction] = []
    for ell in word:
        x, A = check_fraction(x, x, ell)
        acceptances.append(A)
    return sum(x) / 8, x, acceptances


# ---------------------------------------------------------------------------
# Exact low-weight fault enumeration
# ---------------------------------------------------------------------------


def check_syndrome_pair(r: int, s: int, ell: int) -> int | None:
    require(0 <= r < 8 and 0 <= s < 8, "syndrome labels must lie in 0..7")
    require(ell in E_CANONICAL, f"unknown canonical check ell={ell}")
    e = E_CANONICAL[ell]
    if dot(ell, r) != dot(ell, s):
        return None
    return r ^ s ^ (e if dot(ell, s) else 0)


def propagate_tree(leaves: list[int], word: list[int]) -> int | None:
    _validate_check_word(word)
    current = list(leaves)
    expected_leaves = 2 ** len(word)
    require(
        len(current) == expected_leaves,
        f"tree of depth {len(word)} requires {expected_leaves} leaves, obtained {len(current)}",
    )
    require(
        all(type(x) is int and 0 <= x < 8 for x in current),
        "leaf syndromes must be integer labels in 0..7",
    )
    for ell in word:
        if len(current) % 2:
            raise VerificationError("tree level has odd width")
        nxt: list[int] = []
        for j in range(0, len(current), 2):
            out = check_syndrome_pair(current[j], current[j + 1], ell)
            if out is None:
                return None
            nxt.append(out)
        current = nxt
    if len(current) != 1:
        raise VerificationError("word depth does not match leaf count")
    return current[0]


def low_weight_counts(word: list[int], max_weight: int) -> dict[int, tuple[int, int]]:
    _validate_check_word(word)
    require(type(max_weight) is int, "max_weight must be an integer")
    n = 2 ** len(word)
    require(0 <= max_weight <= n, f"max_weight must lie in 0..{n}")
    result: dict[int, tuple[int, int]] = {}
    for weight in range(max_weight + 1):
        correct = 0
        wrong = 0
        visited = 0
        for positions in combinations(range(n), weight):
            for values in product(range(1, 8), repeat=weight):
                visited += 1
                leaves = [0] * n
                for position, value in zip(positions, values):
                    leaves[position] = value
                out = propagate_tree(leaves, word)
                if out is None:
                    continue
                if out == 0:
                    correct += 1
                else:
                    wrong += 1
        expected_visited = sp.binomial(n, weight) * (7 ** weight)
        require(visited == int(expected_visited), "fault enumeration did not visit the complete labeled set")
        result[weight] = (correct, wrong)
    return result


def verify_fault_table() -> None:
    r3 = low_weight_counts([1, 2, 4], 2)
    r4 = low_weight_counts([1, 2, 4, 7], 3)
    require(
        r3 == {0: (1, 0), 1: (0, 0), 2: (20, 48)},
        f"R3 low-weight table mismatch: {r3}",
    )
    require(
        r4 == {0: (1, 0), 1: (0, 0), 2: (40, 0), 3: (32, 384)},
        f"R4 low-weight table mismatch: {r4}",
    )


# ---------------------------------------------------------------------------
# Readout map, correlated shifts, and effective syndrome-shift anatomy
# ---------------------------------------------------------------------------


def verify_readout_error_map() -> None:
    """Verify consequences of Eq. (23) in the notes exactly.

    This function does not reconstruct Eq. (23) from a microscopic braid/fusion
    circuit. It takes the stated classical feed-forward error map
        delta_ff=(m1+m2, m2+m3, 0)
    as input and exhaustively verifies its Hamming-weight consequences and the
    exact small-mu series of the conditional readout error.
    """
    weight_two_deltas: list[tuple[int, int, int]] = []
    for mask in range(16):
        m1 = (mask >> 3) & 1
        m2 = (mask >> 2) & 1
        m3 = (mask >> 1) & 1
        m4 = mask & 1
        weight = m1 + m2 + m3 + m4
        parity = weight & 1
        delta = (m1 ^ m2, m2 ^ m3, 0)
        if weight == 1:
            require(parity == 1, "single readout error should flip checksum parity")
        if weight == 2:
            require(parity == 0, "two readout errors should preserve checksum parity")
            require(delta != (0, 0, 0), "two readout errors unexpectedly cancel feed-forward syndrome")
            weight_two_deltas.append(delta)
        if weight == 4:
            require(parity == 0, "four readout errors should preserve checksum parity")
            require(delta == (0, 0, 0), "four readout errors should cancel feed-forward syndrome")

    expected_multiset = sorted(
        [(0, 1, 0), (1, 1, 0), (1, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
    )
    require(sorted(weight_two_deltas) == expected_multiset, "weight-two readout syndromes mismatch Eq. (25) proof")

    mu = sp.symbols("mu")
    A = (1 - mu) ** 4 + 6 * mu ** 2 * (1 - mu) ** 2 + mu ** 4
    err = sp.cancel(6 * mu ** 2 * (1 - mu) ** 2 / A)
    series = sp.series(err, mu, 0, 4).removeO()
    require(sp.expand(series).coeff(mu, 2) == 6, "readout error has wrong quadratic coefficient")
    require(sp.series(err, mu, 0, 3).removeO() == 6 * mu ** 2, "readout O(mu^2) expansion mismatch")


def verify_common_mode_identity() -> None:
    for ell in range(1, 8):
        for e in range(8):
            if dot(ell, e) != 1:
                continue
            for r in range(8):
                for s in range(8):
                    for d in range(8):
                        acc0 = dot(ell, r) == dot(ell, s)
                        acc1 = dot(ell, r ^ d) == dot(ell, s ^ d)
                        require(acc0 == acc1, "common-mode shift changed acceptance")
                        if acc0:
                            out0 = r ^ s ^ (e if dot(ell, s) else 0)
                            out1 = (r ^ d) ^ (s ^ d) ^ (
                                e if dot(ell, s ^ d) else 0
                            )
                            expected = out0 ^ (e if dot(ell, d) else 0)
                            require(out1 == expected, "common-mode output identity failed")


def propagate_tree_with_node_shift(
    word: list[int],
    shift_level: int,
    shift_index: int,
    shift_syndrome: int,
) -> int | None:
    _validate_check_word(word, allow_empty=False)
    require(type(shift_level) is int, "shift level must be an integer")
    require(type(shift_index) is int, "shift index must be an integer")
    require(type(shift_syndrome) is int, "shift syndrome must be an integer")
    require(0 <= shift_level < len(word), "shift level outside tree")
    require(0 <= shift_index < 2 ** (len(word) - shift_level - 1), "shift index outside level")
    require(1 <= shift_syndrome <= 7, "shift syndrome must be nonzero")
    current = [0] * (2 ** len(word))
    for level, ell in enumerate(word):
        nxt: list[int] = []
        for j in range(0, len(current), 2):
            out = check_syndrome_pair(current[j], current[j + 1], ell)
            if out is None:
                return None
            index = j // 2
            if level == shift_level and index == shift_index:
                out ^= shift_syndrome
            nxt.append(out)
        current = nxt
    require(len(current) == 1, "syndrome-shift tree did not end in one output")
    return current[0]

def effective_syndrome_shift_rows(word: list[int]) -> list[tuple[int, int, int, int]]:
    _validate_check_word(word, allow_empty=False)
    rows: list[tuple[int, int, int, int]] = []
    depth = len(word)
    for level in range(depth):
        locations = 2 ** (depth - level - 1)
        wrong = 0
        rejected = 0
        correct = 0
        for index in range(locations):
            for syndrome in range(1, 8):
                out = propagate_tree_with_node_shift(word, level, index, syndrome)
                if out is None:
                    rejected += 1
                elif out == 0:
                    correct += 1
                else:
                    wrong += 1
        require(correct == 0, f"single effective syndrome shift corrected unexpectedly at level {level+1}")
        require(wrong + rejected == locations * 7, "syndrome-shift accounting incomplete")
        rows.append((level + 1, locations, wrong, rejected))
    return rows


def verify_effective_syndrome_shift_table() -> None:
    rows = effective_syndrome_shift_rows([1, 2, 4, 7])
    expected = [
        (1, 8, 0, 56),
        (2, 4, 4, 24),
        (3, 2, 6, 8),
        (4, 1, 7, 0),
    ]
    require(rows == expected, f"effective syndrome-shift table mismatch: {rows}")
    require(sum(row[2] for row in rows) == 17, "expected 17 accepted-wrong effective syndrome shifts")
    require(sum(row[3] for row in rows) == 88, "expected 88 rejected effective syndrome shifts")


def _expect_verification_error(label: str, action: Callable[[], object]) -> None:
    """Require a deliberately malformed input to fail with VerificationError."""
    try:
        action()
    except VerificationError:
        return
    except Exception as exc:
        raise VerificationError(
            f"invalid-input guard {label!r} raised {type(exc).__name__} instead of VerificationError"
        ) from exc
    raise VerificationError(f"invalid-input guard {label!r} accepted malformed input")


def verify_invalid_input_guards() -> int:
    """Exercise fail-fast guards with representative malformed tree inputs."""
    word = [1, 2, 4, 7]
    cases: list[tuple[str, Callable[[], object]]] = [
        ("negative shift level", lambda: propagate_tree_with_node_shift(word, -1, 0, 1)),
        ("shift level past root", lambda: propagate_tree_with_node_shift(word, len(word), 0, 1)),
        ("noninteger shift level", lambda: propagate_tree_with_node_shift(word, 0.0, 0, 1)),
        ("noninteger shift index", lambda: propagate_tree_with_node_shift(word, 0, 0.0, 1)),
        ("boolean shift syndrome", lambda: propagate_tree_with_node_shift(word, 0, 0, True)),
        ("zero shift syndrome", lambda: propagate_tree_with_node_shift(word, 0, 0, 0)),
        ("shift syndrome above F_2^3", lambda: propagate_tree_with_node_shift(word, 0, 0, 8)),
        ("empty shifted tree word", lambda: propagate_tree_with_node_shift([], 0, 0, 1)),
        ("noncanonical zero check", lambda: propagate_tree_with_node_shift([1, 2, 0, 7], 0, 0, 1)),
        ("noncanonical check above range", lambda: propagate_tree_with_node_shift([1, 2, 8, 7], 0, 0, 1)),
        ("boolean check label", lambda: propagate_tree_with_node_shift([1, 2, 4, True], 0, 0, 1)),
        ("wrong leaf count", lambda: propagate_tree([0] * 15, word)),
        ("negative leaf syndrome", lambda: propagate_tree([-1] + [0] * 15, word)),
        ("leaf syndrome above F_2^3", lambda: propagate_tree([8] + [0] * 15, word)),
        ("boolean leaf syndrome", lambda: propagate_tree([True] + [0] * 15, word)),
        ("negative max_weight", lambda: low_weight_counts([1], -1)),
        ("max_weight above leaf count", lambda: low_weight_counts([1], 3)),
        ("noninteger max_weight", lambda: low_weight_counts([1], 1.0)),
    ]
    for level in range(len(word)):
        locations = 2 ** (len(word) - level - 1)
        cases.extend(
            [
                (
                    f"negative shift index at level {level}",
                    lambda level=level: propagate_tree_with_node_shift(word, level, -1, 1),
                ),
                (
                    f"first out-of-range shift index at level {level}",
                    lambda level=level, locations=locations: propagate_tree_with_node_shift(
                        word, level, locations, 1
                    ),
                ),
            ]
        )

    for label, action in cases:
        _expect_verification_error(label, action)
    return len(cases)


# ---------------------------------------------------------------------------
# Exact R3/R4 maps, root structure, and local series
# ---------------------------------------------------------------------------

G7_COEFFS = [
    1458176,
    -884736,
    548864,
    -368640,
    138496,
    -67808,
    5440,
    -6249,
]

G15_COEFFS = [
    8505108987904,
    -16854525411328,
    28467580108800,
    -27066279919616,
    18851022503936,
    -10509663338496,
    4884383924224,
    -2048766509056,
    716242124800,
    -238242988032,
    58392420352,
    -20373018624,
    4691993600,
    -2122800736,
    263381416,
    -149949201,
]


def check_sym(x: list[sp.Expr], ell: int) -> list[sp.Expr]:
    require(ell in E_CANONICAL, f"unknown canonical check ell={ell}")
    e = E_CANONICAL[ell]
    den = 1 + x[ell] ** 2
    y: list[sp.Expr] = [sp.Integer(1)] + [sp.Integer(0)] * 7
    for v in range(1, 8):
        if dot(v, e) == 0:
            y[v] = sp.cancel((x[v] ** 2 + x[v ^ ell] ** 2) / den)
        else:
            y[v] = sp.cancel(2 * x[v] * x[v ^ ell] / den)
    return y


def exact_round_map(word: list[int]) -> tuple[sp.Symbol, sp.Expr]:
    _validate_check_word(word)
    F = sp.symbols("F")
    x0 = (8 * F - 1) / 7
    x: list[sp.Expr] = [sp.Integer(1)] + [x0] * 7
    for ell in word:
        x = check_sym(x, ell)
    return F, sp.cancel(sum(x) / 8)


def _verify_factorized_round(
    word: list[int],
    coeffs: list[int],
    lower: Fraction,
    upper: Fraction,
    sample_F: Fraction,
    label: str,
) -> tuple[sp.Symbol, sp.Expr, sp.Expr]:
    F, R = exact_round_map(word)
    g = sp.Poly.from_list(coeffs, gens=F).as_expr()
    numerator, denominator = sp.fraction(sp.factor(R - F))
    expected = -(F - 1) * (8 * F - 1) * g
    require(sp.expand(numerator - expected) == 0, f"{label} numerator factorization mismatch")

    den_poly = sp.Poly(denominator, F)
    require(den_poly.count_roots(0, 1) == 0, f"{label} denominator has a root in [0,1]")
    require(den_poly.eval(Fraction(0)) > 0, f"{label} denominator not positive at F=0")
    require(den_poly.eval(Fraction(1, 2)) > 0, f"{label} denominator not positive at F=1/2")
    require(den_poly.eval(Fraction(1)) > 0, f"{label} denominator not positive at F=1")

    g_poly = sp.Poly(g, F)
    require(g_poly.count_roots(0, 1) == 1, f"{label} g polynomial should have exactly one root in (0,1)")
    require(Fraction(0) < lower < upper < Fraction(1), f"{label} invalid root bracket")
    require(sp.sign(g_poly.eval(lower)) != sp.sign(g_poly.eval(upper)), f"{label} root bracket has no sign change")
    require(g_poly.eval(0) < 0 and g_poly.eval(1) > 0, f"{label} endpoint signs unexpected")

    sample_expr = sp.cancel(R.subs(F, sp.Rational(sample_F.numerator, sample_F.denominator)) - sp.Rational(sample_F.numerator, sample_F.denominator))
    require(sample_expr > 0, f"{label} does not improve the exact sample point F={sample_F}")
    return F, R, g


def verify_round_maps() -> None:
    # R4 theorem and bracket from the notes.
    F4, R4, _ = _verify_factorized_round(
        [1, 2, 4, 7],
        G15_COEFFS,
        Fraction(1329775015825, 2199023255552),
        Fraction(664887507913, 1099511627776),
        Fraction(3, 4),
        "R4",
    )
    exact_r4_sample = sp.cancel(R4.subs(F4, sp.Rational(3, 4)))
    require(
        exact_r4_sample == sp.Rational(1466150375091, 1556750188228),
        "R4(3/4) does not match the exact value quoted in the notes",
    )

    # R3 / Bravyi recurrence threshold in the present parametrization.
    F3, R3, _ = _verify_factorized_round(
        [1, 2, 4],
        G7_COEFFS,
        Fraction("0.61560748080966942"),
        Fraction("0.61560748080966943"),
        Fraction(3, 4),
        "R3",
    )

    eps = sp.symbols("eps")
    series4 = sp.series(1 - R4.subs(F4, 1 - eps), eps, 0, 6).removeO()
    expected4 = (
        sp.Rational(384, 343) * eps ** 3
        + sp.Rational(11520, 2401) * eps ** 4
        + sp.Rational(218112, 16807) * eps ** 5
    )
    require(sp.expand(series4 - expected4) == 0, "R4 Taylor series mismatch")

    series3 = sp.series(1 - R3.subs(F3, 1 - eps), eps, 0, 4).removeO()
    expected3 = sp.Rational(48, 49) * eps ** 2 + sp.Rational(960, 343) * eps ** 3
    require(sp.expand(series3 - expected3) == 0, "R3 Taylor series mismatch")


# ---------------------------------------------------------------------------
# W14 exact certificate
# ---------------------------------------------------------------------------


def minimal_unsigned_big_endian(n: int) -> bytes:
    if n <= 0:
        raise ValueError("expected a positive integer")
    return n.to_bytes((n.bit_length() + 7) // 8, "big")


def verify_w14() -> Fraction:
    W14 = [6, 1, 5, 1, 2, 7, 4, 5, 2, 1, 6, 5, 2, 1]
    out, _, _ = apply_word_fraction(Fraction(3, 5), W14)
    require(out > Fraction(609, 1000), "W14 output does not exceed 0.609")
    require(out.numerator.bit_length() == 72648, "W14 numerator bit length mismatch")
    require(out.denominator.bit_length() == 72648, "W14 denominator bit length mismatch")
    digest = sha256(
        minimal_unsigned_big_endian(out.numerator)
        + b"/"
        + minimal_unsigned_big_endian(out.denominator)
    ).hexdigest()
    require(
        digest == "6d203367282e4715b11351dd959dbae87ca967b5603bbc6dc82a10d2f596fa38",
        f"W14 SHA-256 mismatch: {digest}",
    )
    return out


# ---------------------------------------------------------------------------
# W25 exact rational interval certificate
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Interval:
    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @staticmethod
    def point(x: Fraction | int) -> "Interval":
        z = Fraction(x)
        return Interval(z, z)

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __mul__(self, other: "Interval") -> "Interval":
        values = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return Interval(min(values), max(values))

    def __truediv__(self, other: "Interval") -> "Interval":
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval denominator contains zero")
        values = (
            self.lo / other.lo,
            self.lo / other.hi,
            self.hi / other.lo,
            self.hi / other.hi,
        )
        return Interval(min(values), max(values))

    def round_outward(self, denominator: int) -> "Interval":
        if denominator <= 0:
            raise ValueError("rounding denominator must be positive")
        lo_scaled = self.lo * denominator
        hi_scaled = self.hi * denominator
        lo_int = lo_scaled.numerator // lo_scaled.denominator  # floor
        hi_int = -((-hi_scaled.numerator) // hi_scaled.denominator)  # ceil
        return Interval(Fraction(lo_int, denominator), Fraction(hi_int, denominator))

    def contains(self, q: Fraction) -> bool:
        return self.lo <= q <= self.hi


def verify_interval_arithmetic() -> None:
    """Sanity-check interval containment, including negative endpoints."""
    tests = [
        Interval(Fraction(-3, 2), Fraction(5, 4)),
        Interval(Fraction(1, 7), Fraction(5, 9)),
        Interval(Fraction(-11, 13), Fraction(-2, 17)),
    ]
    for a in tests:
        for b in tests:
            s = a + b
            p = a * b
            for xa in (a.lo, a.hi):
                for xb in (b.lo, b.hi):
                    require(s.contains(xa + xb), "interval addition lost an endpoint combination")
                    require(p.contains(xa * xb), "interval multiplication lost an endpoint combination")
            if not (b.lo <= 0 <= b.hi):
                q = a / b
                for xa in (a.lo, a.hi):
                    for xb in (b.lo, b.hi):
                        require(q.contains(xa / xb), "interval division lost an endpoint combination")
        rounded = a.round_outward(10 ** 6)
        require(rounded.lo <= a.lo and a.hi <= rounded.hi, "outward rounding lost original interval")


def check_interval(
    x: list[Interval],
    ell: int,
    grid_denominator: int,
) -> list[Interval]:
    require(ell in E_CANONICAL, f"unknown canonical check ell={ell}")
    e = E_CANONICAL[ell]
    one = Interval.point(1)
    two = Interval.point(2)
    den = one + x[ell] * x[ell]
    require(den.lo > 0, "interval check denominator is not strictly positive")
    y = [Interval.point(1)] + [Interval.point(0)] * 7
    for v in range(1, 8):
        if dot(v, e) == 0:
            num = x[v] * x[v] + x[v ^ ell] * x[v ^ ell]
        else:
            num = two * x[v] * x[v ^ ell]
        y[v] = (num / den).round_outward(grid_denominator)
    return y


def decimal_string(q: Fraction, digits: int, *, upper: bool = False) -> str:
    if digits < 0:
        raise ValueError("digits must be nonnegative")
    scale = 10 ** digits
    scaled = q * scale
    scaled_int = -((-scaled.numerator) // scaled.denominator) if upper else scaled.numerator // scaled.denominator
    sign = "-" if scaled_int < 0 else ""
    magnitude = abs(scaled_int)
    integer, fractional = divmod(magnitude, scale)
    if digits == 0:
        return f"{sign}{integer}"
    return f"{sign}{integer}.{fractional:0{digits}d}"


def verify_w25() -> tuple[Interval, Interval]:
    W25 = [6, 1, 5, 1, 2, 7, 4, 5, 2] + [1, 6, 5, 2] * 4
    F0 = Fraction(5997, 10000)
    grid = 10 ** 80
    x0 = (8 * F0 - 1) / 7
    x = [Interval.point(1)] + [Interval.point(x0)] * 7

    # A point trajectory in exact Fraction arithmetic becomes computationally
    # impractical after repeated squaring because numerator/denominator sizes
    # grow explosively. Rational outward interval truncation keeps endpoints
    # manageable while retaining a rigorous enclosure.
    acc_product = Interval.point(1)
    for ell in W25:
        A = (Interval.point(1) + x[ell] * x[ell]) / Interval.point(2)
        acc_product = (acc_product * A).round_outward(grid)
        x = check_interval(x, ell, grid)

    fidelity = Interval.point(0)
    for z in x:
        fidelity = fidelity + z
    fidelity = (fidelity / Interval.point(8)).round_outward(grid)

    lower_display = Fraction(
        60550725738718132769784559683258983421783519651344,
        10 ** 50,
    )
    upper_display = Fraction(
        60550725738718132769784559683258983421783519651346,
        10 ** 50,
    )
    require(lower_display < fidelity.lo, "W25 lower display bound is not strictly below certified interval")
    require(fidelity.hi < upper_display, "W25 upper display bound is not strictly above certified interval")

    Fstar_upper = Fraction(664887507913, 1099511627776)
    require(fidelity.lo > Fstar_upper, "W25 certificate does not enter the proven R4 basin")

    pooled_yield = acc_product / Interval.point(2 ** len(W25))
    require(pooled_yield.lo > 0, "W25 pooled-yield lower bound is not positive")
    require(Fraction(decimal_string(fidelity.lo, 80)) <= fidelity.lo, "displayed fidelity lower bound is not outward")
    require(Fraction(decimal_string(fidelity.hi, 80, upper=True)) >= fidelity.hi, "displayed fidelity upper bound is not outward")
    require(Fraction(decimal_string(pooled_yield.lo, 20)) <= pooled_yield.lo, "displayed yield lower bound is not outward")
    require(Fraction(decimal_string(pooled_yield.hi, 20, upper=True)) >= pooled_yield.hi, "displayed yield upper bound is not outward")
    return fidelity, pooled_yield


# ---------------------------------------------------------------------------
# High-precision effective post-check syndrome-noise checks (not exact certificates)
# ---------------------------------------------------------------------------


def fraction_to_mpf(q: Fraction) -> mp.mpf:
    return mp.mpf(q.numerator) / mp.mpf(q.denominator)


def check_noisy_mp(
    x: list[mp.mpf],
    ell: int,
    lam: mp.mpf,
) -> tuple[list[mp.mpf], mp.mpf]:
    require(ell in E_CANONICAL, f"unknown canonical check ell={ell}")
    # Numerical root finders may probe a tiny analytic continuation outside
    # the physical lambda interval. Final reported roots/crossings are checked
    # to lie in the physical region.
    e = E_CANONICAL[ell]
    den = 1 + x[ell] ** 2
    require(den > 0, "noisy check denominator is nonpositive")
    kappa = 1 - mp.mpf(8) * lam / 7
    y = [mp.mpf(1)] + [mp.mpf(0)] * 7
    for vv in range(1, 8):
        if dot(vv, e) == 0:
            ideal = (x[vv] ** 2 + x[vv ^ ell] ** 2) / den
        else:
            ideal = 2 * x[vv] * x[vv ^ ell] / den
        y[vv] = kappa * ideal
    return y, den / 2


def word_noisy_mp(Fv: mp.mpf, lam: mp.mpf, word: list[int]) -> mp.mpf:
    _validate_check_word(word)
    # Root finders may evaluate the analytic continuation slightly outside [0,1].
    # Physical fixed points are validated after convergence.
    x0 = (8 * Fv - 1) / 7
    x = [mp.mpf(1)] + [x0] * 7
    for ell in word:
        x, _ = check_noisy_mp(x, ell, lam)
    return mp.fsum(x) / 8


def r4_noisy_mp(Fv: mp.mpf, lam: mp.mpf) -> mp.mpf:
    return word_noisy_mp(Fv, lam, [1, 2, 4, 7])


def noisy_roots(lam: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    f = lambda z: r4_noisy_mp(z, lam) - z
    lower_guess = mp.mpf("0.6047116657215318756") + 2 * lam
    upper_guess = 1 - mp.mpf("2.43") * lam
    lower = mp.findroot(
        f,
        (lower_guess - mp.mpf("0.01"), lower_guess + mp.mpf("0.01")),
    )
    if lam == 0:
        upper = mp.mpf(1)
    else:
        upper = mp.findroot(
            f,
            (upper_guess - mp.mpf("0.02"), upper_guess + mp.mpf("0.02")),
        )
    residual_tol = mp.mpf("1e-45")
    require(abs(f(lower)) < residual_tol, f"lower fixed-point residual too large at lambda={lam}")
    require(abs(f(upper)) < residual_tol, f"upper fixed-point residual too large at lambda={lam}")
    require(mp.mpf(0) <= lower <= upper <= mp.mpf(1), "fixed points outside physical interval or misordered")
    return lower, upper


def verify_noisy_recurrence_numbers(fseed_representative: Fraction) -> None:
    mp.mp.dps = 80
    expected = {
        "0.0001": ("0.604885431446266501", "0.999757042328097016"),
        "0.001": ("0.606459520811496560", "0.997561265472276269"),
        "0.005": ("0.613690033554464149", "0.987589998096974921"),
        "0.01": ("0.623332228826214012", "0.974571486065393376"),
        "0.02": ("0.645211515251388839", "0.946095042634381636"),
    }
    tol = mp.mpf("5e-16")
    for lam_s, (lo_s, hi_s) in expected.items():
        lam = mp.mpf(lam_s)
        require(mp.mpf(0) <= lam <= mp.mpf(1), "tabulated lambda outside physical range")
        lo, hi = noisy_roots(lam)
        require(abs(lo - mp.mpf(lo_s)) < tol, f"lower noisy root mismatch at lambda={lam_s}")
        require(abs(hi - mp.mpf(hi_s)) < tol, f"upper noisy root mismatch at lambda={lam_s}")
        d_lo = mp.diff(lambda z: r4_noisy_mp(z, lam), lo)
        d_hi = mp.diff(lambda z: r4_noisy_mp(z, lam), hi)
        require(abs(d_lo) > 1, f"lower fixed point is not numerically unstable at lambda={lam_s}")
        require(abs(d_hi) < 1, f"upper fixed point is not numerically attracting at lambda={lam_s}")

    Fstar = mp.mpf("0.6047116657215318756")
    dlam = mp.diff(lambda ll: r4_noisy_mp(Fstar, ll), mp.mpf(0))
    dF = mp.diff(lambda zz: r4_noisy_mp(zz, mp.mpf(0)), Fstar)
    shift = -dlam / (dF - 1)
    require(abs(shift - mp.mpf("1.7365360404235246271")) < mp.mpf("2e-18"), "unstable-root shift mismatch")

    eq1 = lambda FF, ll: r4_noisy_mp(FF, ll) - FF
    eq2 = lambda FF, ll: mp.diff(lambda zz: r4_noisy_mp(zz, ll), FF) - 1
    Fc, lc = mp.findroot((eq1, eq2), (mp.mpf("0.787"), mp.mpf("0.0475")))
    require(abs(Fc - mp.mpf("0.7870464318037645000")) < mp.mpf("2e-18"), "saddle-node fidelity mismatch")
    require(abs(lc - mp.mpf("0.04748530115982306916")) < mp.mpf("2e-18"), "saddle-node lambda mismatch")
    require(abs(eq1(Fc, lc)) < mp.mpf("1e-45"), "saddle-node fixed-point residual too large")
    require(abs(eq2(Fc, lc)) < mp.mpf("1e-40"), "saddle-node tangency residual too large")

    W14 = [6, 1, 5, 1, 2, 7, 4, 5, 2, 1, 6, 5, 2, 1]
    W25 = [6, 1, 5, 1, 2, 7, 4, 5, 2] + [1, 6, 5, 2] * 4

    def lower_root(ll: mp.mpf) -> mp.mpf:
        return noisy_roots(ll)[0]

    h25 = lambda ll: word_noisy_mp(mp.mpf("0.5997"), ll, W25) - lower_root(ll)
    h14 = lambda ll: word_noisy_mp(mp.mpf("0.6"), ll, W14) - lower_root(ll)
    l25 = mp.findroot(h25, (mp.mpf("0.8e-7"), mp.mpf("1.2e-7")))
    l14 = mp.findroot(h14, (mp.mpf("2e-5"), mp.mpf("3e-5")))
    require(mp.mpf(0) < l25 < mp.mpf(1), "W25 crossing outside physical lambda range")
    require(abs(l25 - mp.mpf("1.0070215932877671555e-7")) < mp.mpf("2e-22"), "W25 post-check syndrome-noise crossing mismatch")
    require(mp.mpf(0) < l14 < mp.mpf(1), "W14 crossing outside physical lambda range")
    require(abs(l14 - mp.mpf("2.5168632113128632438e-5")) < mp.mpf("2e-20"), "W14 post-check syndrome-noise crossing mismatch")
    require(abs(h25(l25)) < mp.mpf("1e-40"), "W25 crossing residual too large")
    require(abs(h14(l14)) < mp.mpf("1e-40"), "W14 crossing residual too large")

    # The W25 midpoint is only a high-precision numerical representative.
    # The final faulty-tree calculation is conditional: the factor (1-t) Z R
    # requires a classical fault flag independent of the ideal trajectory, or
    # equivalently a mixture containing the ideal instrument with weight >= 1-t.
    # Converting the tree probability to eta also assumes 15 iid locations.
    Fseed = fraction_to_mpf(fseed_representative)
    x0 = (8 * Fseed - 1) / 7
    x = [mp.mpf(1)] + [x0] * 7
    acceptances: list[mp.mpf] = []
    for ell in [1, 2, 4, 7]:
        x, A = check_noisy_mp(x, ell, mp.mpf(0))
        acceptances.append(A)
    R = mp.fsum(x) / 8
    Z = (
        acceptances[0] ** 8
        * acceptances[1] ** 4
        * acceptances[2] ** 2
        * acceptances[3]
    )
    tmax = Z * (R - Fseed) / (Fseed + Z * (R - Fseed))
    eta = 1 - (1 - tmax) ** (mp.mpf(1) / 15)
    require(abs(R - mp.mpf("0.60757261096713259305")) < mp.mpf("2e-20"), "conditional faulty-tree R value mismatch")
    require(abs(Z - mp.mpf("0.0008060439641114767003")) < mp.mpf("2e-22"), "conditional faulty-tree Z value mismatch")
    require(abs(tmax - mp.mpf("2.7493662371752501e-6")) < mp.mpf("2e-20"), "conditional faulty-tree probability value mismatch")
    require(abs(eta - mp.mpf("1.8329131764811310e-7")) < mp.mpf("2e-21"), "conditional iid-location fault value mismatch")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _emit(message: str, quiet: bool) -> None:
    if not quiet:
        print(message)


def run_verification(*, exact_only: bool = False, quiet: bool = False) -> tuple[Interval, Interval]:
    verify_environment()
    _emit(f"[env] {environment_string()}", quiet)
    if not __debug__:
        _emit("[env] Python optimized mode (-O) detected; explicit verifier checks remain active", quiet)

    verify_basic_f2_and_labels()
    _emit("[ok] F_2^3 encoding and Majorana label table", quiet)

    verify_affine_structure()
    _emit("[ok] affine group, normalizer, and all 28 conjugate checks", quiet)

    verify_fourier_formula()
    _emit("[ok] generalized Fourier formula; rational regression + exhaustive basis verification", quiet)

    verify_round_maps()
    _emit("[ok] R3/R4 factorization, root structure, exact sample signs, and Taylor laws", quiet)

    verify_fault_table()
    _emit("[ok] exact low-weight leaf-fault counts for R3 and R4", quiet)

    verify_readout_error_map()
    verify_common_mode_identity()
    verify_effective_syndrome_shift_table()
    _emit("[ok] readout-map consequences, common-mode identity, and effective syndrome-shift table", quiet)

    guard_cases = verify_invalid_input_guards()
    _emit(f"[ok] invalid-input safeguard suite ({guard_cases} malformed cases rejected)", quiet)

    verify_interval_arithmetic()
    _emit("[ok] rational interval arithmetic sanity checks", quiet)

    w14 = verify_w14()
    _emit("[ok] W14 exact Fraction certificate; SHA-256 regression/reproducibility check", quiet)

    fidelity, pooled_yield = verify_w25()
    _emit("[ok] W25 exact rational interval certificate", quiet)
    _emit("     fidelity lower = " + decimal_string(fidelity.lo, 80), quiet)
    _emit("     fidelity upper = " + decimal_string(fidelity.hi, 80, upper=True), quiet)
    _emit("     pooled-yield lower = " + decimal_string(pooled_yield.lo, 20), quiet)
    _emit("     pooled-yield upper = " + decimal_string(pooled_yield.hi, 20, upper=True), quiet)

    if not exact_only:
        # The midpoint is used only as a numerical representative of the
        # rigorously certified W25 interval. It is not asserted to be the
        # exact W25 output. The interval width is far below the precision of
        # the numerical quantities checked in the post-check syndrome-noise section.
        fseed_representative = (fidelity.lo + fidelity.hi) / 2
        verify_noisy_recurrence_numbers(fseed_representative)
        _emit("[ok] high-precision effective post-check syndrome-noise fixed points and robustness numbers", quiet)
    else:
        _emit("[skip] high-precision numerical post-check syndrome-noise checks (--exact-only)", quiet)

    # Keep w14 live to make accidental removal of the certificate conspicuous to
    # static readers and to guarantee that it was computed before success.
    require(w14 > 0, "internal W14 liveness check failed")
    _emit("ALL CHECKS PASSED", quiet)
    return fidelity, pooled_yield


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Independent verifier for the Majorana non-Gaussian ancilla purification notes. "
            "The full mode runs exact certificates and high-precision numerical observations."
        )
    )
    parser.add_argument(
        "--exact-only",
        action="store_true",
        help="run the exact/certificate routines only; skip the mpmath-based post-check syndrome-noise section",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="suppress progress output; failures still raise a non-zero exit status",
    )
    parser.add_argument(
        "--show-environment",
        action="store_true",
        help="print interpreter/package/platform versions and exit",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        verify_environment()
        if args.show_environment:
            print(environment_string())
            return 0
        run_verification(exact_only=args.exact_only, quiet=args.quiet)
        return 0
    except (VerificationError, ValueError, ZeroDivisionError, ArithmeticError) as exc:
        print(f"VERIFICATION FAILED: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
