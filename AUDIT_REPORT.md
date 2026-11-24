# Scientific code audit — 27 October 2026

## Verdict

**PASS for the specifically claimed code-backed results of the paper**, under the exact syndrome-model assumptions written in the manuscript. The supplied verifiers and the new release reproducer ran successfully under Python 3.13.5 / SymPy 1.14.0 / mpmath 1.3.0. This does not constitute a full independent simulation of the Majorana circuit, or a proof that the schedules are optimal.

## What was tested

1. **v1.2 extended verifier:** tested full and `--exact-only`; both passed. Also tested the full verifier under `python -O`; it still passed all checks, including the 26 intentional malformed-input cases. The earlier original 793-line user verifier, supplied separately and byte-identical to the original bundled in the user ZIP, also passed in ordinary and `python -O` modes. The v1.2 safeguards address the risk of Python disabling `assert` statements.
2. **Independent exact consistency:** the probability-space reconstruction of R3 and R4 agrees with the Fourier recursion, including fault counts, symbolic factorization, Sturm root counting, Taylor coefficients and the W25 rational interval certificate. W14 and H12 rational regression hashes match the manuscript.
3. **H12 endpoints and costs:** unnormalized integer-weight recursion agrees exactly with normalized probability and Fourier implementations; the closed cost expression agrees with forward retries and reverse demand propagation; stage-by-stage records were generated.
4. **Controlled fresh-input comparison:** recomputed H12 and the all-self-pairing control at raw F = 0.5995, 0.5997 and 0.6000 using exact integer weights. At 0.5995 the respective fidelities are 0.605026344233701... and 0.551387370616520...; the expected raw-copy costs are approximately 1.05193 million and 2.32654 million.
5. **Even-parity twirl corollary:** exhaustive exact algebra check of all 64 matrix units in the eight-dimensional even sector, the eight dephasing characters, and the 168 linear frame permutations. This supports, but does not replace, the manuscript's proof of Gaussian implementability.
6. **Numerical effective noise:** noisy fixed points, fold, four local bootstrap noise crossings, one-post-check-fault counts, and conditional faulty-tree numerical bounds agree across implementations. Crossing values are explicitly treated as numerical and local.
7. **Matched end-to-end costs:** checked raw F = 0.65 and 0.9; added full-pipeline endpoint F = 0.5995 (eight rounds after prefix, 2.77516e20 raw copies approximately) and 0.5997 (five rounds after prefix, 4.52768e14 approximately), with extra 80/110-digit crosschecks.

## Changes to the code for this release

- Kept the user-supplied **extended verifier v1.2** and **heterogeneous verifier** unchanged.
- Retained the independent review and full H12 recheck from the previously prepared reproducibility set.
- Replaced the `assert` checks in the one-stage controlled comparison with explicit failures that remain active under `python -O`; added a direct cost ordering check.
- Added `--output` to the comparison script so it writes its JSON to a chosen location without modifying the released source folder, and extended it to include both near-threshold end-to-end factory costs.
- Added an exact `GL(3,2)` / dephasing algebra checker for the arbitrary-even-state corollary.
- Added `run_all.py`, pinned `requirements.txt`, per-script logs, a coverage matrix, output consistency checks, and a machine-readable summary with SHA-256 source hashes.

## Boundaries of the audit

The certificate for a *particular* H12 schedule is exact. Discovery of an optimal schedule, mathematical optimality of 0.5995, or exhaustive search through all possible Gaussian protocols is not established.

The code verifies the effective syndrome transformation taken from the original published Bravyi comparison. It does **not** independently compile and verify all the microscopic braid, parity-readout, and feed-forward Kraus operators, nor does it simulate device routing, state storage, leakage, or noisy hardware. An accepted state with arbitrary off-diagonal coherences can be treated by the manuscript's free Gaussian twirl proof, with finite-dimensional algebra checked here, rather than by extending every syndrome recursion directly to arbitrary coherences.

The four-mode convex-Gaussian concurrence and FLO no-go argument invoke established mathematical statements and an analytic proof in the manuscript; they are not consequences of this executable archive alone. Likewise, the asymptotic resource exponent follows analytically from the exact local recurrence, not from a finite numerical fit.

The effective operation-noise roots/fold and the low-fidelity crossings are high-precision numerical observations and are not given interval-certified global uniqueness or microscopic hardware-threshold interpretations.

## Reproducibility status

The archive is self-contained **as a code supplement**, with no need to fetch additional scripts. The package does not contain the original exploratory search logs nor the literature sources; those are not needed to reproduce its theorem-level numerical certificates. It should be made publicly available by the authors or submitted as supplementary material before describing it as an archived public deposit.
