# Majorana ancilla purification — executable verification package

Companion code for **“Reducing the cost of ancilla purification for Majorana quantum computation”**.

Repository: https://github.com/chylonczyk/Reducing-the-cost-of-ancilla-purification-for-Majorana-quantum-computation

This archive reproduces the *syndrome-model* mathematical results and numerical comparisons described in the manuscript. The current uploaded preprint was checked against this package; its source hash and the new test record are in [RECHECK_CURRENT_PREPRINT.md](RECHECK_CURRENT_PREPRINT.md). The eight Python sources, including the runner, and the pinned requirements are unchanged from the supplied archive.

No private files, external data, network access, or hardware access are required to run the calculations after installing the dependencies. Installing the dependencies may require network access.

## Reproduce all claimed computational checks

After extracting the ZIP, change into its `majorana_verification/` directory.
When installing from the GitHub repository, run from the repository root,
where `README.md`, `run_all.py`, and `requirements.txt` are located.


Requires **Python 3.10+**. The supplied archive and this recheck were exercised with Python **3.13.5**, SymPy **1.14.0**, and mpmath **1.3.0** on Linux. Other supported Python versions were not exercised in this recheck.

```bash
python -m pip install -r requirements.txt
python run_all.py --output-dir verification_results --optimized-guard-check
```

An exit code of **0**, followed by `ALL REPRODUCIBILITY CHECKS PASSED`, is the intended success condition. The command writes one log per script, `verification_summary.json`, the stage-by-stage H12 CSV data, and supporting numerical JSON files under `verification_results/`. A failed scientific assertion or inconsistent reproducibility result exits nonzero. All scientific verification conditions in the primary verifier remain active under Python's `-O` optimized mode.

There is also a faster way to run the exact checks of the original balanced recurrence *alone*, without its numerical noise continuation:

```bash
python verify_modelB_gate_noise_extended.py --exact-only
```

For individual scripts:

```bash
python verify_modelB_gate_noise_extended.py
python independent_review.py --json verification_results/independent_review.json
python verify_heterogeneous_bootstrap.py
python verify_h12_recheck.py --output verification_results/h12_recheck --reference verify_heterogeneous_bootstrap.py
python verify_one_stage_control.py
python verify_even_sector_twirl.py
python verify_compare.py --output verification_results/comparison_results.json
```

Before passing paths to `independent_review.py --json`, create their parent directory if necessary (`mkdir -p verification_results`); the all-in-one runner does this automatically. Running `verify_heterogeneous_bootstrap.py` directly requires `mpmath` for the optional numerical section; the exact certificate uses `fractions.Fraction`. The all-in-one runner *requires* the numerical section and will not silently count a skipped numerical check as a full pass.

The supplied **`reference_results/`** gives an example of expected JSON/CSV outputs from a passing run; run the scripts yourself to obtain a new independent record on your machine.

## Files and provenance

| Script | Role | Method |
|---|---|---|
| `verify_modelB_gate_noise_extended.py` | Original extended verifier **v1.2**, including 26 malformed-input safeguard tests, affine frame checks, exact R3/R4 polynomial checks, low-weight fault counts, W14/W25 certificates, readout checks, and numerical syndrome-noise roots | Exact `Fraction`/SymPy and rational intervals; separate mpmath numerical tests |
| `independent_review.py` | Independently reconstructs accepted syndrome enumerators, recurrence factorization/Sturm counts, and W25 probability-space interval certificate; checks an adversarial faulty-branch counterexample | Exact integer/rational and symbolic |
| `verify_heterogeneous_bootstrap.py` | Computes H12 in independently normalized probability and Fourier coordinates; checks exact endpoint at 0.5995; evaluates local noisy crossings | Exact `Fraction`; numerical mpmath |
| `verify_h12_recheck.py` | Independent *unnormalized integer-weight* H12 recheck, stage CSVs, cross-implementation agreement at every stage, expected retry and reverse-demand costs, SHA-256 fingerprint, and numerical noise crossings | Exact integer/Fraction; separate numerical mpmath |
| `verify_one_stage_control.py` | H12 versus otherwise identical all-self-pairing control (replace only `R_7` with `S_7`) at F=0.5995, 0.5997, 0.6000 | Exact integer/Fraction; replaced disableable `assert` with always-active checks |
| `verify_even_sector_twirl.py` | Tests stabilizer dephasing and `GL(3,2)` averaging on **all 64 matrix-unit basis elements** in the eight-dimensional even sector | Exact `Fraction`, standard library |
| `verify_compare.py` | Independent pooled/retry finite-accuracy factory costs, R3/R4 one-operation-fault counts, and complete H12+R4 costs at both low-fidelity endpoints | Numerically high precision `mpmath`; cross-precision checks |
| `run_all.py` | Executes scripts, checks that expected tests *actually ran*, compares their JSON outputs against stated targets, records logs and SHA-256 source hashes | Python standard library, strict nonzero exit on failure |

The scientific-verification sources include the extended verifier v1.2, the normalized Fourier/probability heterogeneous verifier, the unnormalized-integer H12 recheck, and separate recurrence, control, twirl, and comparison scripts. All eight Python files, including `run_all.py`, are byte-for-byte unchanged in this packaging refresh. The earlier development audit is retained as [AUDIT_REPORT.md](AUDIT_REPORT.md); the new source-alignment and execution record is [RECHECK_CURRENT_PREPRINT.md](RECHECK_CURRENT_PREPRINT.md).

The source code used to *discover* the words W25/H12 is not included. Verifying the fixed constructions in the paper does **not** require re-running that search.

Some messages in `independent_review.py` refer to historical issues (the missing ideal-branch assumption, a missing zero fixed point, and finite-batch versus asymptotic yield). These are intentional regression/counterexample checks. Their successful execution does not mean that the corrected manuscript still makes those mistakes.

## Results and mathematical scope

Please consult [CLAIM_COVERAGE.md](CLAIM_COVERAGE.md) for a claim-by-claim map to individual scripts, and [AUDIT_REPORT.md](AUDIT_REPORT.md) for the final code audit.

- The R3/R4 recurrence, low-weight accepted-fault enumerators, exact factorization and rational W14/W25/H12 endpoint certificates are checked **exactly** within the stated syndrome model.
- The 3-, 4-, 5-, and 8-round resource-count claims are **numerical** calculations of **prescribed** sequences, with independent retries/pooling and free isotropization. These are neither global optimality results nor circuit-level hardware costs.
- The effective post-check syndrome-noise fixed points, folds and local crossing values are **numerical**, without rigorous global interval enclosures. One-fault combinatorics (including the coefficient `17/7`) is exact within the syndrome-fault model.
- The arbitrary even-parity input corollary follows from a free Gaussian stabilizer dephasing/twirl *as proved in the manuscript*. The new checker verifies the finite-dimensional algebra on a basis; it does **not** implement the complete physical braid circuit.
- The physical Bravyi comparison instrument, Gaussian compilation into a particular Majorana hardware architecture, full parity preparation, and the four-mode convex-Gaussian no-go theorem are analytic or imported physical/theoretical ingredients, **not** independently derived by these numerical scripts.
- The cost model assumes independent raw blocks and independent retries, perfect storage, ideal parity-sector selection (if required), and the paper's stated syndrome-level operation-noise model. No noisy storage, correlated hardware faults, spatial routing, or complete fault-tolerant device is simulated.

The known mathematical limitation is intentional: passing scripts certifies the model-level computations and checks algebraic identities. It is not an independent proof of the microscopic braid-and-fusion instrument or a demonstration of operational feasibility.

## Publishing and code availability

For the initial GitHub upload, copy the **contents** of the extracted
`majorana_verification/` directory to the repository root, rather than uploading only
an opaque ZIP. This makes the README and source files immediately visible.

The archive can also be distributed as `majorana_verification.zip`, as named in the
manuscript. Do not claim public availability until the files have actually been
uploaded. [DATA_AVAILABILITY.tex](DATA_AVAILABILITY.tex) gives replacement text to use
after publication of the repository contents.

For a fixed submission snapshot, retain the archive and its checksums, and record the
actual commit or release identifier when one has been created. No release identifier,
DOI, publication status, or author approval has been invented in this package.

No software licence has been selected in this refresh. The authors should agree on a
licence and add its file before the public release; no copyright ownership has been
assigned here. The manuscript author list and manuscript text have not been modified.

## Integrity and working directories

`CHECKSUMS.sha256` covers all distributed files except the checksum list itself.
On systems with `sha256sum`, run `sha256sum -c CHECKSUMS.sha256` before executing the
scripts. The runner records fresh source hashes in its output summary.

Use a new output directory for each run. The unmodified runner does not clear old
outputs automatically. Its nonzero exit code, not the mere presence of an older
`verification_summary.json`, determines failure.

The runner has a 240-second timeout per child script. On an unusually slow machine,
a timeout is not a failed mathematical certificate: run the individual commands
above to distinguish a timeout from a failed scientific check.

Generated outputs and Python caches are excluded by `.gitignore`. Tracked
`reference_results/` are examples from a passing run, not inputs to the certificate
calculations. Current recheck evidence is recorded separately under `recheck_results/`.
