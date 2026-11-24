# Recheck against the current preprint

## Verdict

**PASS: the supplied `majorana_verification.zip` covers the computational
claims in the current preprint's Data and code availability paragraph.**
No mathematical algorithm, certificate endpoint, numerical result, or Python
source needed changing for that purpose.

This is a recheck of the archive actually uploaded with this request, not an
inference from earlier conversation summaries. The matching manuscript is the
uploaded text source `Wklejony tekst(10).txt`.

- Manuscript title: *Reducing the cost of ancilla purification for Majorana quantum computation*
- Uploaded manuscript SHA-256: `5b16db2ee6f1e43d419e4452796a0b32df7879e7362589c01e8941233e62f355`
- Original uploaded ZIP SHA-256: `db7926bbb6917a53f37f08f470f6dc0914df819c21874f887c2612b82320a705`
- Recheck date: 2026-11-24

## Executions performed

The documented full command passed:
```text
python run_all.py --output-dir <fresh-directory> --optimized-guard-check
```

It was then rerun in a clean virtual environment exposing only the already
installed pinned SymPy and mpmath packages (copied into the environment by local
package-directory links; no package download was performed). The clean run passed
with no script stderr. The full runner was also exercised with `PYTHONOPTIMIZE=1`,
so all child scripts ran in optimized mode, not only the extended verifier.

Separate normal and `-O` executions of the extended verifier's `--exact-only`
mode passed as well.

Environment exercised: Python 3.13.5, SymPy 1.14.0, mpmath 1.3.0, Linux.
This recheck does not certify every Python/platform version supported by the
README, and it is not a test of a new network installation of the dependencies.

The seven scientific scripts completed successfully:
1. `verify_modelB_gate_noise_extended.py` (v1.2)
2. `independent_review.py`
3. `verify_heterogeneous_bootstrap.py`
4. `verify_h12_recheck.py`
5. `verify_one_stage_control.py`
6. `verify_even_sector_twirl.py`
7. `verify_compare.py`

All nine entries in the original `CHECKSUMS.sha256` matched. An AST scan found no
disableable Python `assert` statements in any of the eight Python files.

The regenerated comparison JSON and both stage CSVs match the bundled reference
files byte-for-byte. The independent-recurrence and H12 recheck JSON files match
in all scientific fields (environment metadata may differ).

## Coverage of the requested statement

| Required content | Present and rerun | Qualification |
|---|---|---|
| Integer and rational certificates | Extended, heterogeneous, and integer-weight recheck scripts | Exact arithmetic for the fixed noiseless constructions |
| Independent recurrence checks | `independent_review.py` and integer-polynomial H12 recheck | Direct syndrome enumeration independent of the Fourier recursion |
| Controlled fourth-stage comparison | `verify_one_stage_control.py` | Changes only the fresh-input stage to self-pairing, leaving the word fixed |
| Numerical resource calculations | `verify_compare.py` | Specified pooled/retry schedules, not optimized factory costs |
| Effective-noise calculations | Extended, heterogeneous, H12 recheck, and comparison scripts | High-precision numerical roots and derivatives; not certified global noisy basins |
| Dependencies and reproducibility commands | README, pinned requirements, runner | All necessary Python files are present; no private data inputs |
| General even-input corollary | `verify_even_sector_twirl.py` | Exact dephasing/twirl algebra; Gaussian implementation remains an analytical argument |

## Alignment with the current manuscript

Eighteen additional source-alignment check groups passed, including:

- all twelve H12 stage rows and their displayed acceptance probabilities;
- all seventeen R4 accepted-fault enumerator rows and the sixteen g15 coefficients;
- the H12 endpoint enclosures at 0.5995 and 0.5997;
- the 25520/25521-bit rational fingerprint at 0.5995;
- exact retry costs and yield enclosures;
- the common-input yield ratio and its printed enclosing interval;
- both controlled fourth-stage comparisons;
- all six fixed-accuracy cost-table rows, stopping rounds, and final errors;
- both near-endpoint complete-factory examples;
- all four displayed effective-noise local crossings;
- the 2304-leaf / 2303-node / 9212-measurement accounting;
- the conditioned depolarizing conversion to 1602/2551.

All 33 equation/table/corollary labels referenced in `CLAIM_COVERAGE.md` still
exist in the uploaded manuscript.

Representative reproduced values:

```text
H12, raw F=0.5995:
  0.605026344233701236200599895568 < Fout
  Fout < 0.605026344233701236200599895569
  1051931.361812434 < expected raw-copy cost < 1051931.361812435

Fourth-stage all-self-pairing control, raw F=0.5995:
  Fout = 0.551387370616520...  (below the R4 basin)

Same-input yield ratio, raw F=0.5997:
  H12 yield / W25 yield = 16968397.494075350489047... .

Numerical local operation-noise crossings:
  H12, F=0.5997: 1.29318072951395676241...e-4
  H12, F=0.5995: 5.65511716055257587648...e-6
```

The ratio above uses the freshly recomputed narrow W25 interval. The wider
printed ratio interval in the manuscript is also correct.

The exact and numerical evidence must not be conflated. Integer/Fraction
evaluation and outward rational intervals certify the stated noiseless
endpoints. The noise crossings remain numerical observations even when stable
under precision changes and reproduced in two coordinate systems.

## Changes made in the refreshed distribution

**No changes were made to any `.py` file or to `requirements.txt`.**
Every such file is byte-for-byte identical to the original uploaded archive.

Only packaging and documentation were refreshed:
- README instructions now say which directory to run from and identify the
  current manuscript by its audit record rather than an ambiguous draft date.
- The README includes the intended repository address and initial-upload guidance.
- `.gitignore` excludes generated local outputs and caches.
- `recheck_results/` contains fresh execution and manuscript-alignment evidence.
- `DATA_AVAILABILITY.tex` supplies proposed public-repository wording, to be used
  only after the files have actually been uploaded.
- `CHECKSUMS.sha256` now covers all distributed files except itself.

The older `AUDIT_REPORT.md` is retained as historical documentation; the present
file records the new comparison with the latest uploaded manuscript.

## What should be changed before submission

The uploaded manuscript still says that the ZIP is available from the authors
upon request. Once the repository contains these files, replace that phrase with
the actual public-repository location. Do not claim an already archived release
or DOI until one exists.

For the repository layout, upload the **contents** of the extracted
`majorana_verification/` folder at repository root. The ZIP may additionally be
provided as the supplementary archive with the exact filename
`majorana_verification.zip`.

Agree the software licence with the author group. This packaging refresh does
not select a licence, assign copyright, change authorship, or modify the preprint.

For verification runs, use a fresh output directory. The unmodified runner has
a 240-second timeout per child script and does not remove old summaries when
reusing a directory. A nonzero process exit must not be overridden by an old
PASS file. Individual script commands can diagnose a timeout on slower hardware.

## Limits

This is an audit of the released computational evidence and its correspondence
to the manuscript. It is not a new proof of every analytical statement.

The archive does not independently simulate Bravyi's complete braid-and-fusion
instrument, prove the imported four-mode concurrence theorem, optimize over all
Gaussian protocols, or estimate a full device's hardware threshold. Those
limits are stated in the README and coverage matrix and are compatible with the
quoted code-availability claim.

The search procedure used to discover the fixed schedules is not included.
It is not needed to reproduce the fixed constructions or certificates claimed
in this manuscript.

Neither the repository nor the manuscript was uploaded or modified remotely
during this recheck.
