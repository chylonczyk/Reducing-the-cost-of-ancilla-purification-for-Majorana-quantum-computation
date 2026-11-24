#!/usr/bin/env python3
"""Run the complete reproducibility audit from a clean working directory.

A nonzero exit status indicates a failed check. Logs and JSON/CSV outputs
are written to --output-dir. No network access or experimental data needed.
"""
from __future__ import annotations
import argparse
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SOURCES = (
    'verify_modelB_gate_noise_extended.py',
    'independent_review.py',
    'verify_heterogeneous_bootstrap.py',
    'verify_h12_recheck.py',
    'verify_one_stage_control.py',
    'verify_even_sector_twirl.py',
    'verify_compare.py',
)


def demand(test, why):
    if not test:
        raise RuntimeError(f'REPRODUCIBILITY CHECK FAILED: {why}')


def run_script(output, label, args, must_contain, timeout=240):
    command = [sys.executable, str(ROOT/label), *map(str, args)]
    print('RUN', ' '.join(command), flush=True)
    cp = subprocess.run(command, cwd=output, capture_output=True, text=True,
                        encoding='utf-8', errors='replace', timeout=timeout,
                        check=False)
    name = Path(label).stem
    log = output / f'{name}.log'
    log.write_text('COMMAND: ' + ' '.join(command) + '\n\n'
                   + cp.stdout + '\n[STDERR]\n' + cp.stderr,
                   encoding='utf-8')
    demand(cp.returncode == 0, f'{label} exit code {cp.returncode}; see {log}')
    for token in must_contain:
        demand(token in cp.stdout, f'{label} did not emit {token!r}; see {log}')
    demand('[skip]' not in cp.stdout, f'{label} skipped verification; see {log}')
    print('PASS', label, flush=True)
    return cp.stdout


def verify_results(out):
    r = json.loads((out/'h12_recheck'/'recheck_results.json').read_text())
    c = json.loads((out/'comparison_results.json').read_text())
    demand(r['basis_pairs_verified'] == 1792, 'H12 independent check count')
    demand(r['R4']['Sturm_variations'] == [8, 7], 'R4 Sturm sign variations')
    from fractions import Fraction as F
    fstar_hi = F(664887507913, 1099511627776)
    bounds = r['exact']['0.5995']['Fout']
    demand(F(bounds[0]) > fstar_hi, 'exact H12 lower output does not exceed exact R4 upper bound')
    demand(F(bounds[0]) < F(bounds[1]), 'H12 interval reversed')
    demand(r['exact']['0.5995']['leaves'] == 2304, 'H12 accepted-tree leaves')
    demand(r['exact']['0.5995']['nodes'] == 2303, 'H12 accepted-tree nodes')
    cost = r['exact']['0.5995']['cost']
    demand(F(cost[0]) <= F('1051931.3618124345') <= F(cost[1]), 'H12 cost enclosure')
    ratio = r['yield_ratio_from_recomputed_W25_interval']
    demand(Decimal(ratio['lower']) > Decimal('16968300')
           and Decimal(ratio['upper']) < Decimal('16968500'), 'same-F yield ratio')
    for key in ['H12_5997', 'H12_5995', 'W25_5997', 'W14_6000']:
        demand(key in r['noise_numerical'], f'missing noisy local crossing {key}')
        demand(r['noise_numerical'][key]['signs_at_90_110_percent'] == ['+', '-'],
               f'missing crossing sign bracket {key}')
    for family in ('R3','R4'):
        demand(c['single_postcheck_faults'][family]['first_order_infidelity_coefficient'] == '17/7',
               f'first-order noise coefficient {family}')
    expected_rounds = {'.65':{'R3':6,'R4':4,'H12_then_R4':1},
                       '.9':{'R3':3,'R4':2,'H12_then_R4':0}}
    for raw, variants in expected_rounds.items():
        for k,n in variants.items():
            row = c['finite_accuracy_numerical'][raw][k]
            demand(row['rounds_after_prefix'] == n, f'finite-accuracy stopping round: {raw}, {k}')
            demand(Decimal(row['final_infidelity']) <= Decimal('1e-6'),
                   f'finite-accuracy target failure: {raw}, {k}')
    for raw, rounds, target_cost in (('.5995',8,Decimal('2.77516e20')),
                                     ('.5997',5,Decimal('4.52768e14'))):
        row=c['endpoint_factories_numerical'][raw]
        demand(row['rounds_after_prefix']==rounds, f'near-endpoint stopping for {raw}')
        demand(Decimal(row['final_infidelity'])<Decimal('1e-6'),
               f'near-endpoint accuracy for {raw}')
        demand(abs(Decimal(row['expected_raw_copies'])/target_cost-1) < Decimal('2e-5'),
               f'near-endpoint cost for {raw}')
    return dict(
        H12_exact_fidelity_interval=bounds,
        H12_expected_copy_cost_interval=cost,
        W25_vs_H12_yield_ratio_interval=ratio,
        first_order_noisy_coefficient='17/7 for both R3 and R4',
        checked_exact_fidelity_above_Fstar_upper=True,
        noise_crossings='four local numerical crossings checked, NOT interval-certified global thresholds',
        finite_accuracy_rounds=expected_rounds,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', default='verification_results', type=Path,
                    help='where to write logs and machine-readable checks')
    ap.add_argument('--optimized-guard-check', action='store_true',
                    help='also run the main v1.2 verifier with python -O')
    args = ap.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    try:
        import sympy, mpmath
        versions = {'python':platform.python_version(),'sympy':sympy.__version__,
                    'mpmath':mpmath.__version__}
        demand(tuple(map(int,platform.python_version_tuple())) >= (3,10,0),
               'Python>=3.10 required')
        run_script(out,SOURCES[0],[],['ALL CHECKS PASSED',
            'invalid-input safeguard suite (26 malformed cases rejected)',
            'high-precision effective post-check syndrome-noise'])
        if args.optimized_guard_check:
            cp = subprocess.run([sys.executable,'-O',str(ROOT/SOURCES[0])],cwd=out,
                                capture_output=True,text=True,timeout=240)
            (out/'extended_optimized.log').write_text(cp.stdout+'\n'+cp.stderr)
            demand(cp.returncode == 0 and 'ALL CHECKS PASSED' in cp.stdout
                   and 'invalid-input safeguard suite (26 malformed cases rejected)' in cp.stdout,
                   'optimized -O main verifier failed')
        run_script(out,SOURCES[1],['--json',out/'independent_review.json'],
                   ['ALL INDEPENDENT CHECKS PASSED'])
        run_script(out,SOURCES[2],[],['ALL HETEROGENEOUS-BOOTSTRAP CHECKS PASSED',
                   'NUMERICAL ONLY (not part of the exact theorem)'])
        run_script(out,SOURCES[3],['--output',out/'h12_recheck','--reference',ROOT/SOURCES[2]],
                   ['ALL INDEPENDENT RECHECKS PASSED'])
        run_script(out,SOURCES[4],[],['ALL CONTROL COMPARISONS PASSED'])
        run_script(out,SOURCES[5],[],['ALL EVEN-SECTOR TWIRL CHECKS PASSED'])
        run_script(out,SOURCES[6],['--output',out/'comparison_results.json'],
                   ['ALL COMPARISON CHECKS PASSED'])
        result = verify_results(out)
        hashes = {file:hashlib.sha256((ROOT/file).read_bytes()).hexdigest() for file in SOURCES}
        summary={'status':'PASS', 'environment':versions,
                 'script_sha256':hashes,'scientific_check_results':result,
                 'limits':['Published Bravyi braid-and-fusion instrument is an assumption, not fully re-simulated',
                           'The Gaussian implementation of affine syndrome transformations is analytic',
                           'Convex-Gaussian obstruction imports an existing theorem',
                           'Search algorithm/optimality not reproduced',
                           'Postselection source parity, independent inputs, and ideal storage are assumed',
                           'Numeric effective operation-noise roots are not global certified device thresholds']}
        (out/'verification_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print('ALL REPRODUCIBILITY CHECKS PASSED; see',out/'verification_summary.json',flush=True)
    except Exception as exc:
        (out/'verification_failed.txt').write_text(f'{type(exc).__name__}: {exc}\n')
        print('FAILED',type(exc).__name__,str(exc),file=sys.stderr,flush=True)
        raise SystemExit(1) from exc


if __name__=='__main__':
    main()
