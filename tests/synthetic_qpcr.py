"""Write a synthetic QuantStudio-style export with known truth, for testing qpcr_quant.py."""
import numpy as np, pandas as pd

SLOPE, ICPT = -3.373, 33.81
UNIT, BG = 1.4e5, 1.2e4          # copies per bead, free-oligo background per well
ROWS = 'ABCDEFGHIJKLMNOP'


def ct(copies):
    return SLOPE * np.log10(copies) + ICPT


def make(path, seed=0):
    rng = np.random.default_rng(seed)
    rec = []
    top = 10e-6 / 1e6 * 2e-6 * 6.02214076e23
    for i, rr in enumerate(['AB', 'CD', 'EF', 'GH', 'IJ', 'KL']):
        for r in rr:
            rec.append((f'{r}9', ct(top / 10 ** i) + rng.normal(0, 0.05), 85.5))
    lam = 0.25
    k = rng.poisson(lam, 96)
    for j, (c, r) in enumerate([(c, r) for c in range(1, 7) for r in ROWS]):
        n = BG + k[j] * UNIT * rng.lognormal(0, 0.15)
        rec.append((f'{r}{c}', ct(n) + rng.normal(0, 0.05), 85.5))
    for r in ROWS:
        rec.append((f'{r}7', ct(5 * BG + rng.poisson(1.0) * UNIT + 1) , 85.5))
    for r in 'ACEGIKMO':
        rec.append((f'{r}8', np.nan if r in 'AC' else 37.0, 79.0))
    rec += [('M9', 38.0, 78.0), ('N9', np.nan, np.nan)]
    true_bulk = {'S1': 1.9e5, 'S2': 1.6e5, 'S3': 2.8e5}
    counts = {'S1': 955, 'S2': 1295, 'S3': 995}
    for smp, rows in [('S1', 'ABCD'), ('S2', 'EFGH'), ('S3', 'IJKL')]:
        beads = counts[smp] * 100 * 5
        total = true_bulk[smp] * beads
        for r, fold in zip(rows, [1e3, 1e3, 1e4, 1e4]):
            rec.append((f'{r}10', ct(total / fold / 200 * 2), 85.5))
    res = pd.DataFrame(rec, columns=['Well Position', 'CT', 'Tm'])
    res.insert(0, 'Well', range(1, len(res) + 1))
    pre = pd.DataFrame([['Block Type', '384-Well Block'], ['Experiment Name', 'synthetic'], [None, None]])
    with pd.ExcelWriter(path) as xw:
        for name, df in [('Results', res[['Well', 'Well Position', 'CT']]), ('Melt Curve Result', res[['Well', 'Well Position', 'Tm']])]:
            pre.to_excel(xw, sheet_name=name, header=False, index=False)
            df.to_excel(xw, sheet_name=name, startrow=len(pre), index=False)
    return dict(slope=SLOPE, intercept=ICPT, unit=UNIT, bg=BG, lam=lam, bulk=true_bulk, k=k)
