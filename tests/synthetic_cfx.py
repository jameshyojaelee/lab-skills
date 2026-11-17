"""Write a synthetic Bio-Rad CFX Maestro export for a KAPA-style library quant run with known answers."""
import csv, os
import numpy as np

SLOPE, ICPT = -3.40, 11.5       # Cq at 1 pM ~ 11.5
STD_PM = [20, 2, 0.2, 0.02, 0.002, 0.0002]
STD_BP = 452
# name: (undiluted nM size-corrected, size bp)
LIBS = {'L1': (10.0, 400), 'L2': (2.5, 450), 'L3': (40.0, 380)}


def cq(pm):
    return SLOPE * np.log10(pm) + ICPT


def write(d, seed=0):
    rng = np.random.default_rng(seed)
    os.makedirs(d, exist_ok=True)
    rows = []
    for i, p in enumerate(STD_PM):
        rows.append((f"{'ABCDEF'[i]}12", cq(p)))
    for c, (name, (nM, bp)) in enumerate(LIBS.items(), start=1):
        und_pM = nM * 1000 * bp / STD_BP          # uncorrected pM the assay sees
        for rws, fold in (('ABC', 1e4), ('DEF', 1e5)):
            for r in rws:
                rows.append((f'{r}{c:02d}', cq(und_pM / fold) + rng.normal(0, 0.03)))
    rows.append(('G12', float('nan')))
    with open(os.path.join(d, 'run -  Quantification Cq Results.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['', 'Well', 'Fluor', 'Target', 'Content', 'Sample', 'Biological Set Name', 'Cq'])
        for well, v in rows:
            m = well[0] + f'{int(well[1:]):02d}'
            w.writerow(['', m, 'SYBR', '', 'Unkn', '', '', 'NaN' if v != v else f'{v:.4f}'])
    return LIBS


if __name__ == '__main__':
    import sys
    write(sys.argv[1] if len(sys.argv) > 1 else 'synthetic_cfx')
