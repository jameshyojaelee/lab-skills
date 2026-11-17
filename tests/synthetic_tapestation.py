"""Write a synthetic TapeStation D1000 HS per-well export (same file names, headers, encoding)."""
import csv, os
import numpy as np

LADDER = [25, 50, 100, 200, 300, 400, 500, 700, 1000, 1500]
N = 760
STEM = '2026-01-01 - 10-00-00-HSD1000'
# index of each ladder size in the trace (log-linear plus slight curvature)
POS = np.array([107, 217, 307, 386, 438, 478, 513, 572, 618, 651], float)

# well: (library peak bp, library % area, dimer bp or None, dimer % area, region % of total, molarity pM)
SAMPLES = {
    'A1': (400, 95.0, None, 0.0, 88.0, 14000.0),
    'B1': (420, 70.0, 150, 25.0, 62.0, 2100.0),
    'C1': (380, 60.0, None, 0.0, 41.0, 40.0),
}


def idx_of(bp):
    return float(np.interp(np.log10(bp), np.log10(LADDER), POS))


def trace(peaks):
    x = np.arange(N, dtype=float)
    y = np.zeros(N)
    for bp, h in peaks:
        y += h * np.exp(-0.5 * ((x - idx_of(bp)) / 4.0) ** 2)
    return y


def write(run_dir):
    os.makedirs(run_dir, exist_ok=True)
    fn = STEM.replace('-00-00-', '.00.00.').replace('-HSD', '.HSD')
    P = lambda s: os.path.join(run_dir, f'{STEM}_{s}')
    with open(P('compactPeakTable.csv'), 'w', newline='', encoding='cp1252') as fh:
        w = csv.writer(fh)
        w.writerow(['FileName', 'Well', 'Sample Description', 'Size [bp]', 'Calibrated Conc. [pg/µl]',
                    'Assigned Conc. [pg/µl]', 'Peak Molarity [pmol/l]', '% Integrated Area', 'Peak Comment', 'Observations'])
        for i, s in enumerate(LADDER):
            obs = 'Lower Marker' if i == 0 else 'Upper Marker' if s == 1500 else ''
            w.writerow([fn, 'EL1', 'Electronic Ladder', s, 300, 250 if s == 1500 else '', 1000, '' if obs else 12.5, '', obs])
        for well, (bp, pa, dbp, dpa, _, mol) in SAMPLES.items():
            w.writerow([fn, well, '', 25, 350, '', 21000, '', '', 'Lower Marker'])
            if dbp:
                w.writerow([fn, well, '', dbp, 120, '', 1200, dpa, '', ''])
            w.writerow([fn, well, '', bp, 3000, '', mol, pa, '', ''])
            w.writerow([fn, well, '', 1500, 250, 250, 256, '', '', 'Upper Marker'])
    with open(P('compactRegionTable.csv'), 'w', newline='', encoding='cp1252') as fh:
        w = csv.writer(fh)
        w.writerow(['FileName', 'WellId', 'Sample Description', 'From [bp]', 'To [bp]', 'Average Size [bp]',
                    'Conc. [pg/µl]', 'Region Molarity [pmol/l]', '% of Total', 'Region Comment'])
        for well, (bp, pa, dbp, dpa, rt, mol) in SAMPLES.items():
            w.writerow([fn, well, '', 300, 600, bp + 10, 3000, mol, rt, ''])
    with open(P('runProperties.csv'), 'w', newline='', encoding='cp1252') as fh:
        w = csv.writer(fh)
        w.writerow(['Analysis Software Version', 'Filename', 'Assay', 'Run End Date'])
        w.writerow(['5.1', fn, 'D1000 High  Sensitivity', '1 Jan 2026, 10:10'])
    cols = {'EL1: Electronic Ladder': trace([(s, 300) for s in LADDER])}
    for well, (bp, pa, dbp, dpa, _, _) in SAMPLES.items():
        pk = [(25, 120), (1500, 40), (bp, 6 * pa)]
        if dbp:
            pk.append((dbp, 6 * dpa))
        cols[f'{well}: '] = trace(pk)
    with open(P('Electropherogram.csv'), 'w', newline='', encoding='cp1252') as fh:
        w = csv.writer(fh)
        w.writerow(list(cols))
        for i in range(N):
            w.writerow([f'{c[i]:.4g}' for c in cols.values()])
    with open(os.path.join(run_dir, 'sample_sheet.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['Well', 'Sample_Name', 'Group', 'Sample_Type'])
        for i, well in enumerate(SAMPLES):
            w.writerow([well, f'Library {i + 1}', 'Libraries', 'sample'])
    return SAMPLES


if __name__ == '__main__':
    import sys
    write(sys.argv[1] if len(sys.argv) > 1 else 'synthetic_run')
