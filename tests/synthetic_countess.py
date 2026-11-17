"""Write synthetic Countess 3 FL exports (CSV + Cloud.json per read) with known answers."""
import csv, json, os
import numpy as np

HEADER = ['Count ID', 'Session ID', 'Sample name', 'Date & Time', 'Count mode', 'Type', '1:1 Stain corrected',
          'Pre-Dilution corrected', 'Total concentration', 'Total cells counted', 'Live concentration',
          'Live cells counted', 'Dead concentration', 'Dead cells counted', 'Viability (%)',
          'Live average size (um)', 'Dead average size (um)', 'Cube 1 name', 'Cube 1 concentration', 'Cube 1 (%)',
          'Cube 1 cells counted', 'Aggregation (%)']
# sample: (total counted, live fraction, fluorescent count, concentration per counted cell)
READS = [('cellsA_rep1', 400, 0.90, 120), ('cellsA_rep2', 380, 0.92, 110), ('cellsB_rep1', 500, 0.80, 25)]
PER_CELL = 2500.0   # cells/ml per counted cell


def write(root, seed=0):
    rng = np.random.default_rng(seed)
    for i, (name, n, live, k) in enumerate(READS):
        d = os.path.join(root, f'2026-01-01T10-0{i}-00_{name}')
        os.makedirs(d, exist_ok=True)
        nl = int(round(n * live)); nd = n - nl
        row = [1000 + i, 1, name, f'2026-01-01T10:0{i}:00', 'Trypan blue', 'FL', 'True', 'False',
               n * PER_CELL, n, nl * PER_CELL, nl, nd * PER_CELL, nd, round(100 * nl / n, 2), 14.5, 12.0,
               'GFP', k * PER_CELL, round(100 * k / n, 2), k, 5.0]
        with open(os.path.join(d, f'{name}.csv'), 'w', newline='') as fh:
            w = csv.writer(fh); w.writerow(HEADER); w.writerow(row)
        pos = np.zeros(n, bool); pos[rng.choice(n, k, replace=False)] = True
        bright = np.where(pos, rng.normal(150, 20, n), rng.normal(60, 10, n)).clip(0, 255)
        cells = [dict(Size=14.5, IsFluorescent=[bool(p)], FluorescentBrightness=[float(b)], IsValid=True,
                      ClumpIndex=None) for p, b in zip(pos, bright)]
        js = dict(CellTypeResults=[dict(CellType='Brightfield', Cube=None, Cells=cells, NormalizedLightIntensity=72.0),
                                   dict(CellType='Epi1', Cube=dict(Name='GFP'), Cells=[], NormalizedLightIntensity=90.0)])
        json.dump(js, open(os.path.join(d, f'{name}_Cloud.json'), 'w'))
    return READS


if __name__ == '__main__':
    import sys
    write(sys.argv[1] if len(sys.argv) > 1 else 'synthetic_countess')
