#!/usr/bin/env python3
"""Countess 3 / 3 FL exports -> per-read table, replicate summary, fluorescence with CIs, seeding math.

    python3 cell_count.py PATH [PATH ...] [--group-regex REGEX] [--seed-cells N] [--outdir DIR]

PATH: a Countess export folder (one per read: <name>.csv, <name>_Cloud.json, images, .fcs),
a parent folder (searched recursively), or a CSV file. Every CSV whose header starts with
'Count ID' is read; a sibling <name>_Cloud.json, when present, adds per-cell data.

--group-regex: groups reads into replicates by sample name; the first capture group is the group
  (default: the sample name with a trailing _rep<N>/_<N>/-<N> removed).
--seed-cells N: also report the µl of suspension that holds N live cells.

Outputs: reads.csv, groups.csv, summary.json, counts.png, fluorescence.png (if any fluorescence),
brightness.png (if Cloud.json present).
"""
import argparse, csv, glob, json, math, os, re, sys
import numpy as np


def wilson(k, n, z=1.96):
    if not n:
        return float('nan'), float('nan')
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float('nan')


def find_csvs(paths):
    out = []
    for p in paths:
        cands = [p] if os.path.isfile(p) else glob.glob(os.path.join(p, '**', '*.csv'), recursive=True)
        for c in cands:
            try:
                with open(c, newline='', encoding='utf-8-sig') as fh:
                    if fh.readline().startswith('Count ID'):
                        out.append(c)
            except UnicodeDecodeError:
                pass
    return sorted(set(out))


def cloud(csv_path):
    j = csv_path[:-4] + '_Cloud.json'
    if not os.path.exists(j):
        return None
    d = json.load(open(j))
    res = {}
    for ct in d.get('CellTypeResults', []):
        cells = [c for c in ct.get('Cells', []) if c.get('IsValid', True)]
        res[ct['CellType']] = dict(cube=(ct.get('Cube') or {}).get('Name'), cells=cells,
                                   light=ct.get('NormalizedLightIntensity'))
    bf = res.get('Brightfield')
    if not bf:
        return None
    cells = bf['cells']
    out = dict(n_objects=len(cells),
               n_clumped=sum(c.get('ClumpIndex') is not None for c in cells),
               size_um=[c.get('Size') for c in cells])
    if cells and cells[0].get('FluorescentBrightness'):
        b = np.array([c['FluorescentBrightness'][0] for c in cells], float)
        f = np.array([bool(c['IsFluorescent'][0]) for c in cells])
        free = np.array([c.get('ClumpIndex') is None for c in cells])
        out.update(brightness=b.tolist(), positive=f.tolist(),
                   k_pos=int(f.sum()), k_pos_nonclumped=int(f[free].sum()), n_nonclumped=int(free.sum()),
                   effective_gate=float(b[f].min()) if f.any() else None,
                   cube=next((v['cube'] for v in res.values() if v['cube']), None))
    return out


def read_one(path, default_group):
    with open(path, newline='', encoding='utf-8-sig') as fh:
        rows = list(csv.DictReader(fh))
    recs = []
    for r in rows:
        name = r.get('Sample name', '').strip() or os.path.basename(path)[:-4]
        rec = dict(file=path, count_id=r.get('Count ID'), sample=name, time=r.get('Date & Time'),
                   mode=r.get('Count mode'),
                   total_per_ml=fnum(r.get('Total concentration')), live_per_ml=fnum(r.get('Live concentration')),
                   dead_per_ml=fnum(r.get('Dead concentration')), total_counted=fnum(r.get('Total cells counted')),
                   live_counted=fnum(r.get('Live cells counted')), dead_counted=fnum(r.get('Dead cells counted')),
                   viability_pct=fnum(r.get('Viability (%)')), live_size_um=fnum(r.get('Live average size (um)')),
                   aggregation_pct=fnum(r.get('Aggregation (%)')),
                   stain_corrected=r.get('1:1 Stain corrected'), predilution_corrected=r.get('Pre-Dilution corrected'))
        if rec['live_counted'] == 0 and rec['dead_counted'] == 0 and rec['total_counted'] > 0:
            # no viability stain (e.g. brightfield-only count): live/dead are not measured, not zero
            rec['viability_pct'] = rec['live_per_ml'] = rec['dead_per_ml'] = float('nan')
            rec['live_size_um'] = float('nan')
            rec['viability_measured'] = False
        else:
            rec['viability_measured'] = True
        for cube in ('Cube 1', 'Cube 2'):
            nm = (r.get(f'{cube} name') or '').strip()
            if nm:
                k = fnum(r.get(f'{cube} cells counted')); n = rec['total_counted']
                lo, hi = wilson(k, n) if n == n and k == k else (float('nan'),) * 2
                rec.update({f'{cube.lower().replace(" ", "")}_name': nm, f'{cube.lower().replace(" ", "")}_pct': fnum(r.get(f'{cube} (%)')),
                            f'{cube.lower().replace(" ", "")}_counted': k,
                            f'{cube.lower().replace(" ", "")}_ci_lo': lo, f'{cube.lower().replace(" ", "")}_ci_hi': hi})
        c = cloud(path)
        if c:
            rec['objects_in_json'] = c['n_objects']; rec['clumped_objects'] = c['n_clumped']
            if 'k_pos' in c:
                rec['json_pos'] = c['k_pos']
                rec['pos_pct_nonclumped'] = 100 * c['k_pos_nonclumped'] / c['n_nonclumped'] if c['n_nonclumped'] else float('nan')
                rec['effective_gate'] = c['effective_gate']
            rec['_cloud'] = c
        rec['group'] = default_group(name)
        recs.append(rec)
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paths', nargs='+'); ap.add_argument('--group-regex'); ap.add_argument('--seed-cells', type=float)
    ap.add_argument('--outdir', default='cell_count_out'); ap.add_argument('--no-plots', action='store_true')
    a = ap.parse_args()
    run(a.paths, a.group_regex, a.seed_cells, a.outdir, not a.no_plots, verbose=True)


def run(paths, group_regex=None, seed_cells=None, outdir='cell_count_out', plots=True, verbose=False):
    if group_regex:
        rx = re.compile(group_regex)
        grp = lambda s: (rx.search(s).group(1) if rx.search(s) else s)
    else:
        grp = lambda s: re.sub(r'([_-](rep)?\d+)$', '', s, flags=re.I)
    recs = []
    for c in find_csvs(paths):
        recs += read_one(c, grp)
    if not recs:
        raise SystemExit('no Countess CSV (header starting with "Count ID") found')
    recs.sort(key=lambda r: (r['time'] or '', str(r['count_id'])))
    os.makedirs(outdir, exist_ok=True)
    warn = []
    for r in recs:
        if not r['viability_measured']:
            warn.append(f"{r['sample']}: no viability stain in this read; live/dead and viability not measured"
                        + ("; seeding uses total cells" if seed_cells else ''))
        if seed_cells:
            base = r['live_per_ml'] if r['viability_measured'] else r['total_per_ml']
            r['ul_for_seed'] = seed_cells / base * 1000 if base > 0 else float('nan')
        if r['total_counted'] == r['total_counted'] and r['total_counted'] < 100:
            warn.append(f"{r['sample']}: only {r['total_counted']:.0f} cells counted (counting error ~{100/math.sqrt(max(r['total_counted'],1)):.0f}%)")
        if r['aggregation_pct'] == r['aggregation_pct'] and r['aggregation_pct'] > 20:
            warn.append(f"{r['sample']}: aggregation {r['aggregation_pct']:.0f}%; concentration is underestimated and fluorescence % mixes clumps")
        if 'json_pos' in r and 'cube1_counted' in r and r['json_pos'] != r['cube1_counted']:
            warn.append(f"{r['sample']}: Cloud.json positives {r['json_pos']} != CSV {r['cube1_counted']:.0f}")
        if r['total_per_ml'] == r['total_per_ml'] and not (1e5 <= r['total_per_ml'] <= 1e7):
            warn.append(f"{r['sample']}: {r['total_per_ml']:.3g} cells/ml is outside the ~1e5-1e7 accurate range; dilute or concentrate")
    keys = [k for k in dict.fromkeys(k for r in recs for k in r) if not k.startswith('_')]
    with open(os.path.join(outdir, 'reads.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction='ignore'); w.writeheader(); w.writerows(recs)
    groups = []
    for g in dict.fromkeys(r['group'] for r in recs):
        rr = [r for r in recs if r['group'] == g]
        def ms(k):
            v = np.array([r.get(k, np.nan) for r in rr], float)
            v = v[np.isfinite(v)]
            return (float(v.mean()) if len(v) else float('nan'), float(v.std(ddof=1) / v.mean() * 100) if len(v) > 1 and v.mean() else float('nan'))
        row = dict(group=g, n_reads=len(rr))
        for k in ('total_per_ml', 'live_per_ml', 'viability_pct', 'live_size_um', 'cube1_pct', 'aggregation_pct'):
            m, cv = ms(k); row[k] = m; row[k + '_cv_pct'] = cv
        if 'cube1_counted' in rr[0]:
            kk = sum(r.get('cube1_counted', 0) for r in rr); nn = sum(r['total_counted'] for r in rr)
            row['cube1_pooled_pct'] = 100 * kk / nn if nn else float('nan')
            row['cube1_pooled_ci_lo'], row['cube1_pooled_ci_hi'] = wilson(kk, nn)
        if seed_cells:
            base = row['live_per_ml'] if row['live_per_ml'] == row['live_per_ml'] else row['total_per_ml']
            row['ul_for_seed'] = seed_cells / base * 1000 if base > 0 else float('nan')
        groups.append(row)
    with open(os.path.join(outdir, 'groups.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(dict.fromkeys(k for g in groups for k in g))); w.writeheader(); w.writerows(groups)
    json.dump(dict(n_reads=len(recs), groups=groups, warnings=warn, seed_cells=seed_cells),
              open(os.path.join(outdir, 'summary.json'), 'w'), indent=1, default=float)
    if plots:
        plot(recs, outdir)
    if verbose:
        for g in groups:
            fl = f"  {g.get('cube1_pooled_pct', float('nan')):.1f}% fluor (95% CI {g.get('cube1_pooled_ci_lo', float('nan')):.1f}-{g.get('cube1_pooled_ci_hi', float('nan')):.1f})" if 'cube1_pooled_pct' in g else ''
            sd = f"  {g['ul_for_seed']:.1f} µl for {seed_cells:g} cells (live if stained)" if seed_cells else ''
            print(f"{g['group']:32s} n={g['n_reads']}  total {g['total_per_ml']:.3g}/ml  live {g['live_per_ml']:.3g}/ml  "
                  f"viab {g['viability_pct']:.1f}%{fl}{sd}")
        for x in warn:
            print('WARNING', x)
    return recs, groups, warn


sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', '..', '_shared'))
import labstyle


def plot(recs, outdir):
    """Rules: _shared/FIGURES.md."""
    x = np.arange(len(recs)); lab = [r['sample'] for r in recs]
    kind = 'single' if len(recs) <= 5 else 'double'
    stained = [r for r in recs if r.get('viability_measured', True)]
    if stained:
        fig, ax = labstyle.figure('double', ncols=2)
        ax[0].bar(x - 0.2, [r['live_per_ml'] for r in recs], 0.4, color=labstyle.STAGES[0], label='Live')
        ax[0].bar(x + 0.2, [r['dead_per_ml'] for r in recs], 0.4, color=labstyle.STAGES[2], label='Dead')
        ax[0].set_ylabel('Cells per ml'); ax[0].legend(loc='upper right'); ax[0].set_title('Concentration')
        ax[1].bar(x, [r['viability_pct'] for r in recs], 0.6, color=labstyle.BLUE)
        ax[1].set_ylabel('Viability (%)'); ax[1].set_ylim(0, 100); ax[1].set_title('Viability')
        for a in ax:
            a.set_xticks(x); a.set_xticklabels(lab, rotation=45, ha='right')
        cap = 'Countess counts per read (cells per ml as exported). Reads without a viability stain show no bars.'
    else:
        fig, ax = labstyle.figure(kind, ncols=1)
        ax.bar(x, [r['total_per_ml'] for r in recs], 0.6, color=labstyle.BLUE)
        ax.set_ylabel('Total objects per ml'); ax.set_xticks(x); ax.set_xticklabels(lab, rotation=45, ha='right')
        cap = 'Countess total objects per ml per read. No viability stain: live/dead not measured.'
    labstyle.save(fig, os.path.join(outdir, 'counts'), caption=cap)
    if any('cube1_pct' in r for r in recs):
        fig, a = labstyle.figure(kind, ncols=1)
        v = np.array([r.get('cube1_pct', np.nan) for r in recs])
        lo = np.array([r.get('cube1_ci_lo', np.nan) for r in recs]); hi = np.array([r.get('cube1_ci_hi', np.nan) for r in recs])
        a.errorbar(x, v, yerr=[v - lo, hi - v], fmt='o', color=labstyle.BLUE, ms=4, label='All objects')
        if any('pos_pct_nonclumped' in r for r in recs):
            a.plot(x + 0.15, [r.get('pos_pct_nonclumped', np.nan) for r in recs], 'o', mfc='white',
                   mec=labstyle.BLUE, ms=4, label='Non-clumped')
            a.legend(loc='best')
        a.set_ylim(0, 100)
        a.set_ylabel(f"{recs[0].get('cube1_name', 'Fluor')}+ objects (%)")
        a.set_xticks(x); a.set_xticklabels(lab, rotation=45, ha='right')
        labstyle.save(fig, os.path.join(outdir, 'fluorescence'), caption=(
            'Fluorescent-positive objects / counted objects per read; error bars are Wilson 95% intervals on '
            'counted objects. Open: positive among non-clumped objects (from Cloud.json).'))
    cl = [r for r in recs if r.get('_cloud', {}).get('brightness')]
    if cl:
        fig, ax = labstyle.figure('single', nrows=len(cl), size=(3.5, min(0.5 * len(cl) + 0.7, 9)), sharex=True)
        ax = np.atleast_1d(ax)
        bins = np.linspace(0, 255, 64)
        for a, r in zip(ax, cl):
            b = np.array(r['_cloud']['brightness']); f = np.array(r['_cloud']['positive'])
            a.hist(b[~f], bins=bins, color=labstyle.STAGES[2]); a.hist(b[f], bins=bins, color=labstyle.STAGES[0])
            if r.get('effective_gate'):
                a.axvline(r['effective_gate'], color=labstyle.INK, lw=0.6, ls=(0, (3, 2)))
            a.set_title(r['sample'], fontsize=6, pad=1)
            a.set_yticks([]); a.spines['left'].set_visible(False)
            if a is not ax[-1]:
                a.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        ax[-1].set_xlabel('Fluorescent brightness per object')
        labstyle.save(fig, os.path.join(outdir, 'brightness'), caption=(
            'Per-object fluorescent brightness from Cloud.json. Dark: instrument-positive; light: negative; dashed '
            'line: lowest brightness called positive (effective gate).'))


if __name__ == '__main__':
    main()
