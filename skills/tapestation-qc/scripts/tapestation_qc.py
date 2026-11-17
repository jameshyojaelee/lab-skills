#!/usr/bin/env python3
"""Agilent TapeStation run -> per-sample QC table, warnings, electropherograms and a labelled gel.

    python3 tapestation_qc.py RUN_DIR [--sheet sample_sheet.csv] [--region 300 600]
                              [--below 200] [--outdir RUN_DIR/tapestation_qc] [--draft-sheet]

RUN_DIR holds the per-well CSV export from TapeStation Analysis software:
  *_compactPeakTable.csv     required
  *_compactRegionTable.csv   optional (region metrics)
  *_Electropherogram.csv     optional (plots, pseudo-gel)
  *_runProperties.csv        optional (assay, date)
Agilent writes these in cp1252 (the pg/µl header has a micro sign).

Sample sheet (optional, recommended): Well,Sample_Name[,anything else]. The export carries no
sample names, so the sheet is the only record of which lane is which. --draft-sheet writes one.

Outputs: qc_table.csv, qc.json, electropherograms.png, gel_labeled.png.
"""
import argparse, csv, json, os, re, sys
import numpy as np

LADDER_WORDS = ('ladder',)


# ------------------------------------------------------------------ reading
def find(run_dir, suffix, required=True):
    hits = sorted(f for f in os.listdir(run_dir) if f.endswith(suffix))
    if not hits:
        if required:
            sys.exit(f'no *{suffix} in {run_dir}')
        return None
    if len(hits) > 1:
        sys.exit(f'more than one *{suffix} in {run_dir}: {hits}')
    return os.path.join(run_dir, hits[0])


def rows(path):
    with open(path, newline='', encoding='cp1252') as fh:
        return list(csv.DictReader(fh))


def col(r, *names):
    """Header spelling varies between assays and software versions; match by prefix."""
    for n in names:
        for k in r:
            if k and k.strip().lower().startswith(n.lower()):
                return r[k]
    return None


def num(v):
    if v is None:
        return None
    v = str(v).strip().replace(',', '')
    try:
        return float(v) if v not in ('', '-', 'NaN') else None
    except ValueError:
        return None


def well_key(w):
    m = re.match(r'([A-H])(\d+)', w.strip().upper())
    return (int(m.group(2)), m.group(1)) if m else (99, w)


def is_ladder(well, desc):
    return well.upper().startswith('EL') or any(x in (desc or '').lower() for x in LADDER_WORDS)


def read_peaks(run_dir):
    out = {}
    for r in rows(find(run_dir, '_compactPeakTable.csv')):
        w = (col(r, 'Well') or '').strip()
        out.setdefault(w, dict(desc=(col(r, 'Sample Description') or '').strip(), peaks=[]))
        out[w]['peaks'].append(dict(
            size_bp=num(col(r, 'Size')), conc_pg_ul=num(col(r, 'Calibrated Conc')),
            molarity_pM=num(col(r, 'Peak Molarity')), pct_area=num(col(r, '% Integrated Area')),
            comment=(col(r, 'Peak Comment') or '').strip(), observation=(col(r, 'Observations') or '').strip()))
    return out


def read_regions(run_dir):
    p = find(run_dir, '_compactRegionTable.csv', required=False)
    out = {}
    if not p:
        return out
    for r in rows(p):
        w = (col(r, 'WellId', 'Well') or '').strip()
        out.setdefault(w, []).append(dict(
            frm=num(col(r, 'From')), to=num(col(r, 'To')), avg_bp=num(col(r, 'Average Size')),
            conc_pg_ul=num(col(r, 'Conc')), molarity_pM=num(col(r, 'Region Molarity')),
            pct_total=num(col(r, '% of Total')), comment=(col(r, 'Region Comment') or '').strip()))
    return out


def read_props(run_dir):
    p = find(run_dir, '_runProperties.csv', required=False)
    if not p:
        return {}
    r = rows(p)
    return {k.strip(): (v or '').strip() for k, v in r[0].items() if k} if r else {}


def read_trace(run_dir):
    p = find(run_dir, '_Electropherogram.csv', required=False)
    if not p:
        return None
    with open(p, newline='', encoding='cp1252') as fh:
        rd = csv.reader(fh)
        hdr = next(rd)
        data = np.array([[float(x or 0) for x in r] for r in rd if r])
    wells = [h.split(':')[0].strip() for h in hdr]
    return {w: data[:, i] for i, w in enumerate(wells)}


def read_sheet(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline='', encoding='utf-8-sig') as fh:
        return {r['Well'].strip(): r for r in csv.DictReader(fh)}


# ------------------------------------------------------------------ sizing
def peak_positions(y, n, min_dist=5):
    from scipy.signal import find_peaks
    pk, pr = find_peaks(y, distance=min_dist, prominence=max(y.max() * 0.02, 1e-9))
    if len(pk) < n:
        return None
    top = pk[np.argsort(pr['prominences'])[::-1][:n]]
    return np.sort(top)


def size_axis(trace, peaks):
    """Map data-point index -> bp for every lane.
    Ladder lane: its peaks are matched in order to the ladder sizes in the peak table.
    Sample lanes: lower and upper markers are aligned to the ladder's (linear in index),
    then the ladder's log-size interpolation is applied."""
    lad = [w for w in trace if is_ladder(w, peaks.get(w, {}).get('desc', ''))]
    if not lad or lad[0] not in peaks:
        return None
    L = lad[0]
    sizes = sorted(p['size_bp'] for p in peaks[L]['peaks'] if p['size_bp'])
    pos = peak_positions(trace[L], len(sizes))
    if pos is None:
        return None
    logs = np.log10(sizes)
    lo_m, hi_m = pos[0], pos[-1]
    axes = {}
    for w, y in trace.items():
        if w == L:
            idx = np.arange(len(y), dtype=float)
        else:
            mk = [p for p in peaks.get(w, {}).get('peaks', []) if 'marker' in p['observation'].lower()]
            from scipy.signal import find_peaks
            pk, pr = find_peaks(y, prominence=max(y.max() * 0.005, 1e-9))
            idx = np.arange(len(y), dtype=float)
            if len(mk) >= 2 and len(pk) >= 2:
                # markers: most prominent peak within +-12% of the trace length of the ladder's markers
                win = 0.12 * len(y)
                def near(c):
                    m = np.abs(pk - c) <= win
                    return pk[m][np.argmax(pr['prominences'][m])] if m.any() else None
                a, b = near(lo_m), near(hi_m)
                if a is not None and b is not None and b > a:
                    idx = lo_m + (idx - a) * (hi_m - lo_m) / (b - a)
        axes[w] = 10 ** np.interp(idx, pos, logs, left=np.nan, right=np.nan)
    return axes


# ------------------------------------------------------------------ QC
def lane_qc(well, info, regs, region, below):
    pk = [p for p in info['peaks'] if 'marker' not in p['observation'].lower()]
    inside = [p for p in pk if 'outside' not in p['comment'].lower()]
    main = max(inside or pk, key=lambda p: p['pct_area'] or 0) if pk else None
    reg = None
    for r in regs or []:
        if r['frm'] == region[0] and r['to'] == region[1]:
            reg = r
    below_pct = sum(p['pct_area'] or 0 for p in pk if (p['size_bp'] or 0) < below)
    rec = dict(well=well,
               main_peak_bp=main and main['size_bp'], main_peak_pct_area=main and main['pct_area'],
               main_peak_conc_pg_ul=main and main['conc_pg_ul'],
               other_peaks='; '.join(f"{p['size_bp']:g} bp ({p['pct_area']:g}%)" for p in pk
                                     if p is not main and p['size_bp'] and p['pct_area']),
               region=f'{region[0]:g}-{region[1]:g}',
               region_avg_bp=reg and reg['avg_bp'], region_conc_pg_ul=reg and reg['conc_pg_ul'],
               region_molarity_nM=reg and reg['molarity_pM'] and reg['molarity_pM'] / 1000,
               region_pct_total=reg and reg['pct_total'], pct_area_below=round(below_pct, 2))
    warn = []
    if not reg:
        have = ', '.join(f"{r['frm']:g}-{r['to']:g}" for r in regs or []) or 'none'
        warn.append(f'{well}: no {region[0]:g}-{region[1]:g} bp region in the export (regions: {have})')
    elif reg['pct_total'] is not None and reg['pct_total'] < 50:
        warn.append(f"{well}: region holds only {reg['pct_total']:.1f}% of the lane")
    if below_pct >= 5:
        warn.append(f'{well}: {below_pct:.1f}% of area below {below:g} bp (primer/adapter dimer?)')
    if any('outside' in p['comment'].lower() for p in pk):
        warn.append(f'{well}: peak(s) outside the sizing range')
    if not pk:
        warn.append(f'{well}: no sample peaks')
    return rec, warn


# ------------------------------------------------------------------ figures (rules: _shared/FIGURES.md)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', '..', '_shared'))
import labstyle


def plot_traces(trace, axes, wells, labels, stem):
    n = len(wells)
    fig, ax = labstyle.figure('single', nrows=n, size=(3.5, min(0.62 * n + 0.6, 9)), sharex=True)
    ax = np.atleast_1d(ax)
    for a, w in zip(ax, wells):
        x = axes[w] if axes else np.arange(len(trace[w]))
        a.plot(x, trace[w], color=labstyle.BLUE, lw=0.8)
        if axes:
            inner = np.isfinite(x) & (x > 35) & (x < 1100)
            if inner.any() and trace[w][inner].max() > 0:
                a.set_ylim(0, trace[w][inner].max() * 1.25)   # markers can run off the top
        a.set_title(labels.get(w, w), fontsize=6, pad=1)
        a.set_yticks([]); a.spines['left'].set_visible(False)
        if a is not ax[-1]:
            a.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    if axes:
        ax[-1].set_xscale('log'); ax[-1].set_xlim(20, 2000)
        ax[-1].set_xticks([25, 50, 100, 200, 500, 1000]); ax[-1].set_xticklabels([25, 50, 100, 200, 500, 1000])
        ax[-1].set_xlabel('Size (bp)')
    else:
        ax[-1].set_xlabel('Data point (no ladder found)')
    labstyle.save(fig, stem, caption='Electropherogram per lane, sized by the electronic ladder after aligning each '
                  'lane to its lower and upper markers; y scaled per lane to the largest non-marker signal.')


def plot_gel(trace, axes, wells, labels, groups, stem):
    """Pseudo-gel rendered from the traces on a common log-bp axis, one intensity scale for all lanes."""
    if not axes:
        return
    grid = np.logspace(np.log10(25), np.log10(1500), 400)
    cols = {}
    for w in wells:
        x, y = axes[w], trace[w]
        ok = np.isfinite(x)
        if ok.sum() >= 3:
            cols[w] = np.interp(grid, x[ok], y[ok])
    inner = (grid > 35) & (grid < 1100)
    vmax = max(c[inner].max() for c in cols.values()) or 1
    fig, ax = labstyle.figure('image', size=(min(0.32 * len(wells) + 1.0, 7.0), 3.0))
    for i, w in enumerate(wells):
        if w in cols:
            ax.imshow(np.clip(cols[w] / vmax, 0, 1)[:, None] ** 0.5, cmap='gray_r', vmin=0, vmax=1, aspect='auto',
                      interpolation='nearest', extent=[i - 0.38, i + 0.38, 0, 1], origin='lower')
    ticks = [25, 50, 100, 200, 300, 500, 700, 1000, 1500]
    ax.set_xlim(-0.6, len(wells) - 0.4); ax.set_ylim(0, 1)
    ax.set_yticks([np.log10(t / 25) / np.log10(1500 / 25) for t in ticks]); ax.set_yticklabels(ticks)
    ax.set_ylabel('Size (bp)')
    ax.set_xticks(range(len(wells)))
    ax.set_xticklabels([labels.get(w, w) for w in wells], rotation=60, ha='right')
    ax.tick_params(axis='x', length=0)
    ax.spines['bottom'].set_visible(False)
    for g in groups:
        idx = [wells.index(w) for w in g['wells'] if w in wells]
        if idx:
            ax.plot([min(idx) - 0.4, max(idx) + 0.4], [1.03, 1.03], color=labstyle.INK, lw=0.6, clip_on=False,
                    transform=ax.get_xaxis_transform())
            ax.text(np.mean(idx), 1.05, g['label'], ha='center', va='bottom', fontsize=6.5, transform=ax.get_xaxis_transform())
    labstyle.save(fig, stem, caption='Pseudo-gel rendered from the TapeStation traces on the ladder size scale; one '
                  'intensity scale for all lanes (set by the brightest sample band, square-root contrast).')


# ------------------------------------------------------------------ main
def run(run_dir, sheet=None, region=(300, 600), below=200, outdir=None, plots=True):
    peaks = read_peaks(run_dir)
    regs = read_regions(run_dir)
    props = read_props(run_dir)
    sh = read_sheet(sheet or os.path.join(run_dir, 'sample_sheet.csv'))
    outdir = outdir or os.path.join(run_dir, 'tapestation_qc')
    os.makedirs(outdir, exist_ok=True)
    wells = sorted((w for w, i in peaks.items() if not is_ladder(w, i['desc'])), key=well_key)
    recs, warns = [], []
    for w in wells:
        rec, wn = lane_qc(w, peaks[w], regs.get(w), region, below)
        s = sh.get(w, {})
        rec = dict(sample=s.get('Sample_Name', '') or peaks[w]['desc'], **rec)
        for k, v in s.items():
            if k not in ('Well', 'Sample_Name') and k:
                rec[k] = v
        recs.append(rec); warns += wn
    with open(os.path.join(outdir, 'qc_table.csv'), 'w', newline='') as fh:
        keys = list(dict.fromkeys(k for r in recs for k in r))
        w_ = csv.DictWriter(fh, fieldnames=keys); w_.writeheader(); w_.writerows(recs)
    out = dict(run_dir=run_dir, assay=props.get('Assay'), run_end=props.get('Run End Date'),
               software=props.get('Analysis Software Version'), region=list(region), below_bp=below,
               samples=recs, warnings=warns, sample_sheet=bool(sh))
    if not sh:
        out['warnings'].insert(0, 'no sample sheet: lanes are identified by well only')
    trace = read_trace(run_dir)
    axes = size_axis(trace, peaks) if trace else None
    if trace and axes:
        # size check: main peak from the trace vs the instrument's peak table
        from scipy.signal import find_peaks
        chk = {}
        for r in recs:
            w, want = r['well'], r['main_peak_bp']
            if w in trace and want and (r['main_peak_pct_area'] or 0) >= 20:
                x, y = axes[w], trace[w]
                ok = np.isfinite(x) & (x > 35) & (x < 1100)
                if ok.any():
                    chk[w] = dict(table_bp=want, trace_bp=float(x[ok][np.argmax(y[ok])]))
        out['size_check'] = chk
    json.dump(out, open(os.path.join(outdir, 'qc.json'), 'w'), indent=1, default=str)
    if plots and trace:
        labels = {w: (sh.get(w, {}).get('Sample_Name') or w) for w in trace}
        lad = [w for w in trace if is_ladder(w, peaks.get(w, {}).get('desc', ''))]
        order = lad + wells
        plot_traces(trace, axes, wells, labels, os.path.join(outdir, 'electropherograms'))
        groups = []
        if sh and any('Group' in s for s in sh.values()):
            for g in dict.fromkeys(s.get('Group', '') for s in sh.values()):
                if g:
                    groups.append(dict(label=g, wells=[w for w, s in sh.items() if s.get('Group') == g]))
        plot_gel(trace, axes, order, {**labels, **{l: 'Ladder' for l in lad}}, groups, os.path.join(outdir, 'gel_labeled'))
    return out


def draft_sheet(run_dir):
    p = os.path.join(run_dir, 'sample_sheet.csv')
    if os.path.exists(p):
        sys.exit(f'{p} exists; edit it or remove it first')
    peaks = read_peaks(run_dir)
    wells = sorted((w for w, i in peaks.items() if not is_ladder(w, i['desc'])), key=well_key)
    with open(p, 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['Well', 'Sample_Name', 'Group', 'Sample_Type'])
        for x in wells:
            w.writerow([x, peaks[x]['desc'], '', 'sample'])
    print(f'wrote {p} with {len(wells)} lanes; fill Sample_Name (and Group for gel brackets)')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir'); ap.add_argument('--sheet'); ap.add_argument('--outdir')
    ap.add_argument('--region', nargs=2, type=float, default=[300, 600])
    ap.add_argument('--below', type=float, default=200)
    ap.add_argument('--draft-sheet', action='store_true'); ap.add_argument('--no-plots', action='store_true')
    a = ap.parse_args()
    if a.draft_sheet:
        return draft_sheet(a.run_dir)
    o = run(a.run_dir, a.sheet, tuple(a.region), a.below, a.outdir, not a.no_plots)
    hdr = f"{'well':5s} {'sample':24s} {'peak bp':>8s} {'region bp':>9s} {'nM':>8s} {'%total':>7s} {'<'+str(int(a.below))+'bp%':>8s}"
    print(f"{o['assay']}  {o['run_end']}\n{hdr}")
    for r in o['samples']:
        f = lambda v, fmt: (fmt % v) if v is not None else '-'
        print(f"{r['well']:5s} {str(r['sample'])[:24]:24s} {f(r['main_peak_bp'],'%8.0f')} {f(r['region_avg_bp'],'%9.0f')} "
              f"{f(r['region_molarity_nM'],'%8.3f')} {f(r['region_pct_total'],'%7.1f')} {r['pct_area_below']:8.1f}")
    for w in o['warnings']:
        print('WARNING', w)


if __name__ == '__main__':
    main()
