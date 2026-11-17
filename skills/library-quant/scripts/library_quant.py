#!/usr/bin/env python3
"""qPCR library quantification against DNA standards of known concentration (e.g. KAPA Library
Quantification, Illumina P5/P7 primers), with size correction and a dilution-agreement check.

    python3 library_quant.py config.yaml [--export-dir DIR] [--outdir DIR]

Reads Bio-Rad CFX Maestro CSV exports ("Quantification Cq Results.csv", optional
"Quantification Amplification Results_SYBR.csv", "Melt Curve Peak Results.csv") or a QuantStudio
Results sheet (.xls/.xlsx). Outputs: libraries.csv, standards.csv, summary.json, library_quant.png.

Per library:
  pM in well        = 10^((Cq - intercept) / slope), from the standard curve (log10 pM vs Cq)
  undiluted pM      = pM in well x dilution fold        (per dilution, mean of replicates)
  size-corrected nM = mean undiluted pM x standard size / library size / 1000
  fluorometric nM   = ng/ul x 1e6 / (660 x library size)  (if given), and the qPCR / fluorometric ratio
"""
import argparse, csv, glob, json, os, re, sys
import numpy as np
import yaml


def norm_well(w):
    m = re.fullmatch(r'\s*([A-Pa-p])0*(\d+)\s*', str(w))
    return f'{m.group(1).upper()}{int(m.group(2))}' if m else str(w).strip()


def expand(spec):
    if isinstance(spec, list):
        return [norm_well(w) for s in spec for w in expand(s)] if any('-' in str(s) for s in spec) else [norm_well(s) for s in spec]
    m = re.fullmatch(r'([A-P])(\d+)-([A-P])(\d+)', str(spec).strip())
    if m:
        R = 'ABCDEFGHIJKLMNOP'
        return [f'{r}{c}' for c in range(int(m.group(2)), int(m.group(4)) + 1)
                for r in R[R.index(m.group(1)):R.index(m.group(3)) + 1]]
    return [norm_well(spec)]


def find(d, pattern):
    hits = sorted(glob.glob(os.path.join(d, pattern)))
    return hits[0] if hits else None


def read_cq(cfg, export_dir):
    p = cfg.get('cq_file')
    p = os.path.join(export_dir, p) if p else find(export_dir, '*Quantification Cq Results.csv')
    if p and p.endswith('.csv'):
        out = {}
        with open(p, newline='', encoding='utf-8-sig') as fh:
            for r in csv.DictReader(fh):
                try:
                    out[norm_well(r['Well'])] = float(r['Cq'])
                except (ValueError, KeyError, TypeError):
                    out[norm_well(r.get('Well', ''))] = float('nan')
        return out, p
    import pandas as pd   # QuantStudio .xls/.xlsx
    p = p or find(export_dir, '*.xls*')
    raw = pd.read_excel(p, sheet_name='Results', header=None)
    i = raw.index[raw[0].astype(str).str.strip() == 'Well'][0]
    df = raw.iloc[i + 1:]; df.columns = raw.iloc[i].astype(str).str.strip()
    return {norm_well(w): pd.to_numeric(c, errors='coerce') for w, c in zip(df['Well Position'], df[cfg.get('ct_column', 'CT')])}, p


def read_melt(export_dir):
    p = find(export_dir, '*Melt Curve Peak Results.csv')
    out = {}
    if p:
        with open(p, newline='', encoding='utf-8-sig') as fh:
            for r in csv.DictReader(fh):
                try:
                    out.setdefault(norm_well(r['Well']), []).append(float(r['Melt Temperature']))
                except (ValueError, TypeError):
                    pass
    return out


def read_amp(export_dir):
    p = find(export_dir, '*Quantification Amplification Results*.csv')
    if not p:
        return None
    with open(p, newline='', encoding='utf-8-sig') as fh:
        rd = csv.reader(fh); hdr = next(rd)
        ci = hdr.index('Cycle'); wells = [norm_well(h) for h in hdr[ci + 1:]]
        rows = [r for r in rd if r and r[ci]]
    cyc = np.array([float(r[ci]) for r in rows])
    return cyc, {w: np.array([float(r[ci + 1 + j]) for r in rows]) for j, w in enumerate(wells)}


def run(cfg_path, export_dir=None, outdir=None, plots=True):
    cfg = yaml.safe_load(open(cfg_path))
    base = os.path.dirname(os.path.abspath(cfg_path))
    export_dir = export_dir or os.path.join(base, cfg.get('export_dir', '.'))
    outdir = outdir or os.path.join(base, cfg.get('outdir', 'library_quant_out'))
    os.makedirs(outdir, exist_ok=True)
    cq, cq_path = read_cq(cfg, export_dir)
    st = cfg['standards']
    sw = expand(st['wells']) if 'wells' in st else None
    lv = st['levels'] if 'levels' in st else [dict(pM=float(p), wells=[w]) for p, w in zip(st['pM'], sw)]
    pts = [(float(l['pM']), cq.get(norm_well(w), np.nan)) for l in lv for w in expand(l['wells'])]
    pts = [(p, c) for p, c in pts if np.isfinite(c)]
    x = np.log10([p for p, _ in pts]); y = np.array([c for _, c in pts])
    slope, icpt = np.polyfit(x, y, 1)
    pred = slope * x + icpt
    r2 = float(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    eff = float(10 ** (-1 / slope) - 1)
    fit = dict(slope=float(slope), intercept=float(icpt), r2=r2, efficiency=eff, n=len(pts),
               cq_range=[float(y.min()), float(y.max())], standard_size_bp=float(st.get('size_bp', 452)))
    pm = lambda c: 10 ** ((c - icpt) / slope)
    melt = read_melt(export_dir)
    std_tm = [t for l in lv for w in expand(l['wells']) for t in melt.get(norm_well(w), [])[:1]]
    warn = []
    if not (0.9 <= eff <= 1.1):
        warn.append(f'standard curve efficiency {eff*100:.1f}% outside 90-110%')
    if r2 < 0.99:
        warn.append(f'standard curve R2 {r2:.4f} < 0.99')
    libs = []
    for s in cfg['libraries']:
        size = float(s['size_bp'])
        per = []
        for d in s['dilutions']:
            c = np.array([cq.get(w, np.nan) for w in expand(d['wells'])], float)
            ok = c[np.isfinite(c)]
            mc = float(ok.mean()) if len(ok) else float('nan')
            per.append(dict(fold=float(d['fold']), cq_mean=mc, cq_sd=float(ok.std(ddof=1)) if len(ok) > 1 else float('nan'),
                            n=len(ok), pM_well=pm(mc), pM_undiluted=pm(mc) * float(d['fold']),
                            in_range=bool(fit['cq_range'][0] - 0.5 <= mc <= fit['cq_range'][1] + 0.5) if len(ok) else False))
        use = [p for p in per if p['in_range']] or per
        und = float(np.mean([p['pM_undiluted'] for p in use]))
        nM = und * fit['standard_size_bp'] / size / 1000
        rec = dict(name=s['name'], size_bp=size, nM_size_corrected=nM, pM_undiluted_uncorrected=und,
                   dilutions_used=len(use))
        for i, p in enumerate(per):
            rec.update({f'd{i+1}_fold': p['fold'], f'd{i+1}_cq': p['cq_mean'], f'd{i+1}_cq_sd': p['cq_sd'],
                        f'd{i+1}_pM_undiluted': p['pM_undiluted'], f'd{i+1}_in_range': p['in_range']})
        if len(per) >= 2 and all(np.isfinite([per[0]['cq_mean'], per[1]['cq_mean']])):
            dcq = per[1]['cq_mean'] - per[0]['cq_mean']
            exp = np.log10(per[1]['fold'] / per[0]['fold']) * -slope
            rec['delta_cq'] = dcq; rec['delta_cq_expected'] = float(exp)
            rec['dilution_ratio'] = per[0]['pM_undiluted'] / per[1]['pM_undiluted']
            if abs(dcq - exp) > 0.5:
                warn.append(f"{s['name']}: dilutions differ by {dcq:.2f} cycles, expected {exp:.2f} (inhibition or pipetting)")
        for p in per:
            if p['n'] and np.isfinite(p['cq_sd']) and p['cq_sd'] > 0.2:
                warn.append(f"{s['name']} 1:{p['fold']:g}: replicate Cq SD {p['cq_sd']:.2f} > 0.2")
            if p['n'] and not p['in_range']:
                warn.append(f"{s['name']} 1:{p['fold']:g}: Cq {p['cq_mean']:.2f} outside the standard range, not used")
        if s.get('fluorometric_ng_ul') is not None:
            fl = float(s['fluorometric_ng_ul']) * 1e6 / (660 * size)
            rec['fluorometric_nM'] = fl; rec['qpcr_over_fluorometric'] = nM / fl if fl else float('nan')
        tms = [t for d in s['dilutions'] for w in expand(d['wells']) for t in melt.get(w, [])[:1]]
        if tms:
            rec['melt_tm_mean'] = float(np.mean(tms))
            if any(len(melt.get(w, [])) > 1 for d in s['dilutions'] for w in expand(d['wells'])):
                warn.append(f"{s['name']}: more than one melt peak in some wells")
        libs.append(rec)
    for c in cfg.get('ntc', []) and [cfg['ntc']] or []:
        for w in expand(c):
            if np.isfinite(cq.get(w, np.nan)) and cq[w] < fit['cq_range'][1] + 3:
                warn.append(f'NTC {w} Cq {cq[w]:.1f} close to the lowest standard')
    with open(os.path.join(outdir, 'libraries.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(dict.fromkeys(k for r in libs for k in r))); w.writeheader(); w.writerows(libs)
    with open(os.path.join(outdir, 'standards.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['pM', 'Cq', 'residual'])
        for (p, c), r in zip(pts, y - pred):
            w.writerow([p, c, r])
    json.dump(dict(cq_file=cq_path, standard_curve=fit, standard_tm=std_tm, libraries=libs, warnings=warn),
              open(os.path.join(outdir, 'summary.json'), 'w'), indent=1, default=float)
    if plots:
        plot(cfg, fit, pts, libs, read_amp(export_dir), lv, outdir)
    return fit, libs, warn


sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', '..', '_shared'))
import labstyle


def plot(cfg, fit, pts, libs, amp, lv, outdir):
    """Rules: _shared/FIGURES.md. Two figures: QC (curve + amplification) and concentrations."""
    fig, ax = labstyle.figure('double', ncols=2)
    A, B = ax
    x = np.log10([p for p, _ in pts]); y = [c for _, c in pts]
    xs = np.linspace(x.min() - 0.5, x.max() + 0.5, 20)
    A.plot(xs, fit['slope'] * xs + fit['intercept'], color=labstyle.GREY, lw=0.7)
    A.plot(x, y, 'o', mfc='white', mec=labstyle.INK, ms=4, label='Standards')
    first = True
    for r in libs:
        for i in (1, 2):
            if f'd{i}_cq' in r and np.isfinite(r[f'd{i}_cq']):
                A.plot((r[f'd{i}_cq'] - fit['intercept']) / fit['slope'], r[f'd{i}_cq'], 'o', color=labstyle.BLUE,
                       ms=3, label='Libraries' if first else None)
                first = False
    A.set_xlabel('log$_{10}$ pM in well'); A.set_ylabel('Cq'); A.set_title('Standard curve')
    A.legend(loc='lower left')
    A.text(0.97, 0.95, f"efficiency {fit['efficiency']*100:.1f}%, R² {fit['r2']:.4f}", transform=A.transAxes,
           ha='right', va='top', fontsize=6, color=labstyle.MUTED)
    if amp:
        cyc, d = amp
        sw = {w for l in lv for w in expand(l['wells'])}
        lw = {w for s_ in cfg['libraries'] for dd in s_['dilutions'] for w in expand(dd['wells'])}
        for w, v in d.items():
            if w in sw or w in lw:
                B.plot(cyc, v, color=labstyle.INK if w in sw else labstyle.BLUE, lw=0.5, alpha=0.7)
        B.set_xlabel('Cycle'); B.set_ylabel('RFU')
    else:
        labstyle.empty(B, 'amplification data not in export')
    B.set_title('Amplification')
    labstyle.save(fig, os.path.join(outdir, 'library_quant_qc'), caption=(
        'Standard curve: Cq vs log10 pM of the DNA standards (open), least-squares fit; libraries (blue) placed on '
        'the curve at each dilution. Amplification: black = standards, blue = libraries.'))
    fig, ax = labstyle.figure('single' if len(libs) <= 6 else 'double', ncols=1)
    xi = np.arange(len(libs)); two = any('fluorometric_nM' in r for r in libs); w = 0.38 if two else 0.6
    bars = ax.bar(xi - (w / 2 if two else 0), [r['nM_size_corrected'] for r in libs], w, color=labstyle.BLUE,
                  label='qPCR, size-corrected')
    if two:
        ax.bar(xi + w / 2, [r.get('fluorometric_nM', np.nan) for r in libs], w, color=labstyle.STAGES[2],
               label='Fluorometric')
        ax.legend(loc='upper left')
    ax.set_yscale('log'); ax.set_ylabel('Undiluted library (nM)')
    for b in bars:
        v = b.get_height()
        ax.annotate(f'{v:.0f}' if v >= 10 else f'{v:.2g}', (b.get_x() + b.get_width() / 2, b.get_height()), xytext=(0, 1),
                    textcoords='offset points', ha='center', va='bottom', fontsize=6)
    ax.set_xticks(xi); ax.set_xticklabels([r['name'] for r in libs], rotation=45 if len(libs) > 4 else 0,
                                          ha='right' if len(libs) > 4 else 'center')
    labstyle.save(fig, os.path.join(outdir, 'library_quant'), caption=(
        'Undiluted library concentration. qPCR: mean of in-range dilutions x fold, x standard size / library size. '
        'Fluorometric: ng/µl x 1e6 / (660 x library size).'))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('config'); ap.add_argument('--export-dir'); ap.add_argument('--outdir'); ap.add_argument('--no-plots', action='store_true')
    a = ap.parse_args()
    fit, libs, warn = run(a.config, a.export_dir, a.outdir, not a.no_plots)
    print(f"standard curve: slope {fit['slope']:.3f}  efficiency {fit['efficiency']*100:.1f}%  R2 {fit['r2']:.5f}")
    for r in libs:
        extra = f"  fluorometric {r['fluorometric_nM']:.2f} nM (ratio {r['qpcr_over_fluorometric']:.2f})" if 'fluorometric_nM' in r else ''
        print(f"{r['name']:16s} {r['nM_size_corrected']:9.3f} nM  ({r['size_bp']:.0f} bp){extra}")
    for x in warn:
        print('WARNING', x)


if __name__ == '__main__':
    main()
