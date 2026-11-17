#!/usr/bin/env python3
"""Absolute qPCR quantification from a QuantStudio export and a YAML plate layout.

    python3 qpcr_quant.py config.yaml [--export file.xls] [--outdir DIR]

Supports: standard curve (Ct vs log10 copies), Ct -> copies, melt-peak QC,
single-bead (digital-style) wells with automatic bead/no-bead threshold, background
subtraction and Poisson fit, and bulk release normalised per bead.

Outputs (outdir): wells.csv, standards.csv, single_bead_<group>.csv, poisson_<group>.csv,
bulk.csv, summary.json, figures (qpcr_qc.png, qpcr_results.png).
See ../references/method.md for formulas and ../references/lessons.md for pitfalls.
"""
import argparse, json, os, re, sys
from math import exp, factorial
import numpy as np
import pandas as pd
import yaml

AVOGADRO = 6.02214076e23
ROWS = 'ABCDEFGHIJKLMNOP'


# ------------------------------------------------------------------ input
def read_sheet(path, name):
    """QuantStudio .xls/.xlsx: header row is the first row whose first cell is 'Well'."""
    raw = pd.read_excel(path, sheet_name=name, header=None)
    i = raw.index[raw[0].astype(str).str.strip() == 'Well'][0]
    df = raw.iloc[i + 1:].copy()
    df.columns = raw.iloc[i].astype(str).str.strip()
    return df.dropna(how='all').reset_index(drop=True)


def expand_wells(spec):
    """'A1-P6' (rectangle), 'A9,C9', ['A9','B9'], or 'col:7' / 'rows:ACE col:8'."""
    if isinstance(spec, list):
        out = []
        for s in spec:
            out += expand_wells(s)
        return out
    spec = str(spec).strip()
    m = re.fullmatch(r'([A-P])(\d+)-([A-P])(\d+)', spec)
    if m:
        r0, c0, r1, c1 = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
        return [f'{r}{c}' for c in range(c0, c1 + 1) for r in ROWS[ROWS.index(r0):ROWS.index(r1) + 1]]
    m = re.fullmatch(r'rows:([A-P]+)\s+col:(\d+)', spec)
    if m:
        return [f'{r}{m.group(2)}' for r in m.group(1)]
    if ',' in spec:
        return expand_wells([s.strip() for s in spec.split(',')])
    return [spec]


def load_wells(export, ct_col='CT'):
    r = read_sheet(export, 'Results')
    r['well'] = r['Well Position'].astype(str).str.strip()
    r['Ct'] = pd.to_numeric(r[ct_col], errors='coerce')
    try:
        m = read_sheet(export, 'Melt Curve Result')
        m['well'] = m['Well Position'].astype(str).str.strip()
        m['Tm'] = pd.to_numeric(m['Tm'], errors='coerce')
        r = r.merge(m[['well', 'Tm']].drop_duplicates('well'), on='well', how='left')
    except Exception:
        r['Tm'] = np.nan
    return r[['well', 'Ct', 'Tm']]


def try_sheet(export, name):
    try:
        return read_sheet(export, name)
    except Exception:
        return None


# ------------------------------------------------------------------ core
def standard_copies(st):
    """Copies per well for each standard level.
    Either explicit {levels: [{wells, copies}]} or a stock + serial dilution:
    {stock_molar, stock_dilution, template_ul, fold, levels: [[wells], ...]} (first level = stock_dilution)."""
    out = {}
    if 'top_copies' in st or 'stock_molar' in st:
        top = st.get('top_copies')
        top = None if top is None else float(top)  # PyYAML reads 1e6 (no sign) as a string
        if top is None:
            top = float(st['stock_molar']) / float(st.get('stock_dilution', 1)) * float(st['template_ul']) * 1e-6 * AVOGADRO
        for i, ws in enumerate(st['levels']):
            for w in expand_wells(ws):
                out[w] = top / float(st.get('fold', 10)) ** i
    else:
        for lv in st['levels']:
            for w in expand_wells(lv['wells']):
                out[w] = float(lv['copies'])
    return out


def fit_standards(wells, st):
    cp = standard_copies(st)
    s = wells[wells.well.isin(cp)].copy()
    s['copies'] = s.well.map(cp); s['log10_copies'] = np.log10(s.copies)
    s = s.dropna(subset=['Ct'])
    slope, icpt = np.polyfit(s.log10_copies, s.Ct, 1)
    pred = slope * s.log10_copies + icpt
    r2 = 1 - ((s.Ct - pred) ** 2).sum() / ((s.Ct - s.Ct.mean()) ** 2).sum()
    eff = 10 ** (-1 / slope) - 1
    s['residual'] = s.Ct - pred
    return s, dict(slope=float(slope), intercept=float(icpt), r2=float(r2), efficiency=float(eff),
                   n=int(len(s)), range_copies=[float(s.copies.min()), float(s.copies.max())])


def auto_threshold(ct, min_gap=1.5):
    """Largest gap in sorted Ct values; returns midpoint, or None if no gap >= min_gap."""
    v = np.sort(ct[np.isfinite(ct)])
    if len(v) < 3:
        return None
    gaps = np.diff(v); i = int(np.argmax(gaps))
    return float((v[i] + v[i + 1]) / 2) if gaps[i] >= min_gap else None


def poisson_table(k, lam_planned):
    n = len(k); lam = -np.log(max((k == 0).mean(), 1e-9))
    def exp_counts(l):
        p = [exp(-l) * l ** i / factorial(i) for i in range(3)]
        return np.array(p + [1 - sum(p)]) * n
    obs = np.array([(k == i).sum() for i in range(3)] + [(k >= 3).sum()])
    return pd.DataFrame(dict(beads=['0', '1', '2', '>=3'], observed=obs,
                             expected_planned=exp_counts(lam_planned) if lam_planned else np.nan,
                             expected_fitted=exp_counts(lam))), float(lam)


def single_bead(wells, g, q, unit=None):
    ws = expand_wells(g['wells'])
    d = wells[wells.well.isin(ws)].copy()
    thr = g.get('threshold', 'auto')
    if thr == 'auto':
        thr = auto_threshold(d.Ct.values, g.get('min_gap', 1.5))
        if thr is None:
            sys.exit(f"group {g['name']}: no Ct gap >= min_gap; set threshold by hand")
    d['copies'] = q(d.Ct)
    d['bead'] = d.Ct < thr
    bg = float(d.loc[~d.bead, 'copies'].median()) if (~d.bead).any() else 0.0
    d['copies_minus_bg'] = d.copies - bg
    # unit = copies of one bead; from this block's median unless given (use a low-lambda block)
    unit = unit or float(d.loc[d.bead, 'copies_minus_bg'].median())
    d['beads_est'] = np.where(d.bead, np.maximum(1, np.rint(d.copies_minus_bg / unit)), 0).astype(int)
    d['per_bead'] = np.where(d.bead, d.copies_minus_bg / d.beads_est.clip(lower=1), np.nan)
    per = d.loc[d.bead, 'per_bead']
    tab, lam = poisson_table(d.beads_est.values, g.get('lambda_planned'))
    res = dict(threshold=float(thr), wells=int(len(d)), bead_wells=int(d.bead.sum()),
               background_copies=bg, per_bead_median=float(per.median()) if len(per) else None,
               per_bead_mean=float(per.mean()) if len(per) else None,
               per_bead_cv=float(per.std() / per.mean()) if len(per) > 1 else None,
               lambda_fitted=lam, lambda_planned=g.get('lambda_planned'), unit_copies=unit)
    return d, tab, res


def bulk(wells, g, q):
    rows = []
    for smp in g['samples']:
        for dil in smp['dilutions']:
            for w in expand_wells(dil['wells']):
                ct = wells.loc[wells.well == w, 'Ct']
                ct = float(ct.iloc[0]) if len(ct) else np.nan
                beads = smp.get('beads')
                if beads is None:
                    b = smp['bead_count']
                    beads = b['per_ul'] * b.get('dilution', 1) * b['volume_ul']
                c = q(ct)
                total = c / g['template_ul'] * dil['fold'] * g['eluate_ul']
                rows.append(dict(sample=smp['name'], label=smp.get('label', smp['name']), dilution=dil['fold'], well=w,
                                 Ct=ct, copies=c, eluate_total=total, beads=beads, per_bead=total / beads))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ plots
def plots(cfg, wells, s, fit, sb, bk, outdir, amp=None, melt=None):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans', 'axes.spines.top': False,
                         'axes.spines.right': False, 'legend.frameon': False, 'legend.fontsize': 6.5})
    dark, acc, grey, light = '#2b2b2b', '#c0392b', '#9a9a9a', '#cfcfcf'
    rng = np.random.default_rng(1)
    # ---- QC figure: standard curve, amplification, melt, plate map
    fig, ax = plt.subplots(1, 4, figsize=(12, 2.9))
    A, B, C, D = ax
    x = np.linspace(s.log10_copies.min() - 0.3, s.log10_copies.max() + 0.3, 10)
    A.plot(x, fit['slope'] * x + fit['intercept'], color=grey, lw=0.8); A.scatter(s.log10_copies, s.Ct, s=12, color=dark, zorder=3)
    A.set_xlabel('log$_{10}$ copies / well'); A.set_ylabel('Ct'); A.set_title('A  Standard curve', loc='left')
    A.text(0.97, 0.95, f"efficiency {fit['efficiency']*100:.0f}%\nR² {fit['r2']:.4f}", transform=A.transAxes, ha='right', va='top', fontsize=7)
    if amp is not None:
        amp = amp.copy(); amp['Cycle'] = pd.to_numeric(amp['Cycle']); amp['dRn'] = pd.to_numeric(amp['Delta Rn'])
        for w, sx in amp.groupby('Well Position'):
            B.plot(sx.Cycle, sx.dRn.clip(lower=1e-3), lw=0.4, color=dark if w in set(s.well) else acc, alpha=0.5)
        B.set_yscale('log'); B.set_xlabel('cycle'); B.set_ylabel('ΔRn')
    B.set_title('B  Amplification (dark = standards)', loc='left')
    if melt is not None:
        melt = melt.copy(); melt['T'] = pd.to_numeric(melt['Temperature']); melt['d'] = pd.to_numeric(melt['Derivative'])
        for w, sx in melt.groupby('Well Position'):
            C.plot(sx['T'], sx.d / sx.d.max(), lw=0.4, color=dark if w in set(s.well) else acc, alpha=0.5)
        C.set_xlabel('°C'); C.set_ylabel('−dF/dT (norm.)')
    tm = cfg.get('melt', {}).get('product_tm')
    if tm:
        C.axvline(tm, color=grey, lw=0.6, ls='--')
    C.set_title('C  Melt curves', loc='left')
    M = np.full((16, 24), np.nan)
    for _, w in wells.iterrows():
        m = re.fullmatch(r'([A-P])(\d+)', w.well)
        if m and np.isfinite(w.Ct):
            M[ROWS.index(m.group(1)), int(m.group(2)) - 1] = w.Ct
    cols = np.where(np.isfinite(M).any(0))[0]
    ncol = int(cols.max()) + 1 if len(cols) else 24
    im = D.imshow(M[:, :ncol], cmap='magma', vmin=np.nanmin(M), vmax=40, aspect='auto')
    D.set_xticks(range(ncol)); D.set_xticklabels(range(1, ncol + 1), fontsize=5.5)
    D.set_yticks(range(16)); D.set_yticklabels(list(ROWS), fontsize=5.5); D.tick_params(length=0)
    fig.colorbar(im, ax=D, fraction=0.05, pad=0.02).set_label('Ct'); D.set_title('D  Plate map (Ct)', loc='left')
    fig.tight_layout(); fig.savefig(os.path.join(outdir, 'qpcr_qc.png'), dpi=220, bbox_inches='tight'); plt.close(fig)
    # ---- results figure: Ct by group, Poisson, per bead
    if not sb and bk is None:
        return
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.0))
    E, F, G = ax
    for i, (name, (d, tab, res)) in enumerate(sb.items()):
        col = np.where(d.bead, acc, grey)
        E.scatter(i + rng.uniform(-0.18, 0.18, len(d)), d.Ct, s=7, color=col, lw=0)
        E.hlines(res['threshold'], i - 0.3, i + 0.3, color=dark, lw=0.6, ls='--')
    for j, cg in enumerate(cfg.get('controls', []), start=len(sb)):
        d = wells[wells.well.isin(expand_wells(cg['wells']))]
        E.scatter(j + rng.uniform(-0.18, 0.18, len(d)), d.Ct.fillna(40), s=7, color=grey, lw=0)
    labels = list(sb) + [c['name'] for c in cfg.get('controls', [])]
    E.set_xticks(range(len(labels))); E.set_xticklabels(labels, fontsize=6.5, rotation=20)
    E.set_ylabel('Ct'); E.invert_yaxis(); E.set_title('E  Ct per well (red = bead)', loc='left')
    if sb:
        name, (d, tab, res) = next(iter(sb.items()))
        k = np.arange(4)
        F.bar(k - 0.2, tab.observed, 0.38, color=acc, label='observed')
        if res['lambda_planned']:
            F.bar(k + 0.2, tab.expected_planned, 0.38, color=light, label=f"Poisson λ planned {res['lambda_planned']:g}")
        F.plot(k + 0.2, tab.expected_fitted, '_', ms=10, mew=1.4, color=dark, label=f"Poisson λ fitted {res['lambda_fitted']:.2f}")
        F.set_xticks(k); F.set_xticklabels(tab.beads); F.set_xlabel('beads per well'); F.set_ylabel('wells')
        F.legend(loc='upper right'); F.set_title(f'F  Beads per well, {name}', loc='left')
    xi = 0
    for name, (d, tab, res) in sb.items():
        v = d.loc[d.bead, 'per_bead']
        G.scatter(xi + rng.uniform(-0.15, 0.15, len(v)), v, s=8, color=acc, lw=0); G.hlines(v.median(), xi - 0.25, xi + 0.25, color=dark)
        xi += 1
    ticks = list(sb)
    if bk is not None and len(bk):
        for sn, dd in bk.groupby('sample', sort=False):
            G.scatter(np.full(len(dd), xi) + rng.uniform(-0.12, 0.12, len(dd)), dd.per_bead, s=12, facecolor='white', edgecolor=dark)
            G.hlines(dd.per_bead.mean(), xi - 0.25, xi + 0.25, color=dark); ticks.append(dd.label.iloc[0]); xi += 1
    G.set_yscale('log'); G.set_xticks(range(len(ticks))); G.set_xticklabels(ticks, fontsize=6.5, rotation=20)
    G.set_ylabel('copies per bead'); G.set_title('G  Copies per bead (red single, open bulk)', loc='left')
    fig.tight_layout(); fig.savefig(os.path.join(outdir, 'qpcr_results.png'), dpi=220, bbox_inches='tight'); plt.close(fig)


# ------------------------------------------------------------------ main
def run(cfg_path, export=None, outdir=None, make_plots=True):
    cfg = yaml.safe_load(open(cfg_path))
    base = os.path.dirname(os.path.abspath(cfg_path))
    export = export or os.path.join(base, cfg['export'])
    outdir = outdir or os.path.join(base, cfg.get('outdir', 'qpcr_quant_out'))
    os.makedirs(outdir, exist_ok=True)
    wells = load_wells(export, cfg.get('ct_column', 'CT'))
    s, fit = fit_standards(wells, cfg['standards'])
    q = lambda ct: 10 ** ((np.asarray(ct, float) - fit['intercept']) / fit['slope'])
    wells['copies'] = q(wells.Ct)
    mc = cfg.get('melt', {})
    if mc.get('product_tm') is not None:
        wells['melt_ok'] = (wells.Tm - mc['product_tm']).abs() <= mc.get('tolerance', 1.5)
    s.to_csv(os.path.join(outdir, 'standards.csv'), index=False)
    summary = dict(export=export, standard_curve=fit, single_bead={}, controls={})
    sb = {}
    unit = None
    for g in cfg.get('single_bead', []):
        d, tab, res = single_bead(wells, g, q, unit=unit if g.get('unit_from_first', True) else None)
        unit = unit or res['unit_copies']
        sb[g['name']] = (d, tab, res); summary['single_bead'][g['name']] = res
        d.to_csv(os.path.join(outdir, f"single_bead_{g['name']}.csv"), index=False)
        tab.to_csv(os.path.join(outdir, f"poisson_{g['name']}.csv"), index=False)
    bk = None
    if cfg.get('bulk'):
        bk = bulk(wells, cfg['bulk'], q)
        bk.to_csv(os.path.join(outdir, 'bulk.csv'), index=False)
        agg = bk.groupby('sample').per_bead.agg(['mean', 'min', 'max', 'count'])
        summary['bulk'] = {k: {c: float(v) for c, v in r.items()} for k, r in agg.iterrows()}
        dd = bk.groupby(['sample', 'dilution']).Ct.mean().unstack()
        if dd.shape[1] == 2:  # dilution check: expected delta Ct = log10(fold ratio) / -slope ... per 10x ~ 3.32
            f0, f1 = sorted(dd.columns)
            exp_d = np.log10(f1 / f0) * -fit['slope']
            summary['bulk_dilution_check'] = {k: dict(delta_ct=float(r[f1] - r[f0]), expected=float(exp_d)) for k, r in dd.iterrows()}
    for cg in cfg.get('controls', []):
        d = wells[wells.well.isin(expand_wells(cg['wells']))]
        summary['controls'][cg['name']] = dict(n=int(len(d)), min_ct=float(d.Ct.min()) if d.Ct.notna().any() else None,
                                               tm=[float(x) for x in d.Tm.dropna()])
    wells.to_csv(os.path.join(outdir, 'wells.csv'), index=False)
    json.dump(summary, open(os.path.join(outdir, 'summary.json'), 'w'), indent=1, default=float)
    if make_plots:
        plots(cfg, wells, s, fit, sb, bk, outdir, try_sheet(export, 'Amplification Data'), try_sheet(export, 'Melt Curve Raw Data'))
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('config'); ap.add_argument('--export'); ap.add_argument('--outdir'); ap.add_argument('--no-plots', action='store_true')
    a = ap.parse_args()
    s = run(a.config, a.export, a.outdir, not a.no_plots)
    f = s['standard_curve']
    print(f"standard curve: slope {f['slope']:.3f}  intercept {f['intercept']:.2f}  R2 {f['r2']:.4f}  efficiency {f['efficiency']*100:.1f}%")
    for k, r in s['single_bead'].items():
        print(f"{k}: threshold Ct {r['threshold']:.2f}, bead wells {r['bead_wells']}/{r['wells']}, "
              f"background {r['background_copies']:.3g}, per bead median {r['per_bead_median']:.3g}, lambda {r['lambda_fitted']:.2f}")
    for k, r in s.get('bulk', {}).items():
        print(f"bulk {k}: per bead mean {r['mean']:.3g} (range {r['min']:.3g}-{r['max']:.3g})")


if __name__ == '__main__':
    main()
