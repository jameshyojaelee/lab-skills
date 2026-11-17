#!/usr/bin/env python3
"""Quantify bands on a denaturing PAGE / agarose gel image from a YAML config.

    python3 gel_quant.py config.yaml [--image other.tif] [--outdir DIR] [--compare second.tif]

Outputs (in outdir): band_fits.csv, metrics.csv, amounts.csv, summary.json,
gel_annotated.png, metrics.png, amounts.png (plots only when configured).

Method (see ../references/method.md for the reasoning):
  1. Lane profile: median across a strip of width 2*halfwidth+1 px, per row.
  2. Background: rolling minimum then rolling mean (window `background_window`), subtracted.
  3. Per-lane shift: the lane's anchor (max or centroid inside a window) minus its reference y.
  4. One Gaussian per expected band, all fitted together; area = height * sigma * sqrt(2*pi).
  5. Molar = area / length (nt), assuming stain signal proportional to length.
  6. Metrics = sum(molar of numerator bands) / sum(molar of denominator bands).
  7. Clipping: pixels >= saturation within +-2 rows of a fitted centre flag that band.
     Bound direction: clipped band in numerator -> value is a lower bound ('>='),
     clipped band only in denominator -> upper bound ('<='), both -> unreliable ('~').
"""
import argparse, csv, json, os, re, sys
import numpy as np
import yaml
from PIL import Image
from scipy.ndimage import uniform_filter1d, minimum_filter1d
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


# ---------------------------------------------------------------- image IO
def read_scn(path):
    """Bio-Rad Image Lab .scn: MIME-like container, ImageData block is little-endian uint16."""
    raw = open(path, 'rb').read()
    i = raw.find(b'Content-Description: ImageData')
    j = raw.find(b'\r\n\r\n', i) + 4
    if j < 4:
        j = raw.find(b'\n\n', i) + 2
    n = int(re.search(rb'Content-Length: (\d+)', raw[raw.rfind(b'Content-Type', 0, i):i]).group(1))
    w = int(re.search(rb'<size_pix width="(\d+)"', raw).group(1))
    h = int(re.search(rb'<size_pix width="\d+" height="(\d+)"', raw).group(1))
    mx = int(re.search(rb'<scanner max_value="(\d+)"', raw).group(1))
    img = np.frombuffer(raw[j:j + n], dtype='<u2').reshape(h, w).astype(float)
    return img, float(mx)


def load_image(path):
    """Returns (2D float array, saturation value)."""
    if path.lower().endswith('.scn'):
        return read_scn(path)
    im = Image.open(path)
    a = np.array(im)
    if a.ndim == 3:
        a = a[..., 0]
    if a.dtype == np.uint8:
        sat = 254.0  # Image Lab 8-bit exports clip at 254; 255 is used for annotation pixels
    elif a.dtype == np.uint16:
        sat = float(a.max()) if a.max() >= 65535 else 65535.0
    else:
        sat = float(np.nanmax(a))
    return a.astype(float), sat


# ---------------------------------------------------------------- helpers
def gsum(y, *p):
    o = np.zeros_like(y, float)
    for k in range(0, len(p), 3):
        o += p[k] * np.exp(-0.5 * ((y - p[k + 1]) / p[k + 2]) ** 2)
    return o


def size_to_y(nt, cal):
    """Log-linear size calibration y = a + b*log10(nt) from >=2 (nt, y) points."""
    pts = np.array(cal, float)
    b, a = np.polyfit(np.log10(pts[:, 0]), pts[:, 1], 1)
    return a + b * np.log10(nt)


def detect_lanes(a, rows, n, min_dist):
    col = a[rows[0]:rows[1], :].mean(0)
    col = uniform_filter1d(col, 5)
    pk, pr = find_peaks(col, distance=min_dist, prominence=1)
    if len(pk) < n:
        sys.exit(f'lane detection found {len(pk)} lanes, config lists {n}; give x positions')
    top = pk[np.argsort(pr['prominences'])[::-1][:n]]
    return sorted(int(x) for x in top)


class Gel:
    def __init__(self, cfg, image=None):
        self.cfg = cfg
        self.path = image or cfg['image']
        self.a, sat = load_image(self.path)
        self.sat = float(cfg.get('saturation', sat))
        p = cfg.get('profile', {})
        self.y0, self.y1 = p.get('rows', [0, self.a.shape[0]])
        self.hw = p.get('halfwidth', 4)
        self.bgw = p.get('background_window', 41)
        self.lanes = cfg['lanes']
        if any('x' not in L for L in self.lanes):
            xs = detect_lanes(self.a, [self.y0, self.y1], len(self.lanes), p.get('min_lane_spacing', 15))
            for L, x in zip(self.lanes, xs):
                L.setdefault('x', x)
        # band table
        cal = cfg.get('calibration')
        self.bands = []
        for b in cfg['bands']:
            b = dict(b)
            if 'y' not in b:
                if not cal:
                    sys.exit(f"band {b['id']} has no y and no calibration")
                b['y'] = float(size_to_y(b['length'], cal))
            self.bands.append(b)
        self.blen = {b['id']: b.get('length') for b in self.bands}

    def profile(self, x):
        p = np.median(self.a[self.y0:self.y1, x - self.hw:x + self.hw + 1], axis=1)
        bg = uniform_filter1d(minimum_filter1d(p, self.bgw), self.bgw)
        return p - bg

    def lane_shift(self, L, d, yy):
        an = dict(self.cfg.get('anchor', {}))
        an.update(L.get('anchor', {}))
        if not an:
            return 0.0
        w = (yy >= an['lo']) & (yy <= an['hi'])
        if an.get('mode', 'max') == 'centroid':
            dd = d[w].clip(0)
            return float((dd * yy[w]).sum() / dd.sum() - an['ref'])
        return float(yy[w][np.argmax(d[w])] - an['ref'])

    def fit(self):
        f = self.cfg.get('fit', {})
        tol, s0, smin, smax = f.get('centre_tolerance', 2), f.get('sigma_init', 2.5), f.get('sigma_min', 1.5), f.get('sigma_max', 4.5)
        r0, r1 = f.get('range', [self.y0, self.y1])
        yy = np.arange(self.y0, self.y1).astype(float)
        out = {}
        for L in self.lanes:
            if L.get('role', 'sample') in ('marker', 'blank', 'skip'):
                continue
            x = int(L['x'])
            d = self.profile(x)
            s = self.lane_shift(L, d, yy)
            p0, lo, hi = [], [], []
            for b in self.bands:
                c = b['y'] + s
                i = int(np.clip(round(c) - self.y0, 0, len(d) - 1))
                p0 += [max(d[i], 0.5), c, s0]; lo += [0, c - tol, smin]; hi += [np.inf, c + tol, smax]
            m = (yy >= r0 + s) & (yy <= r1 + s)
            p, _ = curve_fit(gsum, yy[m], d[m], p0=p0, bounds=(lo, hi), maxfev=50000)
            raw = self.a[self.y0:self.y1, x - self.hw:x + self.hw + 1]
            res = {}
            for j, b in enumerate(self.bands):
                c = p[3 * j + 1]
                k = int(round(c)) - self.y0
                res[b['id']] = dict(area=float(p[3 * j] * p[3 * j + 2] * np.sqrt(2 * np.pi)), y=float(c),
                                    sigma=float(p[3 * j + 2]), clip=int((raw[max(k - 2, 0):k + 3] >= self.sat).sum()))
            res['_shift'] = s
            res['_rms'] = float(np.sqrt(np.mean((d[m] - gsum(yy[m], *p)) ** 2)))
            out[L['name']] = res
        self.fits = out
        return out

    def molar(self, lane, bid):
        L = self.blen[bid]
        if not L:
            sys.exit(f'band {bid} has no length; cannot convert to molar')
        return self.fits[lane][bid]['area'] / L

    def metrics(self):
        rows = []
        for mt in self.cfg.get('metrics', []):
            tags = set(mt.get('lanes_with_tag', []))
            for L in self.lanes:
                n = L['name']
                if n not in self.fits or (tags and not tags & set(L.get('tags', []))):
                    continue
                num = sum(self.molar(n, b) for b in mt['numerator'])
                den = sum(self.molar(n, b) for b in mt['denominator'])
                cn = [b for b in mt['numerator'] if self.fits[n][b]['clip']]
                cd = [b for b in mt['denominator'] if self.fits[n][b]['clip'] and b not in mt['numerator']]
                bound = '~' if cn and cd else '>=' if cn else '<=' if cd else ''
                rows.append(dict(lane=n, metric=mt['name'], value=num / den if den > 0 else float('nan'),
                                 bound=bound, clipped=';'.join(cn + cd)))
        self.metric_rows = rows
        return rows

    def amounts(self):
        am = self.cfg.get('amounts')
        if not am:
            self.amount_rows = []
            return []
        rows = []
        for L in self.lanes:
            n = L['name']
            if n not in self.fits or 'bead_ul' not in L:
                continue
            if am.get('lanes_with_tag') and not set(am['lanes_with_tag']) & set(L.get('tags', [])):
                continue
            for b in am['bands']:
                rows.append(dict(lane=n, band=b, per_unit=self.molar(n, b) / L['bead_ul'], clip=self.fits[n][b]['clip']))
        ref = am.get('reference')
        if ref:
            r = [x['per_unit'] for x in rows if x['lane'] == ref['lane'] and x['band'] == ref['band']][0]
            for x in rows:
                x['relative'] = x['per_unit'] / r
        self.amount_rows = rows
        return rows

    def standards(self):
        """Absolute pmol from standard lanes (role: standard, standard_pmol, standard_band)."""
        st = [L for L in self.lanes if L.get('role') == 'standard' and L['name'] in self.fits]
        if not st:
            return None
        pts = []
        for L in st:
            b = L['standard_band']
            pts.append(dict(lane=L['name'], pmol=L['standard_pmol'], area=self.fits[L['name']][b]['area'],
                            clip=self.fits[L['name']][b]['clip']))
        usable = [p for p in pts if not p['clip']]
        out = dict(points=pts, usable=len(usable))
        if len(usable) >= 2:
            k, c = np.polyfit([p['pmol'] for p in usable], [p['area'] for p in usable], 1)
            out.update(area_per_pmol=float(k), intercept=float(c))
        elif len(usable) == 1:
            out.update(area_per_pmol=usable[0]['area'] / usable[0]['pmol'], intercept=0.0,
                       note='one unclipped standard: line through the origin; add more standards')
        elif pts:
            # saturated standards only give an upper bound for unclipped sample bands
            # each clipped standard underestimates area/pmol; the largest ratio is the tightest lower bound
            k = max(p['area'] / p['pmol'] for p in pts)
            out['area_per_pmol_lower_bound'] = k
            out['note'] = ('standards clipped: true area/pmol >= this value, so sample pmol = area / this value '
                           'is an upper bound (valid for unclipped sample bands only)')
        return out


# ---------------------------------------------------------------- plots
def _style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8.5, 'font.family': 'DejaVu Sans',
                         'axes.spines.top': False, 'axes.spines.right': False})
    return plt


def plot_gel(g, path):
    pc = g.cfg.get('plots', {}).get('gel')
    if not pc:
        return
    plt = _style()
    x0, x1, y0, y1 = pc['crop']
    fig, ax = plt.subplots(figsize=pc.get('size', [5.2, 3.6]), dpi=220)
    sub = g.a[y0:y1, x0:x1]
    vmax = pc.get('vmax', np.percentile(sub, 99.5))
    ax.imshow(np.clip(sub / vmax, 0, 1) ** pc.get('gamma', 0.7), cmap='gray', aspect='auto', extent=[x0, x1, y1, y0])
    for L in g.lanes:
        ax.text(L['x'], y0 - 4, L.get('label', L['name']), ha='center', va='bottom', fontsize=7.5)
    for grp in pc.get('groups', []):
        xs = [L['x'] for L in g.lanes if L['name'] in grp['lanes']]
        ax.plot([min(xs), max(xs)], [y0 - 17] * 2, color='#333', lw=0.8, clip_on=False)
        ax.text(np.mean([min(xs), max(xs)]), y0 - 19, grp['label'], ha='center', va='bottom', fontsize=7.5)
    for an in pc.get('annotate', []):
        ax.annotate(an['text'], (x1, an['y']), (x1 + 8, an['y']), fontsize=7.5, va='center', annotation_clip=False,
                    arrowprops=dict(arrowstyle='-', lw=0.6, color='#444'))
    ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]
    if pc.get('caption'):
        fig.text(0.02, 0.01, pc['caption'], fontsize=6.5)
    fig.savefig(path, bbox_inches='tight', facecolor='white'); plt.close(fig)


def _xpos(order, gaps):
    """x positions with a gap between groups. order: list of lane names; gaps: list of group sizes."""
    xs, x = [], 0.0
    for gsize in gaps:
        for _ in range(gsize):
            xs.append(x); x += 1
        x += 0.6
    return np.array(xs[:len(order)])


def plot_metrics(g, path):
    pc = g.cfg.get('plots', {}).get('metrics')
    if not pc:
        return
    plt = _style()
    lanes, mets = pc['lanes'], pc['metrics']
    xs = _xpos(lanes, pc.get('group_sizes', [len(lanes)]))
    val = {(r['lane'], r['metric']): r for r in g.metric_rows}
    fig, ax = plt.subplots(figsize=pc.get('size', [6.4, 3.3]), dpi=220)
    w = 0.8 / len(mets)
    for i, m in enumerate(mets):
        for x, n in zip(xs, lanes):
            r = val.get((n, m['name']))
            if not r:
                continue
            v = r['value'] * 100
            xx = x + (i - (len(mets) - 1) / 2) * w
            ax.bar(xx, v, w * 0.92, color=m['color'], label=m['label'] if x == xs[0] else None)
            ax.text(xx, v + 1, r['bound'].replace('>=', '≥').replace('<=', '≤') + f'{v:.0f}', ha='center', va='bottom', fontsize=6.5)
    if pc.get('target'):
        ax.axhline(pc['target'], ls='--', lw=0.9, color='#666')
    ax.set_xticks(xs); ax.set_xticklabels([pc.get('tick_labels', {}).get(n, n) for n in lanes])
    for lab in pc.get('group_labels', []):
        idx = [lanes.index(n) for n in lab['lanes']]
        ax.text(xs[idx].mean(), -0.15, lab['text'], transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=8.5, fontweight='bold')
    ax.set_ylim(0, 100); ax.set_ylabel(pc.get('ylabel', 'Efficiency (%, molar)'))
    ax.grid(axis='y', lw=0.4, color='#e6e6e6'); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=7.5, loc='upper center', bbox_to_anchor=(0.5, 1.16), ncol=len(mets))
    if pc.get('caption'):
        fig.text(0.01, -0.17, pc['caption'], fontsize=6.5)
    fig.savefig(path, bbox_inches='tight', facecolor='white'); plt.close(fig)


def plot_amounts(g, path):
    pc = g.cfg.get('plots', {}).get('amounts')
    if not pc or not g.amount_rows:
        return
    plt = _style()
    lanes = pc['lanes']
    xs = _xpos(lanes, pc.get('group_sizes', [len(lanes)]))
    key = 'relative' if 'relative' in g.amount_rows[0] else 'per_unit'
    val = {(r['lane'], r['band']): r for r in g.amount_rows}
    fig, ax = plt.subplots(figsize=pc.get('size', [6.0, 3.4]), dpi=220)
    for x, n in zip(xs, lanes):
        bot, anyclip = 0, False
        for i, s in enumerate(pc['stack']):
            r = val[(n, s['band'])]; v = r[key]; anyclip |= bool(r['clip'])
            ax.bar(x, v, 0.7, bottom=bot, color=s['color'], edgecolor='white', lw=0.8, label=s['label'] if x == xs[0] else None)
            if i == 0:
                ax.text(x, v / 2, ('≥' if r['clip'] else '') + f'{v:.1f}', ha='center', va='center', color='white', fontsize=7.5, fontweight='bold')
            bot += v
        ax.text(x, bot * 1.01, ('≥' if anyclip else '') + f'{bot:.1f}', ha='center', va='bottom', fontsize=7.5)
    ax.set_xticks(xs); ax.set_xticklabels([pc.get('tick_labels', {}).get(n, n) for n in lanes])
    for lab in pc.get('group_labels', []):
        idx = [lanes.index(n) for n in lab['lanes']]
        ax.text(xs[idx].mean(), -0.15, lab['text'], transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=8.5, fontweight='bold')
    ax.set_ylabel(pc.get('ylabel', 'Relative amount'))
    ax.grid(axis='y', lw=0.4, color='#e6e6e6'); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=7.5, loc='upper left')
    if pc.get('caption'):
        fig.text(0.01, -0.16, pc['caption'], fontsize=6.5)
    fig.savefig(path, bbox_inches='tight', facecolor='white'); plt.close(fig)


# ---------------------------------------------------------------- main
def write_csv(path, rows):
    if not rows:
        return
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def run(cfg_path, image=None, outdir=None, compare=None, plots=True):
    cfg = yaml.safe_load(open(cfg_path))
    base = os.path.dirname(os.path.abspath(cfg_path))
    img = image or os.path.join(base, cfg['image'])
    outdir = outdir or os.path.join(base, cfg.get('outdir', 'gel_quant_out'))
    os.makedirs(outdir, exist_ok=True)
    g = Gel(cfg, img)
    g.fit(); g.metrics(); g.amounts(); st = g.standards()
    bf = [dict(lane=n, band=b, **v) for n, r in g.fits.items() for b, v in r.items() if not b.startswith('_')]
    write_csv(os.path.join(outdir, 'band_fits.csv'), bf)
    write_csv(os.path.join(outdir, 'metrics.csv'), g.metric_rows)
    write_csv(os.path.join(outdir, 'amounts.csv'), g.amount_rows)
    if st and 'area_per_pmol' in st:
        for r in bf:
            r['pmol'] = (r['area'] - st['intercept']) / st['area_per_pmol']
    summary = dict(image=img, saturation=g.sat, lanes={L['name']: L['x'] for L in g.lanes},
                   fit_rms={n: r['_rms'] for n, r in g.fits.items()},
                   lane_shift={n: r['_shift'] for n, r in g.fits.items()}, standards=st)
    if compare:
        g2 = Gel(yaml.safe_load(open(cfg_path)), compare); g2.fit(); g2.metrics()
        v2 = {(r['lane'], r['metric']): r for r in g2.metric_rows}
        cmp_rows = [dict(lane=r['lane'], metric=r['metric'], value=r['value'], value_compare=v2[(r['lane'], r['metric'])]['value'],
                         clipped=r['clipped'], clipped_compare=v2[(r['lane'], r['metric'])]['clipped']) for r in g.metric_rows]
        write_csv(os.path.join(outdir, 'compare.csv'), cmp_rows)
        un = [abs(c['value'] - c['value_compare']) for c in cmp_rows if not c['clipped'] and not c['clipped_compare']]
        summary['compare_max_abs_diff_unclipped'] = max(un) if un else None
    json.dump(summary, open(os.path.join(outdir, 'summary.json'), 'w'), indent=1)
    if plots:
        plot_gel(g, os.path.join(outdir, 'gel_annotated.png'))
        plot_metrics(g, os.path.join(outdir, 'metrics.png'))
        plot_amounts(g, os.path.join(outdir, 'amounts.png'))
    return g, summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('config'); ap.add_argument('--image'); ap.add_argument('--outdir'); ap.add_argument('--compare')
    ap.add_argument('--no-plots', action='store_true')
    a = ap.parse_args()
    g, summ = run(a.config, a.image, a.outdir, a.compare, not a.no_plots)
    for r in g.metric_rows:
        print(f"{r['lane']:12s} {r['metric']:14s} {r['bound']:>2s}{r['value']*100:5.1f}%  {r['clipped']}")
    if summ.get('standards'):
        print('standards:', json.dumps(summ['standards']))


if __name__ == '__main__':
    main()
