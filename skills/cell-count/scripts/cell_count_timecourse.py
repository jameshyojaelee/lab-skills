#!/usr/bin/env python3
"""Rebuild a longitudinal Countess plot from a private JSON config.

Config: input_globs (relative to config), arm_regex (one named group 'arm'),
arm_labels (token -> label), channel, start_time (ISO timestamp), title,
repeat_count_ids (optional), outdir (relative to config).
Each export remains a separate observation; no pooling of repeat reads.
"""
import argparse
import csv
from datetime import datetime
import glob
import hashlib
import json
from pathlib import Path
import re

import os
import sys

import numpy as np

import cell_count
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', '..', '_shared'))
import labstyle


def run(config_path, plots=True):
    config_path = Path(config_path).resolve()
    cfg = json.loads(config_path.read_text())
    base = config_path.parent
    out = base / cfg.get('outdir', 'output')
    paths = sorted({p for pattern in cfg['input_globs']
                    for p in glob.glob(str(base / pattern))})
    recs, _, warnings = cell_count.run(paths, outdir=str(out / 'per_read'), plots=False)
    rx = re.compile(cfg['arm_regex'])
    start = datetime.fromisoformat(cfg['start_time'])
    rows = []
    ids = set()
    for r in recs:
        if r['count_id'] in ids:
            raise ValueError('Duplicate Count ID: ' + r['count_id'])
        ids.add(r['count_id'])
        match = rx.search(r['sample'])
        if not match:
            raise ValueError('Unmapped sample: ' + r['sample'])
        arm = cfg['arm_labels'][match.group('arm')]
        raw = next(csv.DictReader(open(r['file'], encoding='utf-8-sig')))
        cube = next((i for i in (1, 2) if r.get(f'cube{i}_name') == cfg['channel']), None)
        if cube is None:
            raise ValueError('Requested channel absent: ' + r['file'])
        n, k = int(r['total_counted']), int(r[f'cube{cube}_counted'])
        cloud = json.loads(Path(r['file']).with_name(Path(r['file']).stem + '_Cloud.json').read_text())
        cells = next(c['Cells'] for c in cloud['CellTypeResults'] if c['CellType'] == 'Brightfield')
        cells = [c for c in cells if c.get('IsValid', True)]
        flags = [bool(c['IsFluorescent'][cube - 1]) for c in cells]
        if len(cells) != n or sum(flags) != k:
            raise ValueError('CSV / JSON count mismatch: ' + r['file'])
        free = [c for c in cells if c.get('ClumpIndex') is None]
        nf = len(free)
        kf = sum(bool(c['IsFluorescent'][cube - 1]) for c in free)
        positive_brightness = [c['FluorescentBrightness'][cube - 1] for c, f in zip(cells, flags) if f]
        lo, hi = cell_count.wilson(k, n)
        flo, fhi = cell_count.wilson(kf, nf)
        # Older Python accepts at most six fractional digits; exports contain seven.
        timestamp = datetime.fromisoformat(re.sub(r'(\.\d{6})\d+', r'\1', r['time']))
        row = dict(arm=arm, count_id=r['count_id'], timestamp=r['time'],
                   date=timestamp.date().isoformat(), hours_post_start=(timestamp - start).total_seconds() / 3600,
                   read_role='repeat' if r['count_id'] in cfg.get('repeat_count_ids', []) else 'primary',
                   channel=cfg['channel'], positive=k, objects=n, positive_pct=100*k/n,
                   ci_low=lo, ci_high=hi, nonclumped_positive=kf, nonclumped_objects=nf,
                   nonclumped_pct=100*kf/nf if nf else None, nonclumped_ci_low=flo, nonclumped_ci_high=fhi,
                   total_per_ml=r['total_per_ml'], viability_measured=r['viability_measured'],
                   aggregation_pct=r['aggregation_pct'],
                   excitation=raw.get(f'Cube {cube} light intensity'),
                   focus=raw.get('Focus value'),
                   lowest_positive_brightness=min(positive_brightness) if positive_brightness else None,
                   stain_corrected=r['stain_corrected'], predilution_corrected=r['predilution_corrected'],
                   source_csv=str(Path(r['file']).resolve()),
                   source_sha256=hashlib.sha256(Path(r['file']).read_bytes()).hexdigest())
        for i in (1, 2):
            row[f'cube{i}_name'] = r.get(f'cube{i}_name', '')
            row[f'cube{i}_counted'] = r.get(f'cube{i}_counted', '')
            row[f'cube{i}_excitation'] = raw.get(f'Cube {i} light intensity', '')
        rows.append(row)
    for arm in cfg['arm_labels'].values():
        settings = {r['excitation'] for r in rows if r['arm'] == arm}
        if len(settings) > 1:
            warnings.append(f'{arm}: {cfg["channel"]} excitation differs across reads: {sorted(settings)}')
    with (out / 'timecourse.csv').open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    (out / 'provenance.json').write_text(json.dumps(dict(config=cfg, warnings=warnings,
        source_count=len(rows), method='Instrument calls recounted from Cloud.json; Wilson 95% intervals. '
        'Intervals describe counted-object sampling only; aggregation, gate uncertainty and biological '
        'variation are not included. Lowest positive brightness is descriptive, not a calibrated gate.'), indent=2))
    if plots:
        plot(rows, cfg, out)
    return rows, warnings


def plot(rows, cfg, out):
    """One panel per arm: % positive over days, Wilson 95% CI. Rules: _shared/FIGURES.md."""
    arms = list(dict.fromkeys(r['arm'] for r in rows))
    kind = 'single' if len(arms) == 1 else 'double' if len(arms) == 2 else 'triple'
    fig, axes = labstyle.figure(kind, ncols=len(arms), sharey=True)
    axes = list(np.atleast_1d(axes))
    days = lambda r: r['hours_post_start'] / 24
    for k, (ax, arm) in enumerate(zip(axes, arms)):
        col = labstyle.COLORS[k % len(labstyle.COLORS)]
        rr = [r for r in rows if r['arm'] == arm]
        prim = sorted((r for r in rr if r['read_role'] == 'primary'), key=days)
        rep = [r for r in rr if r['read_role'] == 'repeat']
        err = lambda g, f, lo, hi: [[r[f] - r[lo] for r in g], [r[hi] - r[f] for r in g]]
        ax.errorbar([days(r) for r in prim], [r['positive_pct'] for r in prim],
                    yerr=err(prim, 'positive_pct', 'ci_low', 'ci_high'), fmt='-o', color=col, ms=4,
                    label='All objects')
        nc = [r for r in prim if r['nonclumped_pct'] is not None]
        ax.errorbar([days(r) + 0.08 for r in nc], [r['nonclumped_pct'] for r in nc],
                    yerr=err(nc, 'nonclumped_pct', 'nonclumped_ci_low', 'nonclumped_ci_high'), fmt='o',
                    mfc='white', color=col, ms=4, label='Non-clumped')
        if rep:
            ax.errorbar([days(r) + 0.16 for r in rep], [r['positive_pct'] for r in rep],
                        yerr=err(rep, 'positive_pct', 'ci_low', 'ci_high'), fmt='D', mfc='white',
                        color=labstyle.GREY, ms=3.5, label='Repeat read')
        ax.set_title(arm)
        ax.set_ylim(0, 100)
        d = [days(r) for r in rr]
        lo, hi = int(np.floor(min(d))), int(np.ceil(max(d)))
        ax.set_xlim(lo - 0.5, hi + 0.5)
        ax.set_xticks(range(lo, hi + 1))
        if k == 0:
            ax.set_ylabel(f"{cfg['channel']}+ objects (%)")
    fig.supxlabel(f"Days after {cfg.get('start_label', 'transduction')}", fontsize=7)
    handles = {}
    for ax in axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            handles.setdefault(l, h)
    axes[0].legend(handles.values(), handles.keys(), loc='upper left')
    caption = (f"{cfg['channel']}-positive objects / counted objects, per Countess read; error bars are Wilson 95% "
               "intervals on counted objects (sampling only; clumping and gate uncertainty excluded). Filled: all "
               "objects, joined across days; open circles: non-clumped objects; open diamonds: repeat reads, not "
               "pooled. Points offset slightly in x for visibility. " + cfg.get('caption_notes', '')).strip()
    labstyle.save(fig, str(out / 'timecourse'), caption=caption)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config')
    args = parser.parse_args()
    data, warnings = run(args.config)
    for r in data:
        print(f'{r["arm"]} {r["date"]} #{r["count_id"]}: {r["positive"]}/{r["objects"]} = {r["positive_pct"]:.2f}%')
    for warning in warnings:
        print('QC:', warning)
