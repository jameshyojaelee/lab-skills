"""Shared figure style for every skill. Rules: _shared/FIGURES.md.

    # in a skill script (works through symlinked skill folders):
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', '..', '_shared'))
    import labstyle
    labstyle.setup()
    fig, ax = labstyle.figure('single')     # or 'double', 'triple', 'full'; nrows=, ncols=
    ... plot with labstyle.COLORS / labstyle.STAGES ...
    labstyle.save(fig, 'out/efficiency', caption='Step 2 = 160/(160+90), molar, top strand.')

save() writes <stem>.pdf (vector master, editable text) and <stem>.png (preview for Notion/Slack),
writes the caption to <stem>.caption.txt, and runs check(). Problems are printed; with
strict=True (or LABSTYLE_STRICT=1) they raise.
"""
import os
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text

# ---------------------------------------------------------------- palette
# Stage / magnitude family (dark = furthest along / most). Same meaning in every figure.
STAGES = ['#0b3d91', '#2a78d6', '#8fb8ea', '#cfe0f6']
# Categorical, fixed order. Assign by entity, never by rank; a filter must not repaint survivors.
COLORS = ['#2a78d6', '#eb6834', '#1b9e77', '#7b5ea7', '#6b6b6b']
BLUE, ORANGE, GREEN, PURPLE, GREY = COLORS
INK, MUTED, RULE = '#222222', '#666666', '#dddddd'
TARGET = '#888888'           # reference / target lines (dashed)

# ---------------------------------------------------------------- sizes (inches)
SIZES = {
    'single': (3.5, 2.6),    # one panel
    'double': (7.0, 2.6),    # two panels in a row
    'triple': (7.0, 2.4),    # three panels in a row
    'full':   (7.0, 4.8),    # 2 x 2 or a tall panel
    'image':  (3.5, 3.0),    # gel / micrograph with labels
}
MAX_W, MAX_H = 7.2, 9.0
MIN_FONT = 6.0


def setup():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
        'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7,
        'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5,
        'axes.titleweight': 'normal', 'axes.titlelocation': 'left',
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'axes.edgecolor': INK, 'axes.labelcolor': INK, 'text.color': INK,
        'xtick.color': INK, 'ytick.color': INK,
        'axes.grid': False, 'grid.color': RULE, 'grid.linewidth': 0.4,
        'lines.linewidth': 1.0, 'lines.markersize': 4, 'errorbar.capsize': 2,
        'legend.frameon': False, 'legend.handlelength': 1.2, 'legend.borderaxespad': 0.2,
        'axes.prop_cycle': matplotlib.cycler(color=COLORS),
        'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03,
        'figure.facecolor': 'white', 'savefig.facecolor': 'white',
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',   # editable text in Illustrator
    })


def figure(kind='single', nrows=1, ncols=None, size=None, **kw):
    """Fixed-size figure with constrained layout. kind sets the size; 'double'/'triple' default to 2/3
    panels in a row unless ncols is given (ncols=1 gives one wide panel)."""
    setup()
    w, h = size or SIZES[kind]
    if ncols is None:
        ncols = {'double': 2, 'triple': 3}.get(kind, 1)
    return plt.subplots(nrows, ncols, figsize=(w, h), layout='constrained', **kw)


def panel_label(ax, letter):
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes, ha='right', va='bottom', fontsize=8, fontweight='bold')


def empty(ax, msg='not in export'):
    """Mark a panel whose data is missing instead of drawing blank axes."""
    ax.set_xticks([]); ax.set_yticks([])
    [sp.set_visible(False) for sp in ax.spines.values()]
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha='center', va='center', fontsize=6.5, color=MUTED)


def target_line(ax, y, label=None):
    ax.axhline(y, ls=(0, (3, 2)), lw=0.7, color=TARGET, zorder=0)
    if label:
        ax.annotate(label, (1, y), xycoords=('axes fraction', 'data'), xytext=(0, 2), textcoords='offset points',
                    ha='right', va='bottom', fontsize=6, color=MUTED)


def bar_labels(ax, bars, fmt='{:.0f}', prefix=None, inside=True):
    """Value labels on bars: inside the top when the bar is tall enough, else just above."""
    ymax = ax.get_ylim()[1]
    for i, b in enumerate(bars):
        v = b.get_height()
        if v != v:
            continue
        txt = (prefix[i] if prefix else '') + fmt.format(v)
        x = b.get_x() + b.get_width() / 2
        if inside and v > 0.12 * ymax:
            ax.text(x, v - 0.02 * ymax, txt, ha='center', va='top', fontsize=6, color='white', fontweight='bold')
        else:
            ax.text(x, v + 0.01 * ymax, txt, ha='center', va='bottom', fontsize=6, color=INK)


def group_axis(ax, positions, labels, groups, y=-0.16):
    """Tick labels per bar plus a bold group label underneath. groups: [(label, [positions])]."""
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='x', length=0)
    for lab, xs in groups:
        ax.text(sum(xs) / len(xs), y, lab, transform=ax.get_xaxis_transform(), ha='center', va='top',
                fontsize=7, fontweight='bold')


# ---------------------------------------------------------------- checks
def _texts(fig):
    """Visible, non-empty texts that are actually drawn (tick labels outside the view limits are skipped)."""
    ticklabels, drawn = set(), []
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                for lab in (tick.label1, tick.label2):
                    ticklabels.add(id(lab))
                    loc = tick.get_loc()
                    if (lab.get_visible() and lab.get_text().strip() and loc is not None
                            and lo - 1e-9 * abs(hi - lo) <= loc <= hi + 1e-9 * abs(hi - lo)):
                        drawn.append(lab)
    for t in fig.findobj(Text):
        if id(t) in ticklabels:
            continue
        if t.get_visible() and t.get_text().strip() and t.get_alpha() != 0:
            drawn.append(t)
    return drawn


def check(fig):
    """Return a list of rule violations (see FIGURES.md)."""
    issues = []
    w, h = fig.get_size_inches()
    if w > MAX_W + 0.01 or h > MAX_H + 0.01:
        issues.append(f'figure is {w:.1f} x {h:.1f} in; max {MAX_W} x {MAX_H}')
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    texts = _texts(fig)
    for t in texts:
        if t.get_fontsize() < MIN_FONT - 0.01:
            issues.append(f'font {t.get_fontsize():.1f} pt < {MIN_FONT} pt: "{t.get_text()[:30]}"')
    figtexts = [t for t in fig.texts if t.get_text().strip()]
    for t in figtexts:
        if '\n' in t.get_text() or len(t.get_text()) > 120:
            issues.append(f'caption-like text inside the figure ("{t.get_text()[:40]}..."); put it in the caption')
    if fig._suptitle is not None and fig._suptitle.get_text().strip() and len(fig.axes) == 1:
        issues.append('figure-wide title on a single panel; use the axis title or the caption')
    import math
    boxes = [(t, t.get_window_extent(r).expanded(0.98, 0.9)) for t in texts]
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            ti, tj = boxes[i][0], boxes[j][0]
            ri, rj = ti.get_rotation() % 180, tj.get_rotation() % 180
            if ri == rj and ri not in (0, 90):
                # parallel rotated labels (e.g. slanted tick labels): their axis-aligned boxes overlap even
                # when the text does not. They collide only if the anchors are closer, perpendicular to
                # the text direction, than the text height.
                (xi, yi), (xj, yj) = ti.get_window_extent(r).p0, tj.get_window_extent(r).p0
                th = math.radians(ri)
                perp = abs(-(xj - xi) * math.sin(th) + (yj - yi) * math.cos(th))
                if perp >= 0.9 * ti.get_fontsize() * fig.dpi / 72:
                    continue
            a, b = boxes[i][1], boxes[j][1]
            ov = max(0, min(a.x1, b.x1) - max(a.x0, b.x0)) * max(0, min(a.y1, b.y1) - max(a.y0, b.y0))
            if ov > 4:
                issues.append(f'text overlap: "{boxes[i][0].get_text()[:20]}" / "{boxes[j][0].get_text()[:20]}"')
    # text on top of data markers
    for ax in fig.axes:
        pts = []
        for ln in ax.get_lines():
            if not ln.get_visible():
                continue
            xy = ln.get_xydata()
            if len(xy) and (ln.get_marker() not in (None, '', 'None', ' ') or len(xy) > 20):
                # markers, and dense traces (vertices are close enough to stand in for the line)
                pts += list(ln.get_transform().transform(xy))
        for c in ax.collections:
            try:
                off = c.get_offsets()
                if len(off) and c.get_visible():
                    pts += list(c.get_offset_transform().transform(off))
            except Exception:
                pass
        for t in ax.texts:
            if not t.get_text().strip():
                continue
            bb = t.get_window_extent(r)
            pad = 2.0 * fig.dpi / 72          # about a marker radius
            for x, y in pts:
                if bb.x0 - pad < x < bb.x1 + pad and bb.y0 - pad < y < bb.y1 + pad:
                    issues.append(f'text on a data point: "{t.get_text()[:20]}"')
                    break
    n_legend_entries = sum(len(ax.get_legend().get_texts()) for ax in fig.axes if ax.get_legend())
    if n_legend_entries > 6:
        issues.append(f'{n_legend_entries} legend entries; split the figure or label directly')
    return sorted(set(issues))


def save(fig, stem, caption=None, png=True, strict=None, quiet=False):
    """Write <stem>.pdf (+ .png, + .caption.txt) after checking the figure against FIGURES.md."""
    strict = strict if strict is not None else os.environ.get('LABSTYLE_STRICT') == '1'
    issues = check(fig)
    if issues:
        msg = f'{os.path.basename(stem)}: ' + '; '.join(issues)
        if strict:
            plt.close(fig)
            raise ValueError('figure rules: ' + msg)
        if not quiet:
            print('FIGURE WARNING', msg)
    os.makedirs(os.path.dirname(os.path.abspath(stem)), exist_ok=True)
    fig.savefig(stem + '.pdf')
    if png:
        fig.savefig(stem + '.png', dpi=300)
    if caption:
        with open(stem + '.caption.txt', 'w') as fh:
            fh.write(textwrap.dedent(caption).strip() + '\n')
    plt.close(fig)
    return issues

