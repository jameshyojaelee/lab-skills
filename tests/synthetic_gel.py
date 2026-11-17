"""Render a synthetic denaturing-gel image with known band amounts, for testing gel_quant.py.

Species (single-stranded, nt): a three-part sequential ligation on beads.
  top strand:    start 50 -> step-1 product 90 -> full length 160, plus free incoming oligo 70
  bottom strand: start 40 -> step-1 product 80 -> full length 120
Band signal = molecules x length (stain scales with length). 12-bit image, saturation 4095.
"""
import numpy as np
from PIL import Image

H, W = 420, 560
SAT = 4095
# size calibration y = A + B log10(nt)
A, B = 520.0, -180.0
def y_of(nt): return A + B * np.log10(nt)

LEN = {'FL': 160, 'FLb': 120, 'S1': 90, 'S1b': 80, 'FREE': 70, 'S0': 50, 'S0b': 40}

# lanes: name, x, molecules per species (arbitrary units), smile shift px
LANES = [
    ('Marker', 40, None, 0),
    ('A-step1', 90,  dict(S1=6, S1b=6, S0=2, S0b=2), 1),
    ('A-final', 140, dict(FL=4.0, FLb=4.0, S1=1.0, S1b=1.0, S0=2.0, S0b=2.0, FREE=1.0), 2),
    ('B-step1', 190, dict(S1=12, S1b=12, S0=6, S0b=6), 2),
    ('B-final', 240, dict(FL=8.0, FLb=8.0, S1=2.0, S1b=2.0, S0=4.0, S0b=4.0, FREE=2.0), 3),
    ('C-final', 290, dict(FL=6.0, FLb=6.0, S1=6.0, S1b=6.0, S0=3.0, S0b=3.0, FREE=1.0), 3),
    ('Std-low', 360, dict(FL=5.0), 2),
    ('Std-high', 410, dict(FL=40.0), 1),   # deliberately saturates
]
STD_PMOL = {'Std-low': 0.01, 'Std-high': 0.08}   # 5 units = 0.01 pmol


def render(path=None, seed=0, amp=6.0):
    rng = np.random.default_rng(seed)
    yy = np.arange(H)[:, None].astype(float)
    xx = np.arange(W)[None, :].astype(float)
    img = 60 + 0.05 * yy + np.zeros((H, W))                 # slow background gradient
    for name, x, sp, shift in LANES:
        lane = np.exp(-0.5 * ((xx - x) / 9.0) ** 4)          # flat-topped lane profile
        img += 25 * lane * np.exp(-yy / 300)                 # smear
        if sp is None:  # marker lane: a few sharp bands
            for nt in (40, 50, 70, 90, 120, 160):
                img += 300 * lane * np.exp(-0.5 * ((yy - y_of(nt)) / 2.0) ** 2)
            continue
        for k, mol in sp.items():
            sig = 2.2 + 0.004 * LEN[k]
            img += amp * mol * LEN[k] / (sig * np.sqrt(2 * np.pi)) * 2.0 * lane * \
                np.exp(-0.5 * ((yy - y_of(LEN[k]) - shift) / sig) ** 2)
    img += rng.normal(0, 4, img.shape)
    img = np.clip(img, 0, SAT).astype(np.uint16)
    if path:
        Image.fromarray(img).save(path)
    return img


def truth(lane):
    sp = dict(LANES[[l[0] for l in LANES].index(lane)][2])
    g = lambda *ks: sum(sp.get(k, 0) for k in ks)
    if 'step1' in lane:
        return dict(step1=g('S1') / g('S1', 'S0'))
    return dict(step1=g('FL', 'S1') / g('FL', 'S1', 'S0'), step2=g('FL') / g('FL', 'S1'),
                step2_bottom=g('FLb') / g('FLb', 'S1b'), overall=g('FL') / g('FL', 'S1', 'S0'))


if __name__ == '__main__':
    import sys
    render(sys.argv[1] if len(sys.argv) > 1 else 'synthetic_gel.tif')
    for n, *_ in LANES[1:6]:
        print(n, {k: round(v, 3) for k, v in truth(n).items()})
