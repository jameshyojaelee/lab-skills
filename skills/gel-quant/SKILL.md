---
name: gel-quant
description: Quantify bands on a gel image (denaturing PAGE, TBE-urea, agarose) and report ligation or product efficiencies, relative amounts per sample, and absolute pmol from standards. Use for any gel quantification, band intensity, ligation yield, or "analyze this gel" request.
---

# Gel band quantification

Turns a gel image plus a YAML description of the gel into band areas, efficiency metrics,
relative amounts, absolute pmol (when standards allow), and three figures. Script:
`scripts/gel_quant.py`. Worked example with known answers: `examples/synthetic/` (image made by `tests/synthetic_gel.py`).

## When to use

- A gel image (Image Lab `.scn`, 16-bit TIFF, or 8-bit TIFF export) needs numbers.
- Someone asks for ligation efficiency, % full length, intermediates, or yield per sample.

## Inputs to collect before running

Ask for anything missing. Do not guess lane identities.

1. **Image file.** Prefer the raw `.scn` (full bit depth). An 8-bit TIFF export is clipped at the
   display maximum, so bright bands read too low. If only TIFFs exist, say so in the report.
2. **Lane map:** which sample is in which well, including empty wells and markers.
3. **Species per lane:** every strand that can appear, with length (nt) and strand (top/bottom).
   For ligation builds: each intermediate, the full-length product, and free incoming oligos.
4. **Amount loaded per lane** (e.g. µl bead stock), if amounts per bead are wanted.
5. **Standards:** lane, species and pmol of each standard.
6. **Controls on the gel** (no-ligase, no-incoming-oligo, size standards). List what was and was
   not run; this goes into the report.

## Steps

1. Copy `configs/template.yaml` (or an example config) next to the image and fill it in.
2. Find lane x positions and band y positions: open the image, or run once without `x` to use
   auto lane detection, then check `summary.json` and the annotated gel. Set each band's expected
   y in a reference lane, or give `calibration` points (nt, y) and band lengths.
3. Set the per-lane `anchor` (gel smile correction): the brightest point of a band present in
   that lane (`mode: max`) or the centroid of a doublet (`mode: centroid`).
4. Include every visible band in `bands`, even ones not used in metrics. Unfitted neighbours
   leak signal into fitted bands.
5. Run: `python3 scripts/gel_quant.py config.yaml [--compare other_exposure.tif]`.
6. Check QC (below) before reporting any number.
7. Report (below).

## QC, required every time

- **Clipping.** `band_fits.csv` column `clip` counts saturated pixels per band. Any clipped band
  makes metrics that use it bounds, not values (`>=` or `<=` in `metrics.csv`). Report them that way.
- **Standards.** Two standards with a 2x amount difference should give about 2x area. Equal areas
  mean saturation. Then pmol is only an upper bound (`summary.json` -> `standards`).
- **Exposure agreement.** With `--compare`, unclipped metrics should agree within a few points
  (`compare_max_abs_diff_unclipped`); a few points between two 8-bit exports is typical.
- **Noise floor.** Look at lanes where a band cannot exist (e.g. full length in a lane that never
  got the last oligo). Their fitted area is the noise floor; ignore shares below it.
- **Sample amount noise.** Tubes split from one batch should give the same total signal per µl.
  Their spread (often 15-30% for bead suspensions) is the error on any per-sample amount comparison.
- **Fit residual.** `summary.json` -> `fit_rms` per lane; a lane much worse than its neighbours
  needs a look (wrong anchor, missing band in the list).

## Metrics

Define each in the config as numerator bands / denominator bands, using molar amounts
(area / length). For a sequential ligation A -> A+B -> A+B+C on the top strand:

| Metric | Formula | Meaning |
| --- | --- | --- |
| Step 1 | (ABC + AB) / (ABC + AB + A) | Share of A that got B |
| Step 2 | ABC / (ABC + AB) | Of those with B, share that also got C |
| Overall | ABC / (ABC + AB + A) | Share of A that reached full length = step 1 x step 2 |

Each step is measured against the molecules that could take part in it.
Exclude free incoming oligos and unassigned bands from all three, and say so.

## Report

Figures: drawn by the script with `_shared/labstyle.py`; any extra figure follows `_shared/FIGURES.md`.
Captions are written to `<figure>.caption.txt`; paste them under the image, not into it.

- **Separate figures** (gel image, efficiency bars, relative amounts), per `../../_shared/FIGURES.md`.
  per sample. Do not combine them, and do not use 100%-stacked plots for the main result.
- **Every caption defines its metric** with the formula, and marks bounds (>=).
- **A short "how quantified" block** for readers outside the analysis: band signal = area under the
  peak (= band volume), area / length = molecules, a table of what each band is, the formulas,
  one worked example lane, the controls present and missing, and the limits.
- Label claims: verified (computed from data), inferred (band identity by size only, mechanisms),
  estimate (absolute amounts borrowed from another experiment).
- Lab notebook entries stay bare: results table, figures, one line on the method. Interpretation
  goes to the meeting page or chat.

See `references/method.md` for the reasoning behind each step and `references/lessons.md` for
pitfalls from past gels.
