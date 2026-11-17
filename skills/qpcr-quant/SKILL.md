---
name: qpcr-quant
description: Absolute qPCR quantification from a QuantStudio export - standard curve, Ct to copies, melt QC, single-bead (digital-style) wells with Poisson check and background subtraction, and bulk release per bead. Use for any qPCR copy-number, copies-per-bead, or standard-curve analysis.
---

# qPCR absolute quantification

Turns a QuantStudio export (.xls/.xlsx with Results, Melt Curve Result, Amplification Data sheets)
plus a YAML plate layout into copies per well, copies per bead (single-bead and bulk), Poisson
checks, QC tables and two figures. Script: `scripts/qpcr_quant.py`. Example layout: `examples/single-bead-bulk/`.

## Inputs to collect before running

1. **Export file** from QuantStudio (Results sheet must have a `Well` header row and a `CT` column;
   set `ct_column` if it differs).
2. **Plate map:** which wells are standards (with copies), samples, bare-bead and NTC controls.
3. **Standards:** stock concentration, dilution to the top standard, template volume, fold steps.
4. **Single-bead blocks:** wells and planned beads per well (lambda).
5. **Bulk releases:** eluate volume, template volume, dilution of each well, and the number of
   beads released (hemocytometer count per µl x dilution x volume taken).
6. **Expected product Tm** from the standard wells.

## Steps

1. Copy `configs/template.yaml` or `examples/single-bead-bulk/config.yaml` and fill in the layout.
2. Run `python3 scripts/qpcr_quant.py config.yaml --export run.xls`.
3. Check QC (below), then report.

## Calculations

- **Standard copies per well** = molar x template volume x Avogadro. Fit Ct = slope x log10(copies)
  + intercept on all standard wells. Efficiency = 10^(-1/slope) - 1.
- **Copies in a well** = 10^((Ct - intercept) / slope). Only inside the standard range.
- **Single-bead wells:** bead vs no-bead split at a Ct threshold (given, or the largest gap in the
  sorted Ct values, `threshold: auto`). Background = median copies of no-bead wells (free oligo in
  the suspension). Copies per bead = (well copies - background) / beads in the well; beads per well
  = round((copies - background) / one-bead unit), unit from the lowest-lambda block.
  Fitted lambda = -ln(fraction of empty wells); compare observed 0/1/2/>=3 counts with Poisson.
- **Bulk:** copies per bead = well copies / template µl x dilution x eluate µl / beads released.
  Two dilutions of the same eluate should differ by log10(fold ratio) x |slope| cycles
  (`bulk_dilution_check` in summary.json).

## QC, required every time

- Standard curve: efficiency 90-110%, R^2 >= 0.99, residuals without a trend.
- Melt: every sample and standard well one peak at the product Tm (+-1.5 °C). Bare beads and NTC
  only primer-dimer (low Tm) or no amplification. Wells with master mix but no template are not NTCs.
- Single-bead: a clear Ct gap between bead and no-bead wells; fitted lambda near planned. If
  bead-positive wells drop column by column in dispensing order, beads settled in the source.
- Background: no-bead wells that amplify mean free product in the suspension. Check it scales with
  the volume of suspension per well. Report it; it matters for printing (free oligos land anywhere).
- Bulk: dilution check within ~0.5 cycle; replicate spread.
- Every copy number must fall inside the standard range.

## Report

Figures: drawn by the script with `_shared/labstyle.py`; any extra figure follows `_shared/FIGURES.md`.
Captions are written to `<figure>.caption.txt`; paste them under the image, not into it.

- Separate a QC figure (standard curve, amplification, melt, plate map) from the results figure
  (Ct per well by group, Poisson, copies per bead).
- Give single-bead median and spread (CV) with n, and bulk per bead with range across dilutions.
- Compare with an independent measure when one exists (e.g. gel band ratio between loadings).
- Label claims: verified (computed), inferred (e.g. settling, source of free oligo), estimate
  (anything relying on bead counts from another experiment).

See `references/method.md` and `references/lessons.md`.
