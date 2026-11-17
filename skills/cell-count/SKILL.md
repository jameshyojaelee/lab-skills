---
name: cell-count
description: Analyze Invitrogen Countess 3 / 3 FL exports (CSV per read, optional Cloud.json per-cell data) - concentration, viability, size, fluorescence % with confidence intervals, replicate summary, aggregation and gate checks, and seeding volumes. Use for any cell count, viability, Countess or fluorescent-fraction request.
---

# Cell counts (Countess 3 / 3 FL)

Script: `scripts/cell_count.py`. Example (synthetic, instrument file format): `examples/synthetic/`.

## Inputs

1. **The Countess export folder(s).** Each read is a folder with `<name>.csv` (one row, header starts
   with `Count ID`), `<name>_Cloud.json` (per-object data), images and `.fcs`. Point the script at the
   folders or a parent folder.
2. **How reads group into replicates:** by default a trailing `_rep1`/`_1`/`-1` is stripped from the
   sample name; otherwise pass `--group-regex` with one capture group.
3. **Target cells to seed**, if a seeding volume is wanted (`--seed-cells`).

## Run

```
python3 scripts/cell_count.py COUNTESS_DIR [--group-regex '(.*)_r\d+'] [--seed-cells 2e5] [--outdir out]
```

Outputs: `reads.csv`, `groups.csv`, `summary.json`, `counts.png`, `fluorescence.png`, `brightness.png`.

## What is reported

- Per read: total / live / dead per ml, counted cells, viability, live size, aggregation, and per cube:
  fluorescent % with a Wilson 95% CI on the counted cells.
- From Cloud.json: positives recounted per object (must equal the CSV count), % positive among
  non-clumped objects, and the effective gate (lowest brightness the instrument called positive).
- Per group: mean and CV across reads, pooled fluorescent % with CI, seeding volume
  (µl = cells wanted / cells per ml x 1000; live cells if a viability stain was used, else total).

## QC, required every time

- **Viability needs a stain.** Brightfield-only counts report live = dead = 0; the script treats those as
  not measured. Never report 0% viability from such a read.
- **Counted cells < 100:** counting error is about 1/sqrt(n); report the CI, not just the %.
- **Aggregation > 20%:** concentration is underestimated and fluorescence % mixes clumps; compare with
  the non-clumped value.
- **Range:** accurate roughly between 1e5 and 1e7 cells/ml; dilute or concentrate outside it.
- **Gates:** compare the effective gate across reads; a moved gate changes % positive by itself.
- **Replicates:** two reads of the same tube should agree within counting error; a larger gap means
  poor mixing.
- Check `1:1 Stain corrected` and `Pre-Dilution corrected` flags before using concentrations.

## Report

Table per group (n reads, live/ml with CV, viability, fluorescent % with CI, seeding µl) and the
figures. Follow `../../_shared/CONVENTIONS.md`. See `references/lessons.md`.
