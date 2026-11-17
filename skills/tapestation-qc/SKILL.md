---
name: tapestation-qc
description: Analyze an Agilent TapeStation run (D1000/HS D1000/D5000, per-well CSV export) - per-sample peak and region size, molarity, % in region, primer/adapter-dimer fraction, warnings, electropherograms and a labelled gel. Use for any TapeStation, library QC, fragment size or Bioanalyzer-style sizing request.
---

# TapeStation QC

Script: `scripts/tapestation_qc.py`. Example run (synthetic, same file format as the instrument):
`examples/synthetic_run/`.

## Inputs

1. **The per-well CSV export** from TapeStation Analysis (File > Export): `*_compactPeakTable.csv`,
   `*_compactRegionTable.csv`, `*_Electropherogram.csv`, `*_runProperties.csv`, plus the PDF report.
   The instrument folder alone (`.HSD1000` + `.pdf`) holds no tables; ask for the export.
2. **A sample sheet** (`sample_sheet.csv`: `Well,Sample_Name[,Group,...]`). The export has no sample
   names. Draft one with `--draft-sheet`, fill it from the run record, and have the person confirm it
   before any number is attributed to a sample. A wrong lane silently moves numbers between samples.
3. **The region window** to report (default 300-600 bp). Use the window the instrument actually applied;
   the script matches it exactly and warns if absent.

## Run

```
python3 scripts/tapestation_qc.py RUN_DIR --draft-sheet          # once, then fill the sheet
python3 scripts/tapestation_qc.py RUN_DIR [--region 300 600] [--below 200]
```

Outputs in `RUN_DIR/tapestation_qc/`: `qc_table.csv`, `qc.json`, `electropherograms.png`,
`gel_labeled.png` (pseudo-gel from the traces, sized by the ladder, lanes named from the sheet,
brackets per `Group`).

## What is reported per sample

| Field | Source |
| --- | --- |
| Main peak bp, % area, pg/µl | largest non-marker peak inside the sizing range (peak table) |
| Other peaks | every other non-marker peak, bp (% area) |
| Region average bp, pg/µl, molarity (nM), % of total | region table, exact window |
| % area below N bp | sum of peak % areas below N (default 200 bp): primer/adapter dimer |

Mass from molarity, when needed: ng/µl = nM x 660 x bp / 1e6 (660 g/mol per bp of dsDNA).

## QC, required every time

- **Warnings** printed by the script: region missing, region < 50% of lane, >= 5% of area below the
  cut-off, peaks outside the sizing range, no sample peaks. Relay them; do not drop them.
- **Size check:** `qc.json` -> `size_check` compares the main peak read from the trace (ladder-sized)
  with the instrument table. Strong lanes should agree within ~3%. A large gap in a strong lane means
  wrong lane mapping or a bad ladder.
- **Lane mapping sanity check:** if another concentration measure exists (Qubit, qPCR), rank both lists
  and confirm the order agrees. Report the mapping as inferred until confirmed.
- **TapeStation molarity is not qPCR molarity.** Do not put it in a field meant for size-corrected qPCR
  (e.g. KAPA) values. Report both separately.

## Report

Figures: drawn by the script with `_shared/labstyle.py`; any extra figure follows `_shared/FIGURES.md`.
Captions are written to `<figure>.caption.txt`; paste them under the image, not into it.

- Table per sample: well, sample, main peak bp, region bp, nM, % in region, % below cut-off.
- The labelled gel and electropherograms as separate figures; no interpretation written on figures.
- Label the region column with its bounds (e.g. "300-600 bp").
- Follow `../../_shared/CONVENTIONS.md`.

See `references/lessons.md`.
