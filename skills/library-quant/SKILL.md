---
name: library-quant
description: qPCR quantification of sequencing libraries against DNA standards (KAPA Library Quantification or similar) - standard curve, size correction, dilution agreement, replicate and range checks, comparison with fluorometric (Qubit) concentration. Use for any library quant, KAPA, or "how much library do I load" request.
---

# Library quantification by qPCR

Script: `scripts/library_quant.py`. Example (synthetic, CFX export format): `examples/synthetic/`.
For absolute copy numbers of a specific target (not library molarity), use `qpcr-quant`.

## Inputs

1. **Instrument export:** Bio-Rad CFX Maestro CSV folder ("Quantification Cq Results.csv"; optional
   amplification and melt peak CSVs) or a QuantStudio .xls/.xlsx Results sheet.
2. **Plate layout:** standard wells with pM, NTC wells, and per library the wells of each dilution.
3. **Standard size** (KAPA DNA standards: 452 bp) and **each library's average size** from TapeStation or
   Bioanalyzer (the region average over the library range, not the main peak if the shape is broad).
4. Optional **fluorometric concentration** (ng/µl, Qubit dsDNA HS) per library for the comparison.

## Run

```
python3 scripts/library_quant.py config.yaml [--export-dir DIR]
```

Outputs: `libraries.csv`, `standards.csv`, `summary.json`, `library_quant.png` (standard curve with
libraries placed on it, amplification curves, concentrations vs fluorometric).

## Calculations

- Standard curve: Cq vs log10 pM in the well; efficiency = 10^(-1/slope) - 1.
- Undiluted pM = pM in well x dilution fold, per dilution (mean of replicates); dilutions outside the
  standard Cq range are not used.
- Size-corrected nM = mean undiluted pM x 452 / library bp / 1000.
- Fluorometric nM = ng/µl x 1e6 / (660 x library bp).

## QC, required every time

- Efficiency 90-110%, R² >= 0.99.
- Replicate Cq SD <= 0.2.
- The two dilutions differ by about |slope| x log10(fold ratio) cycles (3.3 for 10-fold); a smaller gap
  in the less diluted wells means inhibition or carryover.
- Library Cq inside the standard range; NTC well after the lowest standard (or no Cq).
- One melt peak per well; a second peak means adapter or primer dimer that qPCR counts as library.
- qPCR vs fluorometric ratio: libraries with low ratios contain material that does not amplify with the
  P5/P7 primers (incomplete adapters); report it, it changes loading.

## Report

Table: library, size bp, nM (size-corrected), dilution agreement, fluorometric nM and ratio, warnings.
Follow `../../_shared/CONVENTIONS.md`. TapeStation molarity and qPCR molarity are reported separately.
