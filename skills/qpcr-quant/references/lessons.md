# Lessons from past qPCR runs

Single-bead and bulk-release absolute qPCR on a QuantStudio 5 (384-well, SYBR, NEB Luna Universal),
Dynabeads M-270, C-Chip hemocytometer.

- **Polymerase:** use a non-proofreading mix (e.g. NEB Luna Universal). Proofreading enzymes (e.g. KAPA
  HiFi) have 3'->5' exonuclease that can degrade single-stranded templates and primers.
- **Beads settle in the source plate** during multi-pass dispensing: bead-positive wells per column
  fell steadily in dispensing order. Mix the source column before every pass, or refill halfway.
- **Hemocytometer:** C-Chip DHC-F01 (Fuchs-Rosenthal, 0.2 mm deep), large square = 0.2 µl, count at
  least 2 large squares. Undiluted M-270 stock at 1:100 is too dense to count; dilute further.
- **Free target in the bead suspension** is common (no-bead wells well above NTC, scaling with the
  suspension volume per well; seen at 25-30% of signal). Subtract it and report it.
- **Dilution check:** two bulk dilutions 10-fold apart should differ by |slope| (~3.3 cycles). One pair
  differed by 2.4, so that sample is a range. Run each dilution in duplicate.
- **Wells with master mix but no template** give low-Tm artifacts. Exclude them; they are not NTCs.
- **Gel and qPCR agreed** on relative amounts between samples (1.78x vs 1.75x), so a calibrated gel
  band can cross-check qPCR when standards are clean.
- **Expected numbers for a sanity check:** single-bead and bulk copies per bead should agree within
  ~1.5x; fitted lambda within ~40% of planned.
