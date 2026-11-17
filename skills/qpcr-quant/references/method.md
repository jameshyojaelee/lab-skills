# qPCR quantification: method and reasoning

## Standard curve

Copies per well for a standard = concentration (M) x template volume (L) x 6.022e23.
Example: a synthetic target strand at 10 µM diluted 1e6-fold to 10 pM, 2 µl per well
= 10e-12 x 2e-6 x 6.022e23 = 1.20e7 copies, then 10-fold steps to 120 copies.
Fit Ct against log10(copies). A typical good curve: Ct = -3.37 log10(N) + 33.8, efficiency 98%, R^2 0.9997.
A slope of -3.32 is 100% efficiency.

## Single-bead wells

At a planned lambda of 0.2 beads per well, most wells have no bead, so the Ct distribution is
bimodal: a tight no-bead cluster (free oligo in the liquid) and an earlier bead cluster. The
threshold goes in the gap. The no-bead cluster's median is the background added to every well,
so a bead's own copies = well copies - background.

Beads per well: (copies - background) divided by the one-bead unit, rounded, at least 1 for bead
wells. The unit comes from the low-lambda block, where nearly all bead wells hold one bead.

Poisson check: fraction of empty wells = e^(-lambda), so fitted lambda = -ln(empty fraction).
Compare observed counts of 0/1/2/>=3 beads with Poisson at the planned and fitted lambda. A good
match shows each early-Ct well is one bead and the count plus dilution are right.

## Bulk release

Release from a counted number of beads into a known eluate volume (e.g. NaOH then neutralise), dilute,
and run. Copies per bead = (well copies / template µl) x dilution fold x eluate µl / beads.
Run two dilutions 10-fold apart: their Ct difference should be |slope| (about 3.3); a smaller
difference points to inhibition in the less diluted well or pipetting error.

## Background, free oligo

No-bead wells that amplify well above NTC mean free target in the bead suspension. If the
background per well scales with the suspension volume per well, it comes from the suspension, not
the reagents. Report free copies per bead's worth of suspension next to copies on the bead.
