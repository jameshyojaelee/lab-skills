# Lessons from past gels

Practical points from multi-step ligation products run on denaturing PAGE
(Novex 6% TBE-urea, 1.0 mm, XCell SureLock, SYBR Gold, Bio-Rad Gel Doc / Image Lab).

## Imaging
- **Export the raw `.scn`, not only 8-bit TIFFs.** 8-bit exports clip at the display maximum, so bright
  bands and standards read too low even when the camera was not saturated. Two standards 2x apart
  giving equal areas is the sign.
- **TIFF exports carry no exposure time.** File names are save times. Record the exposure separately.
- **Standards must bracket the samples.** On a 6% TBE-urea mini gel with SYBR Gold, 0.125 pmol of a
  160-mer saturated while bead-released samples peaked far lower. Use a series around 0.01-0.1 pmol.
- Keep one exposure where no sample band and no standard is clipped.

## Running the gel
- Fill both chambers to the marks (about 1 L 1x TBE), pre-warm the buffer (~55 °C), pre-run 30 min at
  180 V and rinse wells before loading. An under-filled outer chamber gave skewed lanes and no bands.
- Load <= 10 µl per well in a 1.0 mm 15-well gel; overloaded wells and excess dye ruined a gel.
- Unheated dsDNA ladders are useless for ssDNA sizing and can spill into neighbouring wells. Use
  denatured single-oligo markers of known length (formamide, 95 °C).
- Leave an empty well after a heavily loaded marker lane.

## Analysis
- **Close top/bottom-strand doublets** (e.g. 90 nt top and 80 nt bottom) must be fitted as two peaks.
  It is the least reliable split on the gel.
- **Free incoming oligos stay on beads after washes** (seen at 10-25% of top-strand bands). Exclude them
  from yields and report them; free functional oligos (e.g. amine-modified) behave like product downstream.
- **Unassigned bands** can change a step yield if they are later identified. Report them separately.
- **The same step measured two ways can disagree** (step 1 from step-1-only lanes vs from final lanes).
  Run each step-1 sample next to its matching final sample to test it.
- **Sample amount per lane is noisy** when loading from bead suspensions (15-30% between tubes from
  one batch). Per-sample amount differences smaller than that are not real.
- Biotin, spacer and dU modifications make strands run slower than their length.

## Reporting
- Readers asked "what were the controls?" after a plot without a methods note. Always state the
  controls present and missing next to the numbers.
- Say which metric a percentage is. A last-step efficiency was read as overall yield.
- Efficiency bars per step and overall, plus relative amounts per sample, were clearer than
  100%-stacked composition plots.
