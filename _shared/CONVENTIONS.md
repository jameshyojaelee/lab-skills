# Reporting conventions (all skills)

1. **Define every number.** Each metric is reported with its formula or source column, in the caption
   or table header. "80% efficiency" alone is not a result.
2. **Controls.** List the controls that were run and the ones that were not, next to the numbers.
3. **Bounds and limits.** Values from saturated, clipped, out-of-range or low-count data are reported as
   bounds (>= / <=) or with confidence intervals, never as plain values.
4. **Claim labels.** Verified (computed from the data here), inferred (identity by size, mechanisms,
   lane mappings not yet confirmed), estimate (relies on numbers from another experiment).
5. **Figures.** Raw image (gel, trace) and summary plots are separate figures. No interpretation written
   on figures. One plain caption per figure.
6. **Raw over processed.** Analyze the instrument's raw export at full bit depth where it exists; say so
   when only a processed export was available.
7. **Notebook entries stay bare:** design, volumes, steps, recorded results, figures, one line on the
   analysis method. Interpretation, warnings and decisions go to the discussion or meeting page.
8. **Transcribe, don't fix.** Instrument values are recorded as exported. Missing values are marked
   `TO FILL`, not omitted.
9. **Plain language.** Short sentences, no em dashes, no filler. Call undiluted product "undiluted".
10. **Private details stay private.** Sample names, project codes, database IDs and channel IDs live in
    a private config outside this repo (see `config/README.md`).
