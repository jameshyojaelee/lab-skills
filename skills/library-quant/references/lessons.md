# Library quant lessons

- **CFX exports zero-pad wells in the Cq table** (`A01`) but not in the amplification table (`A1`);
  normalise before matching.
- **Standards singlicate** worked (R² > 0.999) when the six-point series was fresh; run duplicates when
  only one dilution of each library is used.
- **Low qPCR / fluorometric ratios** (around 0.1) flagged libraries whose fluorometric signal came from
  material without both adapters; size from TapeStation alone does not show this.
- **High ratios** (1.3-1.8) are common for well-formed libraries: fluorometric assays under-read when
  the size used is off or the library is partly single-stranded after bead cleanups.
- Use the library's average size over its region, from the same cleanup as the quantified material.
- Dilute in buffer with carrier (10 mM Tris pH 8, 0.05% Tween-20) to avoid losses at 1:10,000+.
