# TapeStation lessons

- **Encoding:** the CSV export is cp1252 (the `pg/µl` header has a micro sign). Reading as UTF-8 fails or
  mangles headers. Header spellings vary between assays; match by prefix.
- **No sample names in the export.** `Sample Description` is usually empty and `sampleTable.csv` holds
  only lot numbers. The sample sheet is the only lane record; keep it with the run.
- **Export folder vs run:** a folder named like the run stem (`YYYY-MM-DD - HH-MM-SS-HSD1000`) is the
  per-well export of an existing run, not a new run. Keep the export files with the run's other files.
- **PDF layout** (TapeStation Analysis 5.x): page 1 gel image and well order, page 2 ladder, then one
  page per well with peak and region tables; last page lists wells. Prefer the CSVs over reading
  numbers off the PDF.
- **Region window:** use the bounds the instrument applied (300-600 is common for short-read
  libraries; some runs used 200-600). Match exactly; nested regions (e.g. 400-532) are separate rows.
- **Thresholds that caught real problems:** region < 50% of lane (most of the lane is not library);
  >= 5% of area below 200 bp (another 0.8x bead cleanup needed).
- **Low-concentration lanes** (region molarity in the tens of pM on HS D1000) are near the floor;
  their sizes are unreliable and the main peak can be noise.
- **Lane boundaries on the instrument gel image:** a column with no dark pixel is a gap between lanes;
  use that instead of guessing when cropping the PDF gel.
- **Copying runs to shared storage:** verify by checksum on both sides, not by a transfer summary
  (rsync's "N files" counts the directory entry). Use `LC_ALL=C sort` on both sides when comparing lists.
