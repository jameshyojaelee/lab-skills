# Countess lessons

- The CSV `Cube N cells counted` equals the number of objects with `IsFluorescent` true in Cloud.json;
  a mismatch means the JSON belongs to a different read.
- Gates are set per read on the instrument. The lowest brightness among positive objects is the
  effective gate; compare it across reads before comparing % positive over time.
- Clumps (`ClumpIndex` not null) inflate aggregation and can carry fluorescence from a neighbour;
  report % positive among non-clumped objects next to the instrument value.
- Brightfield-based counts without trypan blue give no live/dead split.
- Reads of the same sample on the same day are replicates; order them by Count ID.
- For low fluorescent fractions, pool counted cells across replicate reads and give a Wilson interval.
