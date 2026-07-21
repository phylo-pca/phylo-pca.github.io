# PhyloPCA 0.0.0.9000

- `RunTNTParsimony(type = "landmark")` now accepts a list of landmark
  characters (with a matching `dim` vector) to search two or more shape
  characters jointly, e.g. a combined 2D+3D analysis -- TNT's own xread
  format allows any number of `&[landmark Nd]` blocks of differing
  dimensionality in one matrix.
- `scale_fig1_dists`'s `mp`-arm `SPR` column is populated (was `NA`):
  recomputed on Hamilton, guarded by `is.binary()` and a per-cell
  subprocess `timeout` (999/1000 cells; the remaining cell is a
  genuinely intractable pair under `TBRDist::USPRDist()`, left `NA`).
- Initial public release: simulation, alignment, tree-inference and scoring
  functions, and cached results, underlying Smith (2026), a reply to
  Raskin et al. (2026).
