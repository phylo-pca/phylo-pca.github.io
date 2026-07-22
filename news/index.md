# Changelog

## PhyloPCA 0.0.0.9000

- `scale_fig1_dists`‘s `mp`-arm `SPR` is now fully populated
  (1000/1000): the one cell that did not resolve within a 6-hour
  Hamilton budget under
  [`TBRDist::USPRDist()`](https://ms609.github.io/TBRDist/reference/TreeRearrangementDistances.html)
  (tree 176, nC 25 – a genuine NP-hard-worst-case instance, not a bug:
  RF = 20 between two markedly dissimilar 21-tip trees) is filled from
  [`TreeDist::SPRDist()`](https://ms609.github.io/TreeDist/reference/SPRDist.html)’s
  polynomial-time approximation instead, flagged `SPR_exact = FALSE`
  (matching the convention already used for the other arms’ occasional
  approximate fallbacks).
- `RunTNTParsimony(type = "landmark")` now accepts a list of landmark
  characters (with a matching `dim` vector) to search two or more shape
  characters jointly, e.g. a combined 2D+3D analysis – TNT’s own xread
  format allows any number of `&[landmark Nd]` blocks of differing
  dimensionality in one matrix.
- `scale_fig1_dists`’s `mp`-arm `SPR` column is populated (was `NA`):
  recomputed on Hamilton, guarded by
  [`is.binary()`](https://rdrr.io/pkg/ape/man/is.binary.tree.html) and a
  per-cell subprocess `timeout` (999/1000 cells; the remaining cell is a
  genuinely intractable pair under
  [`TBRDist::USPRDist()`](https://ms609.github.io/TBRDist/reference/TreeRearrangementDistances.html),
  left `NA`).
- Initial public release: simulation, alignment, tree-inference and
  scoring functions, and cached results, underlying Smith (2026), a
  reply to Raskin et al. (2026).
