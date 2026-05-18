Changes in version 0.5.1 (2026-05-18)
+ Skip recreation of existing `lazygas/scan`, `lazygas/peakcall`, `lazygas/recalc`, and `lazygas/candidate` GDS nodes on reruns; update peakcall attributes when the peakcall node already exists.
+ Document the `limit_peakcall` argument in `callPeakBlock()`.
+ Fix `buildLazyGas()` when creating a GDS via `create_gds`: create the `annotation/format` folder before writing haplotype or dosage data; store haplotype data under `annotation/format/HAP` (not `EDS`); infer sample and SNP counts from a 2D haplotype matrix correctly.
+ When `create_gds` omits genotype, store an all-zero dummy genotype matrix in the GDS (with a console message) so that `GbsrGenotypeData` validation succeeds; use dosage or haplotype for downstream analyses.

Changes in version 0.4.20 (2026-03-31)
+ Minor bug fix.

Changes in version 0.4.16 (2026-01-23)
+ Add script to notify users if null_formula is specified without fixed_effect.

Changes in version 0.4.13 (2025-12-22)
+ Add argument to add fixed effects in the regression model.

Changes in version 0.4.12 (2025-10-20)
+ Minor bug fix in getGenoPerMarker()

Changes in version 0.4.11 (2025-08-13)
+ Add argument to change the null model in the regression.

Changes in version 0.4.9 (2025-08-02)
+ Bug fix in recalcAssoc(), which was accidentally introduced in the bug fix at 2025-05-16.

Changes in version 0.4.7 (2025-05-16)
+ Minor bug fix in recalcAssoc().

Changes in version 0.4.6 (2025-04-15)
+ Minor bug fix in listCandidate().

Release version 0.4.4 (2025-04-04)

