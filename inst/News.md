Changes in version 0.5.0 (2026-05-16)
+ Association results (scan, peakcall, recalc, candidate) are stored in a companion Parquet dataset (`{gds}.lazygas/`) by default instead of the GDS `lazygas/` subtree.
+ `buildLazyGas()` gains `lazygas_store` (`"parquet"`, `"sqlite"`, `"gds"`, `"auto"`).
+ `importLazyGasResults()` migrates legacy GDS-stored results into the companion store.

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

