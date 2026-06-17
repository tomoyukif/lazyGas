Changes in version 0.6.0 (2026-06-17)
+ Refactor package structure by splitting the monolithic `R/01_functions.R` into focused modules (`scan`, `peakcall`, `recalc`, `plot`, `candidate`, `store`, `interactive`, and related helpers).
+ Expand user-facing docs and examples, including README/vignette updates and additional man pages for the dashboard, variant viewer, and import flow.
+ Add CI and packaging scaffolding (`.github/workflows/R-CMD-check.yaml`, `.Rbuildignore`, `.gitignore`, `LICENSE`) and refresh NAMESPACE/collation metadata.
+ Add and extend tests for full pipeline/integration paths, variant viewer, and candidate-gene search behavior.
+ Bundle demo resources and generated outputs for dashboard/annotation workflows.

+ Add `runLazyGas()` one-command pipeline orchestration with resume support and pipeline history in the companion store.
+ Add GWAS QC helpers: `calcGenomicInflation()`, `plotQQ()`, and `summarizeGWASQC()` (also available in interactive reports via `what = "qq"` / `"qc"`).
+ Add cross-trait peak analysis: `clusterCrossTraitPeaks()`, `summarizeCrossTraitPeaks()`, and `plotMultiTraitOverview()` (distance + peak-marker genotype correlation clustering).
+ Add fine-mapping helpers: `conditionalAssoc()`, `calcCredibleSet()`, and `plotCredibleSet()` (Wakefield ABF credible sets).
+ Add ortholog annotation helpers: `annotateOrthologs()`, `summarizeOrthologMatches()`, and `plotOrthologSummary()`.
+ Add `runLazyGasShiny()` and templates under `inst/shiny/` and `inst/quarto/`.
+ Extend `lazyData()` with `qc`, `multitrait`, `conditional`, `credible_set`, and `pipeline` datasets.
+ Exclude large local development datasets under `R/dev/` from version control (`.gds`, large `.csv`/`.gff`/`.pdf` files).

Changes in version 0.5.0 (2026-06-17)
+ Association results (scan, peakcall, recalc, candidate) are stored in a companion Parquet dataset (`{gds}.lazygas/`) by default instead of the GDS `lazygas/` subtree.
+ `buildLazyGas()` gains `lazygas_store` (`"parquet"`, `"sqlite"`, `"gds"`, `"auto"`).
+ `importLazyGasResults()` migrates legacy GDS-stored results into the companion store.
+ `searchCandidateGenes()` filters and ranks candidate genes by keyword and/or semantic similarity (optional **text2vec**).
+ `getVariantViewerData()` and `plotVariantViewer()` plot gene structure and SnpEff-annotated variants along a peak.
+ `makeInteractiveDashboard()` extends `makeInteractiveSummary()` with clickable candidate `Gene_ID` links to embedded variant viewers.
+ `runDashboardDemo()` runs the full sample-data workflow and writes a demo HTML dashboard; bundled files in `inst/extdata/` (`demo_annotation.gff`, `demo_snpeff.vcf`, `demo_ann.csv`).
+ Static Manhattan plots in HTML reports (`scan_png`) require the suggested package **base64enc**.
+ Fix `assignPvalues()` phenotype indexing when multiple traits are present.
+ Fix haplotype dimension when creating a GDS via `create_gds`.
+ Performance: cache per-chromosome genotypes in peak calling and recalc; cache marker coordinates and SnpEff GDS index; batch companion-store scan metadata updates.

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

