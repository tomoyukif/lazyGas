# lazyGas release history

Newest entries first. Version numbers follow semantic ordering; dates reflect
the release or documentation update for that version.

---

## Changes in version 0.10.0 (2026-10-13) — interactive scores + finemap CS coords

### Interactive summary: phenotype rank scores
+ `makeInteractiveSummary(what = "score")` shows the phenotype-rank table.
  With `what = c("candidate", "score")`, score columns are joined onto the
  candidate list (by `Gene_ID` + `peak_ID` when present). `evidence_json` /
  `query_id` stay out of the HTML display.
+ Inactive `score_*` channels are forced to **0** when ranking and when
  reading `lazyData(..., dataset = "phenotype_rank")` (active sources from
  `phenotypeRank` metadata).
+ Finemap PIP mapping: `.e_enrich_cs_coords()` fills missing Chr/Pos on
  credible-set rows from recalc/peakcall blocks when needed.

### Hardening track note
+ without-LLM I/O + scoring hardening on `master` is complete for this pass.
  Layer B Benchmark runner remains deferred (see
  `ai_assited_genetics/manuscript/lazygas/nollm_lazygas_benchmark_plan.md`).

---

## Changes in version 0.10.0 (2026-10-09) — I/O + scoring hardening (master)

### Multi-trait CS overlap (optional)
+ Add **`compareCrossTraitCredibleSets()`**: after
  `clusterCrossTraitPeaks()`, compare Wakefield ABF credible sets across
  traits and label pairs (`shared-compatible`, `possible-shared`,
  `independent-nearby`, `inconclusive`). Opt-in only — not a default
  `runLazyGas()` step. Results store as `multitrait/cs_overlap`
  (`lazyData(..., dataset = "multitrait", kind = "cs_overlap")`).
  See `?compareCrossTraitCredibleSets` and the vignette multi-trait section.

### Companion store and candidates
+ Companion store is **Parquet only** (`lazygas_store = "parquet"`). Modes
  `sqlite`, `gds`, and `auto` are rejected.
+ `listCandidate()` / `runLazyGas()` no longer take `output_list`. Wide
  candidate tables are **always** written (simple candidate sidecar remains).
+ `inputs_manifest.json` records fingerprints for candidate inputs; mismatch
  against current inputs warns (does not auto-delete). After changing GFF/ann,
  delete the relevant companion section or use `resume = FALSE`.
+ Default HTML from `runLazyGas` summary/dashboard steps goes to
  `{companion}/reports/{pheno}_*.html`.
+ `assignPheno()` accepts `sample_id` as well as `id`/`ID`.
+ Removed narrative / interpretation helpers from this non-LLM build:
  `explainPhenotypeCandidates()`, `answerPhenotypeQuestion()`,
  `reinterpretPeakCandidates()`, `geneSetForReport()`, related report
  helpers, and the explorer Explanation tab. Ranking returns score
  tables; users interpret. HTML dashboards remain available.

### Scoring / ranking (without LLM)
+ `phase1DefaultSources()` is **`snpeff` + `gwas_dist` + `gwas_peak`**.
  The old combined `gwas` channel is split:
  - **gwas_dist**: \((D_{\max}-|d|)/D_{\max}\) per attributed marker; gene =
    max (lead distance proximity within the peak).
  - **gwas_peak**: \(\mathrm{negLog10P}/\mathrm{negLog10P}_{lead}\) per
    attributed marker; gene = max.
+ Annotation, expression, and other channels remain available when passed in
  `sources`.
+ Removed `use_keyword` / `use_semantic` arguments and YAML
  `agents.ranking` flags.
+ Channel raw scores stay on `[0, 1]`; **no** per-peak min–max.
  Inactive `score_*` columns are **0**. YAML default weights include
  finemap / literature / ortholog at **0.50**.
+ Annotation (opt-in): tiered hits via `phenotypeQuery(action_phrases=...,
  noun_phrases=..., words=...)` → scores **1 / 0.5 / 0.25 / 0**.
+ Literature / ortholog (opt-in): undirected path reachability on
  `literature_graph` / `ortholog_graph` (`from`,`to`) → **1 / 0**.
  PubMed count and ortholog “presence floor 0.5” removed. Demo:
  `demo_ortholog_graph.csv`, `demo_literature_graph.csv`.
+ Expression channel: missing tissue/stage query terms score **0**
  (previously a flat 0.25 baseline).
+ Removed the post-normalization **floor 0.1** on channel scores.

### Docs
+ README documents default sources/weights, annotation tiers, and network
  graphs. `llm_report()` belongs to the agentic / LLM-enabled line, not
  this build.

---

## Changes in version 0.9.2 (2026-10-09)

### AI HTML report: fix Genes in candidate region = 0
+ **Genes in candidate region** counts **CS∩GFF** from E-path `gff_windows`
  when available (genes whose GFF windows overlap the peak CS / report
  region). It is distinct from Candidate list rows and from Evidence
  `gene_select` size.
+ Missing rank `Gene_end` no longer treats `Gene_start` as a zero-width
  gene (that falsely reported 0 while the Candidate list still listed
  genes, e.g. Amylose / Wx). README and vignette document the distinction.

---

## Changes in version 0.9.2 (2026-10-08)

### AI HTML report: Gene_ID → Evidence links; region gene count
+ Hypothesis cards and Candidate list Gene_ID cells link to
  `#evidence-gene-<Gene_ID>` when that gene has an Evidence block on the same
  peak; otherwise they stay plain text.
+ Candidate list shows **Genes in candidate region** (CS∩GFF via E-path
  `gff_windows`; see 2026-10-09 fix), distinct from Candidate rows and from
  Evidence gene-set size (`gene_select`).
+ `formatPeakReinterpretation(..., evidence_gene_ids=)` enables the Hypothesis
  card links; Evidence gene headings carry matching HTML `id` anchors.

---

## Changes in version 0.9.2 (2026-09-14)

### recalcAssoc: keep correlated peak groups; narrow with CS
+ `recalcAssoc()` groups redundant peak leads but **retains member peak-block
  markers** only when lead `negLog10P` fold (taller/shorter) is strictly below
  `group_retain_fold` (default 2). Same-chromosome passers merge under the
  representative peak ID; cross-chromosome passers stay as **independent peaks**.
  Relationships remain in `recalc/groups` regardless of fold.
+ With `refine_position = TRUE`, same-chr fold-passing members merge into the
  refined group block; cross-chr fold-passing members remain independent.
+ [calcCredibleSet()] uses these per-peak (same-chr) recalc blocks when
  `recalc = TRUE`. Reports list **Related peaks (grouped)** from `groups`.

### Phase 2 F peak reinterpretation (hypothesis cards)
+ Add `reinterpretPeakCandidates()` / `PeakReinterpretation`: per-peak
  proximity, function (plan-alpha strong/weak gate), impact, expression,
  and composite champions with side-by-side conflicts (no primary label).
+ Distance gate defaults **off**; function seat uses query trait keywords +
  synonym phrases only (related phrases are not strong); empty function seat
  is `NA` ("no function hypothesis").
+ `llm_report()` defaults to `reinterpret = TRUE` and
  `gene_select = "union_all"` (composite Candidate top N ∪ card
  champion/support); HTML gains **Hypothesis cards** before the candidate
  list. Each Evidence gene opens with a **Selected by** note
  (composite ranking / hypothesis cards / both).
+ Companion store section `phenotype_reinterpret/` and
  `lazyData(..., "phenotype_reinterpret")`.

### Phase 2 E pipeline (gene × source evidence)
+ Add working-directory config helpers
  `writeLazyGasExploreConfig()` / `readLazyGasExploreConfig()` and
  `classifyAnnColumns()` (E0; LLM column-name classify when `use_llm=TRUE`).
+ Add gene×source E loop (`R/llm-report-e0.R`, `R/llm-report-e.R`,
  `R/llm-report-id-cache.R`): CS+GFF gene selection (GWAS / QTL), E1–E10
  payloads/templates (E10 = integrated interpretation after Brief summary),
  GO OBO / KEGG REST / Pfam·InterPro ID→name caches,
  InterProScan×SnpEff AA overlap (same protein ID; &lt;80% match warning),
  Yanai τ / relative height, and B5 numeric allowlist grounding.
+ E2 / impact champions count **peak-block-marker-linked** SnpEff only
  (`negLog10P` finite); `llm_report` loads store snpeff and filters per peak
  via `.snpeff_for_peak()`.
+ Per gene, multiple ANN at the **same Chr/Pos/Allele** keep only the
  **worst impact** (HIGH>MODERATE>LOW>MODIFIER) and count as **one** variant;
  report HIGH effect lines are no longer multi-transcript duplicates.
+ `listCandidate()` also writes **simple candidate** Parquet and
  **Gene↔transcript/protein** map (`lazyData(..., "simple_candidate"|"gene_protein_map")`).
+ `llm_report(work_dir=...)` switches evidence to the E path; without
  `work_dir` the legacy peak-batch prose path is unchanged.
+ Ship editable `inst/extdata/expression_coarse_vocab.yaml` for E4 labels.
+ `classifyAnnColumns(use_llm=TRUE)`: normalize common label drift
  (e.g. `class-kegg`→`kegg`) and, if still invalid, run **one LLM scrutiny
  pass** that re-checks the prior JSON against the allowed class list.
+ Shared `.llm_scrutiny_pass()`: constrained LLM steps (E source JSON, B5
  numeric grounding, B2 evidence/JA translation, phenotypeQuery + synonym/
  related expand, annotation relevance) re-send prior output + rules for
  one correction pass before falling back.
+ Credible sets now carry `Chr`/`Pos` from the peak block in
  `calcCredibleSet()`; E path also enriches stored CS via peak-call
  `variant_ID` map so gene selection no longer returns empty evidence.
+ `phenotypeQuery(..., llm_timeout=)` (default 600s): query parse / synonym /
  related expand no longer hard-cap at 90s (Gemma 31B tunnel timeouts).

---

## Changes in version 0.9.2 (2026-07-24)

### Annotation keyword matching
+ Remove **`searchCandidateGenes()`**. Keyword matching helpers remain for
  ranking / `collectGeneEvidence()`.
+ Remove text2vec LSA / semantic cosine scoring from annotation evidence
  (`mode = "semantic"` / `"both"` no longer supported).
+ Rebuild keyword matching: phrase-first terms from `PhenotypeQuery`, English
  stopword filtering, word-boundary / phrase patterns, and
  `matched_keywords` / `unmatched_keywords` in evidence details.
+ Phase 2 A1/A2/A4: `use_llm = TRUE` phrase split + synonym/related expand
  (`synonym_phrases` / `related_phrases` with `from_user`); search terms are
  the union of user / context / synonym / related phrases. Family-level
  unmatched (user phrase with no self/synonym/related hit) is kept for
  reports; Keywords vs Validity prompts no longer use high/low/hit-rate
  language. `use_llm = FALSE` keeps unigram fallback and forces annotation
  LLM relevance off.
+ Update LLM report prompts for matched synonym/related (`from_user`) and
  family unmatched; never print numeric ranking scores or claim a "semantic
  match".
### Explorer UI and workflow polish
+ Reorganize **`runLazyGasExplorer()`** sidebar controls by tab so only
  context-relevant inputs are shown:
  + **GWAS overview**: project loading, phenotype selection, Manhattan mode.
  + **Local LLM**: LLM configuration and Ollama actions.
  + **Ranked genes**: phenotype-query inputs and ranking controls.
  + **Locus / variants** and **Chat / explanation**: tab-specific controls.
+ Add a dedicated **Local LLM** tab in the main panel (before **Ranked genes**)
  that summarizes runtime/model status and diagnostics.
+ Keep phenotype selection anchored to **GWAS overview** as the primary place
  to set active analysis context.

### Progress and interaction behavior
+ Add explicit progress feedback for **Rank candidates**,
  **Generate explanation**, and **Ask** actions.
+ Lock action buttons while long-running rank/chat/explanation tasks are active
  to prevent duplicate submissions and conflicting operations.
+ Remove cancel-button UI after reliability review; retain progress display and
  busy-state locking for stable synchronous execution.

### Packaging and docs
+ Bump package version to **0.9.2**.
+ Refresh **README.md** and vignette wording for updated explorer tabs and
  sidebar behavior.

---

## Changes in version 0.9.1 (2026-07-09)

### Pipeline runner (thin Shiny app)
+ Add **`runLazyGasRunner()`** (`inst/shiny/runner/`) — a separate Shiny app to
  run `buildLazyGas()` → `assignPheno()` → `runLazyGas()` (scan through
  candidate) from filesystem paths, then hand off to **`runLazyGasExplorer()`**
  for exploration and ranking.
+ Explorer sidebar links to **`runLazyGasRunner()`** when companion-store
  results are missing.
+ Bump package version to **0.9.1**; refresh **README.md**, vignette, and
  **`runLazyGasRunner()`** man page.

### Note for upgrades
+ Install Suggested **`shiny`** and **`rtracklayer`** (already required for GFF
  import elsewhere) to use the runner UI:
  `install.packages("shiny")`.

---

## Changes in version 0.9.0 (2026-07-09)

### Phenotype explorer (dashboard plots in Shiny)
+ Embed **GWAS overview** and **Locus / variants** tabs in
  **`runLazyGasExplorer()`**, reusing `plotPheno()`, `plotManhattan()`,
  `plotPeaks()`, `haploPlot()`, and `getVariantViewerData()` — the same
  building blocks as **`makeInteractiveDashboard()`**, without iframe / HTML
  export.
+ **GWAS overview** (available after Load project): phenotype distribution,
  Manhattan (static PNG by default; optional interactive plotly), peakcall,
  recalc peaks, and peak-group table.
+ **Locus / variants** (after ranking + gene selection): haplotype plot for the
  selected peak and on-demand variant viewer when a GFF path is provided
  (demo default: `inst/extdata/demo_annotation.gff`). Stored SnpEff tables are
  preferred; an optional SnpEff GDS path is used only when stored annotations
  are missing.
+ Clear ranking / chat state when the phenotype selection changes so locus
  views stay tied to the active trait.
+ Soft-fail missing companion datasets (scan, peaks, groups, etc.) with UI
  messages instead of hard errors.

### Packaging and documentation
+ Bump package version to **0.9.0**.
+ Declare **`plotly`**, **`htmltools`**, and **`reactable`** in **Imports**
  (already used by interactive HTML reports; now also required for explorer
  plots). **`shiny`** remains in Suggests.
+ Refresh **README.md**, vignette, and `runLazyGasExplorer()` man page for the
  new explorer tabs and dependency layout.
+ Add smoke test that the explorer Shiny app sources the new tab / plot helpers.

### Note for upgrades
+ Install or update Imports: `install.packages(c("plotly", "htmltools", "reactable"))`
  if upgrading an older install that lacked them as hard dependencies.
+ For explorer variant views, set the **GFF path** in the Shiny sidebar (and
  optionally SnpEff GDS). Ranking / Evidence / Chat behave as in v0.7–0.8.

---

## Changes in version 0.8.0 (2026-07-04)

### API cleanup
+ Remove **`annotateOrthologs()`**, **`summarizeOrthologMatches()`**, and
  **`plotOrthologSummary()`**. Ortholog tables remain available as an optional
  evidence channel in **`collectGeneEvidence()`** /
  **`rankPhenotypeCandidates()`** (and explorer demos); dedicated
  user-driven ortholog screening helpers are no longer exported.
+ Remove the legacy Shiny launcher **`runLazyGasShiny()`** and
  `inst/shiny/app.R`. Use **`runLazyGasExplorer()`** (`inst/shiny/explorer/`)
  instead.

### Bug fixes
+ Extend **`assignPvalues()`** with optional **`coef`** and **`any_data`**
  arguments for registering externally computed GWAS results (not merging
  into an existing `scanAssoc()` table). `coef` stores effect columns such as
  `Coef.add`; `any_data` appends extra per-marker columns. Clarify that each
  call writes a new scan table from supplied inputs.
+ Unify **`getGenoPerMarker()`** signature: generic and method both use
  **`marker_index`** (1-based index among valid markers). Legacy
  `marker_id=` is still accepted with a deprecation warning. Positional third
  arguments now correctly select a marker (previously ignored and always
  returned marker 1). Document return shapes: dosage vector; genotype/haplotype
  as `ploidy × nsam` matrices (alleles × samples).

### GWAS effect sizes and visualization
+ Rescale stored scan **`Coef.*` columns to original phenotype units** for
  continuous traits. Phenotypes are still standardized internally for regression;
  coefficients are multiplied by `sd(phenotype)` before writing to the companion
  store. Binary traits are unchanged (GLM log-odds scale). P-values and FDR are
  unaffected.
+ Show additive scan effect (`Coef.add`) and allele direction in **`haploPlot()`**
  titles when scan results are available.
+ Add regression test `tests/testthat/test-scan-coef.R` comparing stored
  coefficients with a manual GLM on non-standardized phenotypes.

### Documentation
+ Bump package version to 0.8.0; refresh `DESCRIPTION` summary text.
+ Rewrite **README.md** for current APIs (`conv_fun`, `runLazyGas()`, explorer
  demos) and add a dedicated **environment setup** section (core and optional
  dependencies, Apptainer/host Ollama, verification steps).
+ Update **vignette** (`lazygas_vignette.Rmd`): environment setup, interpreting
  `Coef.*`, haplotype plot titles, bundled `extdata` table.
+ Remove **`importLazyGasResults()`** and automatic migration from the legacy
  GDS `lazygas/` subtree into the companion store.

### Note for upgrades
+ Re-run **`scanAssoc()`** if you need rescaled `Coef.*` values for results saved
  with lazyGas versions before 0.8.0.
+ Call sites using `annotateOrthologs()` / `summarizeOrthologMatches()` /
  `plotOrthologSummary()` should pass the same table to
  `collectGeneEvidence(..., sources = "ortholog")` or
  `rankPhenotypeCandidates(..., ortholog_table = ...)` instead.
+ Replace `runLazyGasShiny()` with `runLazyGasExplorer()`.

---

## Changes in version 0.7.0 (2026-06-19)

### Phenotype explorer (Phase 4)
+ Add phenotype-guided candidate exploration:
  `phenotypeQuery()`, `collectGeneEvidence()`, `rankPhenotypeCandidates()`,
  and `explainPhenotypeCandidates()`.
+ Add **`runLazyGasExplorer()`** Shiny app (`inst/shiny/explorer/`) and
  **`runExplorerDemo()`** with bundled expression, ortholog, and multi-trait demo
  data (`demo_expression.csv`, `demo_expression_meta.csv`).
+ Extend **`lazyData()`** with `phenotype_rank` and `phenotype_query` datasets;
  companion-store sections `phenotype_query/`, `evidence/cache/`, and
  `phenotype_rank/`.
+ Add **`resolveLazyGasPaths()`** and **`restorePhenoFromStore()`** for reloading
  GDS + companion store and phenotypes in explorer workflows.

### Local LLM (optional)
+ Add official **Apptainer Ollama** setup: `configureLazyGasOllama()`,
  `buildLazyGasOllamaSif()`, `lazyGasOllamaHome()`, `startLocalLLM()`,
  `llmChat()`, and related helpers; assets under `inst/apptainer/`; Shiny
  **Setup Ollama (Apptainer)** button.
+ Host-installed [Ollama](https://ollama.com/) remains supported when `ollama` is
  on `PATH`.

---

## Changes in version 0.6.0 (2026-06-17)

### Package structure and quality
+ Refactor by splitting the monolithic `R/01_functions.R` into focused modules
  (`scan`, `peakcall`, `recalc`, `plot`, `candidate`, `store`, `interactive`,
  and related helpers).
+ Expand user-facing docs and examples; add man pages for dashboard, variant
  viewer, and import flow.
+ Add CI and packaging scaffolding (`.github/workflows/R-CMD-check.yaml`,
  `.Rbuildignore`, `.gitignore`, `LICENSE`); refresh NAMESPACE and Collate.
+ Add and extend tests for full pipeline, variant viewer, and candidate-gene
  search.
+ Bundle demo resources under `inst/extdata/` and demo scripts under `inst/demo/`.
+ Exclude large local development datasets under `R/dev/` from version control.

### Pipeline orchestration
+ Add **`runLazyGas()`** one-command pipeline with **`resume`** support and
  pipeline history in the companion store.

### GWAS QC
+ Add **`calcGenomicInflation()`**, **`plotQQ()`**, and **`summarizeGWASQC()`**
  (also in interactive reports via `what = "qq"` / `"qc"`).

### Multi-trait analysis
+ Add **`clusterCrossTraitPeaks()`**, **`summarizeCrossTraitPeaks()`**, and
  **`plotMultiTraitOverview()`** (genomic distance + peak-marker genotype
  correlation clustering).

### Fine-mapping
+ Add **`conditionalAssoc()`**, **`calcCredibleSet()`**, **`plotCredibleSet()`**,
  **`summarizeConditionalSignals()`**, **`runFineMapping()`**, and
  **`summarizeFineMapping()`** (Wakefield ABF credible sets; interpretation-focused
  summaries documented in the vignette).

### Ortholog evidence (explorer)
+ Bundle `demo_orthologs.csv` for optional ortholog evidence in phenotype
  ranking (`collectGeneEvidence()` / `rankPhenotypeCandidates()`). Standalone
  ortholog-annotation helpers were added earlier in development and later
  removed in favour of the evidence path only (see 0.8.0).

### Demos and UI
+ Add **`runMvpDemo()`** and `inst/demo/run_mvp_demo.R` for multi-trait QC and
  fine-mapping workflows.
+ Add **`runLazyGasShiny()`**; templates under `inst/shiny/` and `inst/quarto/`.

### Data access
+ Extend **`lazyData()`** with `qc`, `multitrait`, `conditional`, `credible_set`,
  `fine_mapping`, and `pipeline` datasets.

---

## Changes in version 0.5.0 (2026-06-17)

### Companion storage (Parquet by default)
+ Store association results (`scan`, `peakcall`, `recalc`, `candidate`, and
  related metadata) in a **companion Parquet dataset** (`{gds_stem}.lazygas/`)
  by default instead of a `lazygas/` subtree inside the GDS file.
+ **`buildLazyGas()`** gains **`lazygas_store`**: `"parquet"` (default),
  `"sqlite"`, `"gds"` (legacy), or `"auto"`.
+ Open GDS connections via **`loadGDS()`** (GBScleanR) in **`buildLazyGas()`**.

### Candidate search and visualization
+ Add **`searchCandidateGenes()`** — keyword and/or semantic similarity filtering
  (optional suggested package **text2vec**).
+ Add **`getVariantViewerData()`** and **`plotVariantViewer()`** for gene structure
  and SnpEff-annotated variants along a peak.
+ Add **`makeInteractiveDashboard()`** — extends **`makeInteractiveSummary()`**
  with clickable candidate **`Gene_ID`** links to embedded variant viewers.
+ Add **`runDashboardDemo()`** and bundled demo annotation files in `inst/extdata/`
  (`demo_annotation.gff`, `demo_snpeff.vcf`, `demo_ann.csv`).

### Fixes and performance
+ Fix **`assignPvalues()`** phenotype indexing when multiple traits are present.
+ Fix haplotype dimension when creating a GDS via **`create_gds`**.
+ Cache per-chromosome genotypes in peak calling and recalc; cache marker
  coordinates and SnpEff GDS index; batch companion-store scan metadata updates.
+ Static Manhattan PNGs in HTML reports (`what = "scan_png"`) require suggested
  package **base64enc**.

---

## Changes in version 0.5.1 (2026-05-18)

+ Skip recreation of existing `lazygas/scan`, `lazygas/peakcall`, `lazygas/recalc`,
  and `lazygas/candidate` GDS nodes on reruns; update peakcall attributes when
  the peakcall node already exists.
+ Document the **`limit_peakcall`** argument in **`callPeakBlock()`**.
+ Fix **`buildLazyGas()`** when creating a GDS via **`create_gds`**: create the
  `annotation/format` folder before writing haplotype or dosage data; store
  haplotype data under `annotation/format/HAP` (not `EDS`); infer sample and SNP
  counts from a 2D haplotype matrix correctly.
+ When **`create_gds`** omits genotype, store an all-zero dummy genotype matrix
  in the GDS (with a console message) so that `GbsrGenotypeData` validation
  succeeds; use dosage or haplotype for downstream analyses.

---

## Changes in version 0.4.20 (2026-03-31)

+ Minor bug fix.

---

## Changes in version 0.4.16 (2026-01-23)

+ Notify users when **`null_formula`** is specified without **`fixed_effect`**.

---

## Changes in version 0.4.13 (2025-12-22)

+ Add **`fixed_effect`** argument to incorporate additional fixed effects in
  **`scanAssoc()`** regression models.

---

## Changes in version 0.4.12 (2025-10-20)

+ Minor bug fix in **`getGenoPerMarker()`**.

---

## Changes in version 0.4.11 (2025-08-13)

+ Add **`null_formula`** argument to **`scanAssoc()`** for custom null models in
  GLM-based association tests.

---

## Changes in version 0.4.9 (2025-08-02)

+ Bug fix in **`recalcAssoc()`**, accidentally introduced in the 2025-05-16 fix.

---

## Changes in version 0.4.7 (2025-05-16)

+ Minor bug fix in **`recalcAssoc()`**.

---

## Changes in version 0.4.6 (2025-04-15)

+ Minor bug fix in **`listCandidate()`**.

---

## Changes in version 0.4.4 (2025-04-04)

First public release.

+ Core GWAS pipeline: **`buildLazyGas()`**, **`assignPheno()`**, **`scanAssoc()`**
  (GLM and mixed-model options), **`callPeakBlock()`**, **`recalcAssoc()`**, and
  **`listCandidate()`** with GFF and optional SnpEff annotation.
+ Genotype formats: dosage, genotype, corrected genotype, and haplotype via
  **`makeConvFun()`** and custom conversion functions.
+ Visualization: Manhattan plots (**`plotManhattan()`**), peak plots
  (**`plotPeaks()`**), haplotype boxplots (**`haploPlot()`**), phenotype plots
  (**`plotPheno()`**).
+ Interactive HTML reports via **`makeInteractiveSummary()`**.
+ SnpEff integration: **`snpeff2gds()`**, **`open_snpeff()`**.
+ External p-values: **`assignPvalues()`** for results from other tools.
+ Peak refinement: optional re-regression around peaks and configurable grouping
  threshold in **`recalcAssoc()`**.
