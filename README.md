# lazyGas

`lazyGas` is an R package for genetic association studies (GWAS/QTL mapping). It
loads genotypes from GDS files (via [GBScleanR](https://github.com/tomoyukif/GBScleanR)),
stores association outputs in a companion Parquet dataset, and provides peak
calling, candidate-gene listing, variant visualization, interactive HTML reports,
cross-trait analysis (including optional credible-set overlap across traits),
fine-mapping helpers, a thin **pipeline runner** Shiny app
(`runLazyGasRunner()`), and a **phenotype explorer** Shiny app (`runLazyGasExplorer()`;
GWAS / locus plots plus optional local LLM support).

## Environment setup

Follow these steps before running analyses or demos.

### 1. System requirements

| Component | Required | Notes |
|-----------|----------|-------|
| R | Yes | R ≥ 4.1 recommended (`tools::R_user_dir()` for cache paths) |
| GNU make, C++11 | Yes | Listed in `SystemRequirements` for dependency builds |
| [GBScleanR](https://github.com/tomoyukif/GBScleanR) | Yes | Hard dependency; install from GitHub before `lazyGas` |
| Network | First install | CRAN/Bioconductor/GitHub downloads; Ollama image pull if using Apptainer LLM |

### 2. Install core R packages

```R
install.packages(c(
  "ggplot2", "arrow", "jsonlite", "dplyr", "cowplot",
  "parameters", "gaston", "vcfR", "plotly", "htmltools",
  "reactable", "devtools"
))
```

Install **GBScleanR** from its GitHub repository (see link above), then install
**lazyGas**:

```R
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("tomoyukif/lazyGas", build_vignettes = TRUE)
```

For a local checkout, load the package without installing:

```R
devtools::load_all("/path/to/lazyGas")
```

### 3. Optional R packages (by feature)

| Package | Used for |
|---------|----------|
| `shiny` | `runLazyGasRunner()` pipeline UI; `runLazyGasExplorer()` phenotype explorer |
| `base64enc` | Static Manhattan PNGs in HTML reports / explorer (`what = "scan_png"` style) |
| `processx` | Background Ollama server from R (agentic / explorer extras) |
| `httr2`, `curl` | Local LLM HTTP calls (optional extras) |
| `rentrez` | Literature evidence in phenotype explorer (optional) |
| `BiocStyle`, `knitr`, `rmarkdown` | Building vignettes |
| `testthat` | Running package tests |

`plotly`, `htmltools`, and `reactable` are **Imports** (interactive reports and
explorer plots). Install as needed for Suggested features, for example:

```R
install.packages(c("shiny", "base64enc", "processx", "httr2", "curl"))
```

### 4. Optional local LLM (phenotype explorer)

Ranking and evidence collection work without an LLM. For natural-language query
parsing and explanations, use either:

**A. Apptainer / Singularity (recommended on HPC; no host `ollama` binary)**

- Install [Apptainer](https://apptainer.org/) or Singularity and ensure
  `apptainer` or `singularity` is on `PATH` (or set `LAZYGAS_APPTAINER_BIN`).
- In R:

```R
library(lazyGas)
configureLazyGasOllama(model = "llama3.2:3b")  # pull image + install CLI wrapper
startLocalLLM(model = "llama3.2:3b")           # serve + pull model
```

Models and cache default to `lazyGasOllamaHome()` (override with
`LAZYGAS_OLLAMA_HOME`). Shell alternative:

```bash
bash inst/apptainer/build_ollama_sif.sh
export LAZYGAS_OLLAMA_SIF="$HOME/.cache/lazygas/ollama/ollama.sif"
```

**B. Host Ollama**

- Install [Ollama](https://ollama.com/) and ensure `ollama` is on `PATH`.
- In R: `startLocalLLM(model = "llama3.2:3b")`.

In the Shiny explorer: **Setup Ollama (Apptainer)** → **Start Ollama** → **Pull model**.

### 5. Verify the installation

```R
library(lazyGas)
out_dir <- tempfile("lazygas_check_")
runDashboardDemo(out_dir = out_dir)
# Expect: {out_dir}/lazygas_dashboard.html
```

Run tests from a source checkout:

```bash
Rscript -e 'devtools::test()'
```

## Vignette

```R
browseVignettes(package = "lazyGas")
```

The vignette walks through the full workflow and repeats the environment setup
with more detail (companion storage, optional LLM, demos).

## Demos

Bundled sample data lives under `inst/extdata/` (`sample.gds`, phenotype CSVs,
demo annotation, SnpEff VCF, expression tables, and ortholog table). After
installing or `devtools::load_all()`, run end-to-end demos without your own files.

**Dashboard demo** (variant viewer + interactive summary):

```R
runDashboardDemo(out_dir = "demo_output")
# -> demo_output/lazygas_dashboard.html
```

**MVP feature demo** (v0.6+ APIs: `runLazyGas()`, GWAS QC, multi-trait peaks,
optional `compareCrossTraitCredibleSets()`, fine-mapping):

```R
runMvpDemo(out_dir = "demo_output/mvp")
# -> demo_output/mvp/lazygas_mvp_report.html
```

**Phenotype explorer demo** (v0.7+ evidence ranking + optional local LLM;
v0.9.0+ GWAS / locus dashboard plots in the explorer; v0.9.1+ separate pipeline
runner):

```R
runExplorerDemo(out_dir = "demo_output/explorer")
# -> demo_output/explorer/phenotype_explorer_report.md

runLazyGasRunner()       # Shiny UI — run scan → candidate (paths on disk)
runLazyGasExplorer()     # Shiny UI — explore / rank (after runner or demo)
# Tabs: GWAS overview, Local LLM, Ranked genes, Evidence, Locus / variants, Chat
```

**Two-app workflow:** use **`runLazyGasRunner()`** once to build companion-store
results, then **`runLazyGasExplorer()`** to explore and rank candidates.

After **Load project**, use **GWAS overview** to set the active phenotype and
inspect phenotype/Manhattan/peak views. Configure optional local LLM settings
on **Local LLM**. Use **Ranked genes** for phenotype-query inputs and
**Rank candidates**. Then select a gene and open **Locus / variants**
(set a GFF path in the sidebar for the variant viewer; demo default is
`inst/extdata/demo_annotation.gff`).

During long-running **Rank candidates** / **Generate explanation** / **Ask**
actions, the explorer shows progress and temporarily disables action buttons to
avoid duplicate submissions.

In the explorer, enter the **filesystem path** to the GDS (not an uploaded copy).
The companion folder `sample.lazygas` must sit in the same directory.

```R
paths <- resolveLazyGasPaths("demo_output/explorer/sample.gds")
lg <- buildLazyGas(paths$gds, companion_path = paths$companion)
lg <- restorePhenoFromStore(lg)
```

From the shell at the package root:

```bash
Rscript inst/demo/run_dashboard_demo.R demo_output
Rscript inst/demo/run_mvp_demo.R demo_output/mvp
Rscript inst/demo/run_explorer_demo.R demo_output/explorer
```

## Usage

### Loading input data

Genotypes stay in the GDS file; association results go to a companion folder
`path/to/your/input.lazygas/` (Parquet, default since v0.5.0):

```R
library(lazyGas)
library(GBScleanR)

gds_fn <- "path/to/your/input.gds"
lg <- buildLazyGas(gds_fn = gds_fn, load_filter = TRUE, overwrite = FALSE)
# Companion store is always Parquet: path/to/your/input.lazygas/
```

### Phenotype data

Sample id column must be `id`, `ID`, or `sample_id` (`id`/`ID` preferred).

```R
pheno <- read.csv("path/to/phenotype_data.csv")
lg <- assignPheno(object = lg, pheno = pheno, rename = "Trait1")
```

### Association scan

Continuous phenotypes are standardized internally for regression; stored
`Coef.*` columns (e.g. `Coef.add`) are rescaled to **original phenotype units**
(change in trait per genotype unit). Binary traits use the GLM log-odds scale.

```R
conv_fun <- makeConvFun(geno_format = "dosage", n_levels = 3)
scanAssoc(
  object = lg,
  formula = "add + dom",
  conv_fun = conv_fun,
  geno_format = "dosage"
)
```

### Peak calling, recalc, and candidates

```R
callPeakBlock(object = lg, signif = 0.05, threshold = 0.8)
recalcAssoc(object = lg, refine_position = FALSE, n_threads = 1L)

gff <- rtracklayer::import.gff("path/to/genes.gff")
listCandidate(object = lg, gff = gff, recalc = TRUE)
```

Or run the core steps in one call:

```R
lg <- runLazyGas(
  object = lg,
  steps = c("scan", "peakcall", "recalc", "candidate"),
  gff = gff,
  formula = "add + dom",
  conv_fun = conv_fun,
  geno_format = "dosage"
)
```

### Extract results

```R
scan <- lazyData(object = lg, dataset = "scan", pheno = "Trait1")
peakcall <- lazyData(object = lg, dataset = "peakcall", pheno = "Trait1")
recalc <- lazyData(object = lg, dataset = "recalc", pheno = "Trait1")
```

`haploPlot()` titles include peak `-log10P` and `Coef.add` (original units for
continuous traits).

### Multi-trait peaks and optional CS overlap

Not part of the default `runLazyGas()` steps — call when you have multiple
phenotypes and want a cross-trait view:

```R
plotMultiTraitOverview(lg, recalc = TRUE)
clusterCrossTraitPeaks(lg, dist_threshold = 500000, r2_threshold = 0.6)

# Optional: label shared vs independent signals via credible-set overlap
cs_ov <- compareCrossTraitCredibleSets(lg, recalc = TRUE)
cs_ov$pairs
lazyData(lg, dataset = "multitrait", kind = "cs_overlap")
```

### Ranking and reports (this release / master)

This build is the **non-LLM** core. Prefer:

- `phenotypeQuery()` + `rankPhenotypeCandidates()` for ranking
- `explainPhenotypeCandidates()` for a Markdown template summary
  (`save = TRUE` writes under `{stem}.lazygas/reports/`)
- `makeInteractiveDashboard()` / `makeInteractiveSummary()`, or
  `runLazyGas(..., steps = c(..., "dashboard"|"summary"))` for HTML
  (default under `{stem}.lazygas/reports/{pheno}_dashboard.html`)

`listCandidate()` always writes both the wide candidate table and the thin
simple-candidate sidecar. The old `output_list` argument is removed.

LLM HTML reports (`llm_report()`, explore `work_dir`, B5 grounding, etc.) live
on the **agentic** branch / LLM-enabled releases—not in this package build.

### Re-running after input changes

`runLazyGas(resume = TRUE)` skips steps whose companion sections already exist.
If you change GFF, annotation, SnpEff, or similar inputs, **delete the affected
companion sections** (e.g. `candidate/`, `snpeff/`) or call with
`resume = FALSE`, then re-run. Successful `listCandidate()` / ranking writes
`inputs_manifest.json` in the companion store; a mismatch against current inputs
emits a **warning** (it does not auto-delete or re-run).

### Companion layout (target)

| Layer | Location |
|-------|----------|
| Genotypes | `.gds` |
| Analysis results | `{stem}.lazygas/` (Parquet only) |
| HTML / Markdown reports | `{stem}.lazygas/reports/` (defaults) |

Eval demos under `ai_assited_genetics/demo/` keep inputs in `demo_data/`,
companions in `stores/`, and run scripts under `run/`.

## Contributing

Issues and pull requests are welcome on
[GitHub](https://github.com/tomoyukif/lazyGas).

## License

GPL-3 — see [LICENSE](LICENSE).
