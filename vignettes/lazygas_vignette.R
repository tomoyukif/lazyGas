## ----message=FALSE, results="hide", warning=FALSE-----------------------------
library(lazyGas)


## ----message=FALSE, results="hide", warning=FALSE-----------------------------
gds_fn <- system.file("extdata", "sample.gds", package = "lazyGas")
lg <- buildLazyGas(gds_fn = gds_fn, load_filter = TRUE, overwrite = FALSE)


## ----eval=FALSE---------------------------------------------------------------
# lg <- importLazyGasResults(lg)


## ----message=FALSE, results="hide", warning=FALSE-----------------------------
# Set a temporary file path for a GDS file.
temp_dir <- tempdir()
sample_gds <- tempfile("sample", temp_dir, ".gds")

# Now we assume genotype and marker data were provided as a matrix or 
# data.frame that have been loaded from CSV files, although here we use data 
# obtained from the GDS file for ease of handling as example usage.
genotype <- getGenotype(object = lg, node = "raw")
dosage <- getGenotype(object = lg, node = "dosage")
haplotype <- getHaplotype(object = lg)
snp.chromosome <- getChromosome(object = lg)
snp.position <- getPosition(object = lg)
snp.allele <- getAllele(object = lg)
sample.id <- getSamID(object = lg)
snp.id <- getMarID(object = lg)

# The following values can be NULL: snp.id, snp.rs.id, and snp.allele
# Either of the following values should be specified at least: genotype, 
# haplotype, or dosage.
create_gds <-  list(genotype = genotype,
                    sample.id = sample.id,
                    snp.id = snp.id,
                    snp.rs.id = NULL,
                    snp.chromosome = paste0("chr", sprintf("%02d", as.numeric(snp.chromosome))),
                    snp.position = snp.position,
                    snp.allele = snp.allele,
                    haplotype = haplotype,
                    dosage = dosage)
lg <- buildLazyGas(gds_fn = sample_gds, 
                   create_gds = create_gds)


## -----------------------------------------------------------------------------
pheno_fn <- system.file("extdata", "pheno.csv", package = "lazyGas")
pheno <- read.csv(file = pheno_fn)


## -----------------------------------------------------------------------------
lg <- assignPheno(object = lg,
                  pheno = pheno,
                  rename = "Fruit weight")


## -----------------------------------------------------------------------------
# getPheno() returns a list of phenotype data. 
pheno <- getPheno(object = lg)
for(i in seq_along(pheno$pheno_names)){
  p <- suppressMessages(plotPheno(object = lg, 
                                  pheno = i,
                                  xlab = pheno$pheno_names[i]))
  print(p)
}


## -----------------------------------------------------------------------------
g <- getGenoPerMarker(object = lg, geno_format = "genotype")
print(g)


## -----------------------------------------------------------------------------
g <- getGenoPerMarker(object = lg, geno_format = "haplotype")
print(g)


## -----------------------------------------------------------------------------
g <- getGenoPerMarker(object = lg, geno_format = "dosage")
print(g)


## -----------------------------------------------------------------------------
# Conversion function for dosage data assuming additive and dominant effects
conv_fun <- function(g){
  add <- g
  dom <- as.numeric(g == 1)
  out <- data.frame(add = add, dom = dom)
  return(out)
}

# For the model matrix above, the formula can be the following.
formula <- "add + dom"


## -----------------------------------------------------------------------------
conv_fun <- function(g){
  add <- colSums(g)
  dom <- as.numeric(g == 1)
  out <- data.frame(add = add, dom = dom)
  return(out)
}

# For the model matrix above, the formula can be the following.
formula <- "add + dom"


## -----------------------------------------------------------------------------
# Conversion function for haplotype data assuming additive and interaction effects
conv_fun <- function(g){
  hap1 <- colSums(g == 1)
  hap2 <- colSums(g == 2)
  hap1_hap2 <- as.numeric(hap1 == 1 & hap2 == 1)
  out <- data.frame(hap2 = hap2, hap1_hap2 = hap1_hap2)
  return(out)
}
formula <- "hap2 + hap1_hap2"



## -----------------------------------------------------------------------------
# The n_levels argument specifies the number of dosage levels.
# For diploid samples, n_levels should be set to 3: 0, 1, and 2 represent 
# null-, single-, and duplex dosages, respectively.
conv_fun <- makeConvFun(geno_format = "dosage", n_levels = 3)
formula <- "add + dom"


## -----------------------------------------------------------------------------
conv_fun


## ----eval=FALSE---------------------------------------------------------------
# # Do not run.
# conv_fun <- function(g){
#   hap1 <- colSums(g == 1)
#   hap2 <- colSums(g == 2)
#   hap3 <- colSums(g == 3)
#   hap4 <- colSums(g == 4)
#   hap1_hap2 <- as.numeric(hap1 == 1 & hap2 == 1)
#   hap1_hap3 <- as.numeric(hap1 == 1 & hap3 == 1)
#   hap1_hap4 <- as.numeric(hap1 == 1 & hap4 == 1)
#   hap2_hap3 <- as.numeric(hap2 == 1 & hap3 == 1)
#   hap2_hap4 <- as.numeric(hap2 == 1 & hap4 == 1)
#   hap3_hap4 <- as.numeric(hap3 == 1 & hap4 == 1)
#   out <- data.frame(hap2 = hap2, hap3 = hap3, hap4 = hap4,
#                     hap1_hap2 = hap1_hap2, hap1_hap3 = hap1_hap3,
#                     hap1_hap4 = hap1_hap4, hap2_hap3 = hap2_hap3,
#                     hap2_hap4 = hap2_hap4, hap3_hap4 = hap3_hap4)
#   return(out)
# }
# formula <- paste0("hap2 + hap3 + hap4 + hap1_hap2 + hap1_hap3",
#                   " + hap1_hap4 + hap2_hap3 + hap2_hap4 + hap3_hap4")


## -----------------------------------------------------------------------------
scanAssoc(object = lg,
          formula = formula,
          conv_fun = conv_fun,
          geno_format = "dosage")


## ----eval = FALSE-------------------------------------------------------------
# scanAssoc(object = lg,
#           geno_format = "genotype",
#           method = "mlm")


## -----------------------------------------------------------------------------
# As an example, extract p values from the LazyGas object and assign them to the object again.
scan <- lazyData(object = lg, dataset = "scan", pheno = "Fruit weight")
lg <- assignPvalues(object = lg, pheno_name = "Fruit weight", 
                    p_values = scan$P.model, 
                    geno_format = "dosage", conv_fun = conv_fun, 
                    formula = formula)


## -----------------------------------------------------------------------------
pheno <- getPheno(object = lg)
for(i in seq_along(pheno$pheno_names)){
  p <- plotManhattan(object = lg, pheno = i)
  print(p)
}

# lazyData() function can be used to extract the results.
scan <- lazyData(object = lg, dataset = "scan", pheno = pheno$pheno_names[1])
scan[which.max(scan$negLog10P), ]


## ----message=FALSE------------------------------------------------------------
callPeakBlock(object = lg, signif = 0.05, threshold = 0.8)


## ----warning=FALSE------------------------------------------------------------
for(i in seq_along(pheno$pheno_names)){
  p <- plotPeaks(object = lg, pheno = i)
  print(p)
}


## -----------------------------------------------------------------------------
recalcAssoc(object = lg, 
            n_threads = 1, 
            refine_position = FALSE, 
            grouping_threshold = 0.05)


## ----warning=FALSE------------------------------------------------------------
for(i in seq_along(pheno$pheno_names)){
  p <- plotPeaks(object = lg, pheno = i, recalc = TRUE)
  print(p)
}


## -----------------------------------------------------------------------------
for(i in seq_along(pheno$pheno_names)){
  out <- haploPlot(object = lg, pheno = i, recalc = FALSE)
  # As the output is a list of ggplot objects, only the first plot is shown here. 
  print(out[[1]]) 
}


## -----------------------------------------------------------------------------
for(i in seq_along(pheno$pheno_names)){
  out <- haploPlot(object = lg, pheno = i, recalc = TRUE)
  # As the output is a list of ggplot objects, only the first plot is shown here. 
  print(out[[1]]) 
}


## ----eval=FALSE---------------------------------------------------------------
# # Do not run.
# gff_fn <- "path/to/gff/gene_annotation.gff"
# gff <- rtracklayer::import.gff(gff_fn)
# 
# snpeff_fn <- "path/to/SnpEff/output/variants.snpeff.vcf"
# gds_fn <- sub("\\.vcf", ".gds", snpeff_fn)
# snpeff2gds(vcf_fn = snpeff_fn, out_fn = gds_fn)
# snpeff_gds <- open_snpeff(gds_fn = gds_fn)
# 
# ann <- read.csv("path/to/gene_annotation.csv")
# names(ann)[1] <- "Gene_ID"
# listCandidate(object = lg, gff = gff, snpeff = snpeff_gds, ann = ann, recalc = TRUE)


## ----message=FALSE, warning=FALSE---------------------------------------------
gff_fn <- system.file("extdata", "demo_annotation.gff", package = "lazyGas")
vcf_fn <- system.file("extdata", "demo_snpeff.vcf", package = "lazyGas")
ann_fn <- system.file("extdata", "demo_ann.csv", package = "lazyGas")

have_demo_ann <- nzchar(gff_fn) && file.exists(gff_fn) &&
  nzchar(vcf_fn) && file.exists(vcf_fn)

if (have_demo_ann) {
  gff <- rtracklayer::import.gff(gff_fn)
  ann <- read.csv(ann_fn, stringsAsFactors = FALSE)
  snpeff_gds_fn <- tempfile(pattern = "lazygas_snpeff_", fileext = ".gds")
  snpeff2gds(vcf_fn = vcf_fn, out_fn = snpeff_gds_fn, verbose = FALSE)
  snpeff_gds <- open_snpeff(snpeff_gds_fn)
  listCandidate(
    object = lg,
    gff = gff,
    snpeff = snpeff_gds,
    ann = ann,
    recalc = TRUE
  )
  candidate <- lazyData(object = lg, dataset = "candidate", pheno = "Fruit weight")
  knitr::kable(candidate[, c("peak_ID", "Gene_ID", "Gene_chr", "negLog10P", "Name")])
} else {
  message("Demo annotation files are not installed; skipping listCandidate() example.")
}


## -----------------------------------------------------------------------------
cand_demo <- data.frame(
  peak_ID = 1L,
  Gene_ID = c("g1", "g2", "g3"),
  Gene_chr = "1",
  Gene_start = 1:3,
  dist2peak = 0,
  negLog10P = 3,
  Description = c(
    "fruit weight development",
    "root hair elongation",
    "cell wall biosynthesis"
  ),
  stringsAsFactors = FALSE
)
searchCandidateGenes(
  candidate = cand_demo,
  query = "fruit weight",
  mode = "keyword",
  keyword_match = "all"
)


## ----eval=FALSE---------------------------------------------------------------
# searchCandidateGenes(object = lg, pheno = "Fruit weight",
#                        query = "fruit ripening carbohydrate",
#                        mode = "both", min_score = 0.15, top_n = 50)


## -----------------------------------------------------------------------------
pheno <- getPheno(object = lg)
sapply(pheno, head)


## -----------------------------------------------------------------------------
scan <- lazyData(object = lg, dataset = "scan", pheno = pheno$pheno_names[1])
head(scan)


## -----------------------------------------------------------------------------
peakcall <- lazyData(object = lg, dataset = "peakcall", pheno = pheno$pheno_names[1])
head(peakcall)


## -----------------------------------------------------------------------------
recalc <- lazyData(object = lg, dataset = "recalc", pheno = pheno$pheno_names[1])
head(recalc)


## -----------------------------------------------------------------------------
groups <- lazyData(object = lg, dataset = "groups", pheno = pheno$pheno_names[1])
head(groups)










## ----eval=FALSE---------------------------------------------------------------
# makeInteractiveSummary(
#   object = lg,
#   pheno = "Fruit weight",
#   out_fn = "sample_summary.html",
#   what = c("scan", "scan_png", "peakcall", "recalc", "groups", "candidate")
# )


## ----eval=FALSE---------------------------------------------------------------
# makeInteractiveDashboard(
#   object = lg,
#   pheno = "Fruit weight",
#   gff = gff,
#   out_fn = "lazygas_dashboard.html",
#   snpeff = snpeff_gds,
#   ann = ann,
#   recalc = TRUE,
#   what = c("scan_png", "peakcall", "recalc", "candidate")
# )


## ----eval=FALSE---------------------------------------------------------------
# runDashboardDemo(out_dir = "demo_output")
# # Opens: demo_output/lazygas_dashboard.html


## Rscript inst/demo/run_dashboard_demo.R demo_output

## -----------------------------------------------------------------------------
sessionInfo()

