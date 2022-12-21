# Load functions
source("~/hdd2/softDevel/lazyGeno/R/00_functions.R")

#######################################################
# Build an object holding genotype and phenotype data #
#######################################################
gds_fn <- "~/hdd2/og_gwas/genetics/input/variants/ogCore_variants.gds"
gds <- loadGDS(gds_fn, load_filter = TRUE)

pheno_fn <- "~/hdd2/og_gwas/withNU/gwas/output/2022/phenotypeQC/filtered_phenotype.Rdata"
load(pheno_fn)
sampleID_pheno <- phe$ID
phe <- subset(phe, select = c(TIL_filt, Plant_height_filt))

obj <- buildLazyGWAS(gds = gds,
                     pheno = phe,
                     sampleID_pheno = sampleID_pheno)
obj    # Confirm the object structure

#######################
# Regression analysis #
#######################
out_dir <- "~/hdd2/softDevel/lazyGeno/gwastest"    # Prepare a direcoty to output files
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Execute regression analysis
out_fn <- file.path(out_dir, "test_")
scangwas <- scanGWAS(obj,
                   out_fn = out_fn,
                   nPC = 2,
                   n_core = 1)

########################
# Draw manhattan plots #
########################
signif <- "x$FDR <= 0.05"
out_fn <- file.path(out_dir, "test_")
out_fmt <- "png"
plotManhattan(scangwas, signif = signif, out_fn = out_fn, out_fmt = out_fmt)


##############
# Call peaks #
##############
signif <- "x$FDR <= 0.05"
rsquare <- 0.6
out_fn <- file.path(out_dir, "test_")
peakblock <- callPeakBlock(scangwas, signif = signif,
                           out_fn = out_fn, rsquare = rsquare)

##################
# Execute snpEff #
##################
args <- "-Xmx8g -jar ~/tools/snpEff/snpEff.jar og.wk21 "
vcf_fn <- "~/hdd2/og_gwas/genetics/input/variants/ogCore_variants.vcf"
out_fn <- "~/hdd2/og_gwas/genetics/input/variants/ogCore_variants.snpEff.vcf"
system2("java", paste0(args, vcf_fn, " > ", out_fn))

###########################
# List up candidate genes #
###########################
annotation_fn <- "~/hdd2/og_gwas/genomics/05_filesForDatabase/combinedGeneAnnotation.tsv"
gff_fn <- "~/hdd2/og_gwas/genomics/05_filesForDatabase/wk21_genome_chr.gff3"
snpeff_fn <- "~/hdd2/og_gwas/genetics/input/variants/ogCore_variants.snpEff.vcf"
out_fn <- file.path(out_dir, "test_")
listCandidate(peakblock,
              annotation_fn = annotation_fn,
              gff_fn = gff_fn,
              snpeff_fn = snpeff_fn,
              out_fn = out_fn)
