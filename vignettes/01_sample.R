#######################################################
# Build an object holding genotype and phenotype data #
#######################################################
geno_fn <- "~/hdd3/gbs/runs/2022_TP72_TP74_NBWRC46_KoshiWCSL/NippWRC46/gbscleanr_out/NBxWRC46_F2_MCPtaggR_GBScleanR_dosage.csv"
geno <- read.csv(geno_fn, header = FALSE)

str(geno)    # Check data structure
geno[1:5, 1:5]    # Check data structure

sampleID_geno <- geno[-(1:2), 1]    # Extract sample id info from the genotype data
marker_chr <- as.character(geno[1, -1])    # Extract marker position info (chromosome)
marker_pos <- as.numeric(geno[2, -1])    # Extract marker position info (basepairs)
geno <- geno[-(1:2), -1]    # Extract genotype info only

pheno_fn <- "~/hdd2/toriba/toriba_f2/pheno/Furuta-san_WRC46xNB_F2_Toriba211227.csv"
pheno <- read.csv(pheno_fn)
str(pheno)    # Check data structure
pheno[1:5, 1:5]    # Check data structure

sampleID_pheno <- pheno$ID    # Extract sample id info from the phenotype data
sampleID_pheno <- paste0("sample", sampleID_pheno)    # Modify format if neccessary
pheno <- pheno[, -1]    # Extract phenotype info

obj <- buildLazyQTL(geno = geno,
                    sampleID_geno = sampleID_geno,
                    sampleID_pheno = sampleID_pheno,
                    pheno = pheno,
                    marker_chr = marker_chr,
                    marker_pos = marker_pos)
obj    # Confirm the object structure

#######################
# Regression analysis #
#######################
out_dir <- "~/hdd2/softDevel/lazyGeno/test"    # Prepare a direcoty to output files
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Make a function to make a data.frame to be used for the regression
# We will use the model y ~  add + dom + e, where y is phenotype, add denotes
# the additive effect, dom represents the dominant effect, and e is the residual
# error term.
make_df <- function(x){
  add <- x    # Set x as additive effect
  dom <- as.numeric(!x %in% c(0, 2))    # Convert x into values indicating dominant effect
  df <- data.frame(add = add, dom = dom)
  return(df)
}

# Execute regression analysis
out_fn <- file.path(out_dir, "test_")
scanqtl <- scanQTL(obj,
                   out_fn = out_fn,
                   formula = "add + dom",
                   makeDF_FUN = make_df)

########################
# Draw manhattan plots #
########################
signif <- "x$FDR <= 0.05"
out_fn <- file.path(out_dir, "test_")
plotManhattan(scanqtl, signif = signif, out_fn = out_fn)


##############
# Call peaks #
##############
signif <- "x$FDR <= 0.05"
rsquare <- 0.6
out_fn <- file.path(out_dir, "test_")
peakblock <- callPeakBlock(scanqtl, signif = signif,
                           out_fn = out_fn, rsquare = rsquare)

##################
# Execute snpEff #
##################
args <- "-Xmx8g -jar ~/tools/snpEff/snpEff.jar nb_rapmsu"
vcf_fn <- "~/hdd3/gbs/runs/2022_TP72_TP74_NBWRC46_KoshiWCSL/NippWRC46/gbscleanr_out/NBxWRC46_F2_MCPtaggR_GBScleanR.vcf"
out_fn <- "~/hdd3/gbs/runs/2022_TP72_TP74_NBWRC46_KoshiWCSL/NippWRC46/gbscleanr_out/NBxWRC46_F2_MCPtaggR_GBScleanR.snpEff.vcf"
system2("java", paste0(args, vcf_fn, " > ", out_fn))

###########################
# List up candidate genes #
###########################
load("~/hdd3/genomeData/rice/cultivar_sativa/nb_combined/nbCombined_genetable.Rdata")
snpeff_fn <- "~/hdd3/gbs/runs/2022_TP72_TP74_NBWRC46_KoshiWCSL/NippWRC46/gbscleanr_out/NBxWRC46_F2_MCPtaggR_GBScleanR.snpEff.vcf"
out_fn <- file.path(out_dir, "test_")
ann <- subset(ann, select = c(GeneID, rep_TxID, chr:end))
names(ann) <- c("GeneID", "TxID", "Chr", "Start", "End")
listCandidate(peakblock,
              annotation_fn = ann,
              snpeff_fn = snpeff_fn,
              out_fn = out_fn)

###################################
# Plot haplotype x phenotype plot #
###################################
out_dir <- "~/hdd2/softDevel/lazyGeno/test"
pvalues_fn <- list.files(out_dir, "scanQTL", full.names = TRUE)

geno_fn <- "~/hdd3/gbs/runs/2022_TP72_TP74_NBWRC46_KoshiWCSL/NippWRC46/gbscleanr_out/NBxWRC46_F2_MCPtaggR_GBScleanR_dosage.csv"
geno <- read.csv(geno_fn, header = FALSE)
sampleID_geno <- geno[-(1:2), 1]    # Extract sample id info from the genotype data
marker_chr <- as.character(geno[1, -1])    # Extract marker position info (chromosome)
marker_pos <- as.numeric(geno[2, -1])    # Extract marker position info (basepairs)
geno <- geno[-(1:2), -1]    # Extract genotype info only
pheno_fn <- "~/hdd2/toriba/toriba_f2/pheno/Furuta-san_WRC46xNB_F2_Toriba211227.csv"
pheno <- read.csv(pheno_fn)
sampleID_pheno <- pheno$ID    # Extract sample id info from the phenotype data
sampleID_pheno <- paste0("sample", sampleID_pheno)    # Modify format if neccessary
pheno <- pheno[, -1]    # Extract phenotype info

qtlscan <- buildQTLscan(geno = geno,
                        pvalues_fn = pvalues_fn,
                        sampleID_geno = sampleID_geno,
                        sampleID_pheno = sampleID_pheno,
                        pheno = pheno,
                        marker_chr = marker_chr,
                        marker_pos = marker_pos)

signif <- "x$FDR <= 0.05"
rsquare <- 0.6
out_fn <- file.path(out_dir, "test_")
peakblock <- callPeakBlock(qtlscan, signif = signif,
                           out_fn = out_fn, rsquare = rsquare)

haplo_fn <- haploPlot(peakblock, out_fn = out_fn, out_fmt = "pdf")
