library(genovisr)
library(SNPRelate)
library(SeqArray)

df <- simGenovis(n_sample = 100)

snpgdsCreateGeno(gds.fn = "inst/extdata/temp.gds",
                 genmat = df$genotype,
                 sample.id = df$sample_info$id,
                 snp.id = seq_along(df$marker_info$id),
                 snp.rs.id = seq_along(df$marker_info$id),
                 snp.chromosome = df$marker_info$chr,
                 snp.position = df$marker_info$pos,
                 snp.allele = paste(df$marker_info$ref_allele,
                                    df$marker_info$alt_allele, sep = "/"),
                 snpfirstdim = FALSE, compress.annotation = "LZMA_RA",
                 compress.geno = "LZMA_RA",
                 other.vars = list(HAP = df$haplotype,
                                   EDS = df$dosage))

seqSNP2GDS(gds.fn = "inst/extdata/temp.gds",
           out.fn = "inst/extdata/sample.gds")

gds <- openfn.gds("inst/extdata/sample.gds", readonly = FALSE)

addfolder.gdsn(node = index.gdsn(gds, "annotation/format"), name = "HAP")
addfolder.gdsn(node = index.gdsn(gds, "annotation/format"), name = "EDS")
add.gdsn(index.gdsn(gds, "annotation/format/HAP"), name = "data",
         val = df$haplotype, storage = "bit2", compress = "LZMA_RA")
add.gdsn(index.gdsn(gds, "annotation/format/EDS"), name = "data",
         val = df$dosage, storage = "bit2", compress = "LZMA_RA")
closefn.gds(gds)
unlink("inst/extdata/temp.gds")


qtl <- sample(x = df$marker_info$id, size = 3)
qtl_index <- which(df$marker_info$id %in% qtl)
df$marker_info[qtl_index, ]
# id chr      pos ref_allele alt_allele
# 3456    4_9644821   4  9644821          G          A
# 3871   4_19751697   4 19751697          G          A
# 11666 12_21496434  12 21496434          G          A

qtl_ds <- df$dosage[, qtl_index]

e <- rnorm(n = nrow(qtl_ds), mean = 0, sd = 1)
q1_add <- qtl_ds[, 1] * 5
q1_dom <- as.numeric(qtl_ds[, 1] == 1) * 5
q2_add <- qtl_ds[, 2] * 3
q2_dom <- as.numeric(qtl_ds[, 2] == 1) * 0
q3_add <- qtl_ds[, 3] * 2
q3_dom <- as.numeric(qtl_ds[, 3] == 1) * 1

pheno <- q1_add + q1_dom + q2_add + q2_dom + q3_add + q3_dom + e
pheno <- data.frame(id = df$sample_info$id, pheno = pheno)
write.csv(pheno, "inst/extdata/pheno.csv", row.names = FALSE)
