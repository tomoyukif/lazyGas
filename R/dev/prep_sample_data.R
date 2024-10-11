library(genovisr)
library(SNPRelate)
library(SeqArray)
library(Biostrings)

nb <- readDNAStringSet("~/hdd3/genomeData/rice/cultivar_sativa/os_nb_rapdb/IRGSP-1.0_genome.fasta")

df <- simGenovis(n_sample = 100, n_chr = 12, chr_len = width(nb) * 1e-6)

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
# 696    1_29869384   1 29869384          G          A
# 5334    6_9583957   6  9583957          G          A
# 11510 12_14220188  12 14220188          G          A

qtl_ds <- df$dosage[, qtl_index]

e <- rnorm(n = nrow(qtl_ds), mean = 0, sd = 1)
q1_add <- qtl_ds[, 1] * 3
q1_dom <- as.numeric(qtl_ds[, 1] == 1) * 3
q2_add <- qtl_ds[, 2] * 2
q2_dom <- as.numeric(qtl_ds[, 2] == 1) * 0
q3_add <- qtl_ds[, 3] * 1
q3_dom <- as.numeric(qtl_ds[, 3] == 1) * 2

st <- 1^2 + 3^2 + 2^2 + 1^2
3^2/st
2^2/st
1^2/st

pheno <- q1_add + q1_dom + q2_add + q2_dom + q3_add + q3_dom + e
pheno <- data.frame(id = df$sample_info$id, pheno = pheno)
write.csv(pheno, "inst/extdata/pheno.csv", row.names = FALSE)




