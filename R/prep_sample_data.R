library(genovisr)
library(SNPRelate)
library(SeqArray)

df <- simGenovis(n_sample = 100)

snpgdsCreateGeno(gds.fn = "inst/extdata/temp.gds",
                 genmat = df$genotype,
                 sample.id = df$sample_info$id,
                 snp.id = df$marker_info$id,
                 snp.rs.id = df$marker_info$id,
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
