library(GBScleanR)
#browseVignettes(package = "GBScleanR")
#closeGDS(gds)


gds <- loadGDS("/Users/okapi//Desktop/GBSCreenR_V2/2022/TP72/disallow_missing/disallow_missing_genotype.gds",
               load_filter = T)

#dosage0 <- getGenotype(gds, node = "dosage")
#er <-getErrorRate(gds)

# gds <- makeScheme(gds, generation = 2, crosstype = "self")
# gds <- makeScheme(object = gds, generation = 5, crosstype = "self")

gds <- initScheme(gds, mating = rbind(1, 2))
gds <- addScheme(gds, crosstype = "selfing")

gds <- estGeno(gds)
er <-getErrorRate(gds)

#dosage1 <- getGenotype(gds, node = "dosage")

gds <- setPloidy(gds, ploidy = 3)


#sum(table(dosage1))

gds <- setFixedParameter(gds, bias = er$bias, er$mismap, parent_geno = T)
gds <- initScheme(gds, mating = rbind(1, 2))
gds <- addScheme(gds, crosstype = "selfing")
gds <- estGeno(gds)


dosage <-getGenotype(gds, "dosage")

table(dosage0)

cor <- getGenotype(gds, "cor", phased = T)
hap <- getHaplotype(gds)
hap[, 1, 1:20]




ID <- read.csv("/Users/okapi/Desktop/緊急事態_再解析/Book2.csv")

sample.id <- getSamID(object = gds)
snp.chromosome <- getChromosome(object = gds)
snp.position <- getPosition(object = gds)

sample.id <- paste0("sample_", formatC(as.numeric(sub("sample", "",sample.id)), width=3,flag="0"))
rownames(dosage) <- sample.id
sample.id <- sample.id[order(sample.id)]
dosage <- dosage[order(rownames(dosage)),]

rownames(dosage) <- paste0("sample_",sprintf("%03d",ID$New_ID))
dosage <- dosage[order(rownames(dosage)),]
sample.id <- paste0("sample_", sprintf("%03d",ID$New_ID))
sample.id <- sample.id[order(sample.id)]

setwd("/Users/okapi/Desktop/GBSCreenR_V2/genotype_for_anueploidy/TP72/3somy/")
length(snp.position)

colnames(dosage) <- paste0(snp.chromosome, "_", snp.position)

write.csv(dosage, "TP72genotype_ploidy3.csv")



pdf("TP72_3somyGenotype.pdf")
for (i in 1:length(rownames(dosage))) {
  p <- plotDosage(gds, coord = c(4, 3),ind = i)
  print(p)
}
dev.off()


