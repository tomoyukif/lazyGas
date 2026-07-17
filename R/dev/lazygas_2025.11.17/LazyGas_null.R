# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
devtools::install_github("tomoyukif/lazyGas", build_vignettes = F)

library(lazyGas)
library(GBScleanR)
library(ggplot2)
library(GenovisR)
library(SNPRelate)
library(SeqArray)
library(Biostrings)
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(rtracklayer)
# library(CMplot)
library(dplyr)



ploidy = 4
pad_and_rbind <- function(dfs) {

  max_cols <- max(sapply(dfs, ncol))
  col_names <- as.character(1:max_cols)
  dfs_aligned <- lapply(dfs, function(df) {
    nc <- ncol(df)
    colnames(df) <- as.character(1:nc)
    if (nc < max_cols) {
      missing_cols <- setdiff(col_names, colnames(df))
      for (col in missing_cols) {
        df[[col]] <- ""
      }
      df <- df[ , col_names]
    }
    return(df)
  })

  do.call(rbind, dfs_aligned)
}

gff_fn <- "R/dev/lazygas_2025.11.17/nb_combined_all.gff"
gff <- rtracklayer::import.gff(gff_fn)

ann <- read.csv("R/dev/lazygas_2025.11.17/orthopair_summary_annotated.csv")
msa <- read.csv("R/dev/lazygas_2025.11.17/orthopair_summary_msa.csv")
ann <- subset(ann, select = c(NB, WK21, NB_original_gene_id, NB_chr, NB_start, WK21_chr, WK21_start,
                              NB_AA_len:WK21_TE, MSU_Note:WK21_PANTHER, NB_Pfam, WK21_Pfam))
msa <- subset(msa, select = c(NB, WK21,  NB_WK21_AA_Del_Length:NB_WK21_AA_Mis_Number))
ann <- left_join(x = ann, y = msa, by = c("NB", "WK21"))
names(ann)[3] <- "Gene_ID"

rap <- import.gff3("R/dev/lazygas_2025.11.17/locus.gff")
rap <- mcols(rap)
names(rap) <- c("source", "type", "score", "phase", "Gene_ID", "Name", "RAP_Note",
                "TranscriptVariants", "Oryzabase_geneSymbol",
                "Oryzabase_geneName", "RAPDB_geneSymbol",
                "RAPDB_geneName", "CGSNL_geneSymbol", "CGSNL_geneName")
rap <- subset(rap, select = c(Gene_ID, RAP_Note, Oryzabase_geneName,
                              RAPDB_geneName, CGSNL_geneName))
rap <- as.data.frame(rap)
ann <- left_join(x = ann, y = rap, by = "Gene_ID")
ann <- subset(ann, select = c(NB:WK21_TE, RAP_Note, MSU_Note:CGSNL_geneSymbol,
                              Oryzabase_geneName:CGSNL_geneName,
                              Oryzabase_Trait_Ontology:NB_WK21_AA_Mis_Number))

anuelist <- read.csv("R/dev/lazygas_2025.11.17/anueploidlist.csv",
                     row.names = 1)
anuelist <- anuelist[-c(60, 98, 191, 192),]
fdr_boarder <- "0.10"
peak_list_allpheno <- NULL

ID <- read.csv("R/dev/lazygas_2025.11.17/Book2.csv")

phe_data <- read.csv("R/dev/lazygas_2025.11.17/Phenotype_緊急事態.csv")
phe_72 <- phe_data[phe_data$Population_ID == "21TP72",]
phe_72$Serial_N <- formatC(phe_72$Serial_N,width=3,flag="0")
phe_72 <- phe_72[-c(98),]
phe_72 <- phe_72[, -c(12, 13, 17)]
phe_vec <- c("Days to Heading", "Shattering", "Awn Length (mm)", "Culm Length (cm)", "Panicle Length (cm)", "Panicle Number",
             "Total Grain Number", "Spikelet Namber per Panicle", "Seed Fertility (%)","Grain Area Size (mm2)", "Grain Length (mm)", "Grain Width (mm)")
phe_vec_2 <- c("dth", "sh", "al", "cl", "pl", "pn", "tgn", "sn", "sf", "ga", "gl", "gw")



dosage_3somy <- read.csv("R/dev/lazygas_2025.11.17/TP72genotype_ploidy3.csv",
                         row.names = 1)
dosage_3somy <- dosage_3somy * 20
dosage_3somy <- dosage_3somy[-c(60, 98, 191, 192),]
colnames(dosage_3somy) <- 1:length(colnames(dosage_3somy))
dosage_3somy <- as.matrix(dosage_3somy)

dosage_5somy <- read.csv("R/dev/lazygas_2025.11.17/TP72genotype_ploidy5.csv",
                         row.names = 1)
dosage_5somy <- dosage_5somy * 12
dosage_5somy <- dosage_5somy[-c(60, 98, 191, 192),]
colnames(dosage_5somy) <- 1:length(colnames(dosage_5somy))
dosage_5somy <- as.matrix(dosage_5somy)


gds <- loadGDS("R/dev/lazygas_2025.11.17/genotype_disallow_missing.gds",
               load_filter = T, ploidy = 4)
#dosage2 <- getGenotype(object = gds, node = "dosage")

genotype <- getGenotype(object = gds, node = "raw")
dosage <- getGenotype(object = gds, node = "dosage")
haplotype <- getHaplotype(gds)
snp.chromosome <- getChromosome(object = gds)
snp.position <- getPosition(object = gds)
snp.allele <- getAllele(object = gds)
sample.id <- getSamID(object = gds)
snp.id <- getMarID(object = gds)
#dosage[is.na(dosage)] <- 63


#sample.idとdosage
sample.id <- paste0("sample_", formatC(as.numeric(sub("sample", "",sample.id)), width=3,flag="0"))
rownames(dosage) <- sample.id
sample.id <- sample.id[order(sample.id)]
dosage <- dosage[order(rownames(dosage)),]

#TP72
rownames(dosage) <- paste0("sample_",sprintf("%03d",ID$New_ID))
dosage <- dosage[order(rownames(dosage)),]
sample.id <- paste0("sample_", sprintf("%03d",ID$New_ID))
sample.id <- sample.id[order(sample.id)]
sample.id <- sample.id[-c(60, 98, 191, 192)]
dosage <- dosage[-c(60, 98, 191, 192),]
dosage <- dosage * 15

#2025不要ならコメントアウト
####################################################################################
#2025の個体を追加
anuelist25 <- read.csv("R/dev/lazygas_2025.11.17/anuelist_2025.csv",
                       row.names = 1)
rownames(anuelist25) <- anuelist25$ind_name
anuelist25 <- anuelist25[18:32,-1]

dosage_3somy25 <- read.csv("R/dev/lazygas_2025.11.17/2025genotype_ploidy3.csv",
                           row.names = 1)
dosage_3somy25 <- dosage_3somy25 * 20
dosage_3somy25 <- dosage_3somy25[18:32,]
dosage_3somy25 <- dosage_3somy25[order(rownames(dosage_3somy25)),]
colnames(dosage_3somy25) <- 1:length(colnames(dosage_3somy25))
dosage_3somy25 <- as.matrix(dosage_3somy25)

dosage_5somy25 <- read.csv("R/dev/lazygas_2025.11.17/2025genotype_ploidy5.csv",
                           row.names = 1)
dosage_5somy25 <- dosage_5somy25 * 12
dosage_5somy25 <- dosage_5somy25[18:32,]
dosage_5somy25 <- dosage_5somy25[order(rownames(dosage_5somy25)),]
colnames(dosage_5somy25) <- 1:length(colnames(dosage_5somy25))
dosage_5somy25 <- as.matrix(dosage_5somy25)


dosage2025 <- read.csv("R/dev/lazygas_2025.11.17/2025genotype_okaestimatetetraploid.csv",
                       row.names = 1)
adddosage <- dosage2025[18:32,]
adddosage <- adddosage[order(rownames(adddosage)),]
adddosage <- adddosage * 15

addmarker <- colnames(adddosage)
addchr <- sub("_.*", "", addmarker)


for (i in 1:length(rownames(anuelist25))) {
  if(sum(!is.na(anuelist25[i,])) == 0) {
  } else {
    for (j in 1:length(colnames(anuelist25))) {
      anue_j <- anuelist25[i, j]
      if (is.na(anue_j) == 0) {
        if (anue_j == 3) {
          chr_j <- which(addchr == paste0("chr", sprintf("%02d", j)))
          adddosage[i, chr_j] <- dosage_3somy25[i, chr_j]
        } else if (anue_j == 5) {
          chr_j <- which(addchr == paste0("chr", sprintf("%02d", j)))
          adddosage[i, chr_j] <- dosage_5somy25[i, chr_j]
        } else {
          stop("stop!!!")
        }
      } else {
      }
    }
  }
}

marker2025 <- colnames(adddosage)
marker2025 <- data.frame(snp = marker2025,
                         chr = sub("_.*", "", marker2025),
                         pos = as.numeric(sub(".*_", "", marker2025)))
marker72 <- data.frame(snp = paste0(snp.chromosome, "_", snp.position),
                       chr = snp.chromosome,
                       pos = snp.position)

i <- 1
SNPref <- NULL
for (i in 1:12) {
  snp25_i <- filter(marker2025, chr == paste0("chr", sprintf("%02d", i)))
  snp72_i <- filter(marker72, chr == paste0("chr", sprintf("%02d", i)))

  SNPref_j <- NULL
  for (j in 1:length(rownames(snp72_i))) {
    ref_j <-  which.min(abs(snp25_i[, "pos"] - snp72_i[, "pos"][j]))

    if(is.vector(ref_j) == T){
      ref_j <- ref_j[1]
    } else{
      ref_j <- ref_j
    }
    ref_j <- snp25_i$pos[ref_j]
    SNPref_j <- c(SNPref_j, ref_j)
  }
  length(SNPref_j)
  SNPref_j <- cbind(paste0("chr", sprintf("%02d", i)),
                    SNPref_j)
  SNPref <- rbind(SNPref, SNPref_j)
}
length(rownames(SNPref))
SNPref <- as.data.frame(SNPref)
SNPref$SNP <- paste0(SNPref$V1, "_", SNPref$SNPref_j)
length(unique(SNPref$SNP))

compSNP <- data.frame(newsnp = marker72$snp,
                      oldsnp = SNPref$SNP)

adddosage <- as.data.frame(t(adddosage))
adddosage$oldsnp <- rownames(adddosage)

adddosage <- left_join(x = compSNP, y = adddosage, by = "oldsnp")
rownames(adddosage) <- adddosage$newsnp
adddosage <- adddosage[,-c(1, 2)]
adddosage <- as.data.frame(t(adddosage))
rownames(adddosage) <- paste0("sample_", sub("21TP072_", "", rownames(adddosage)))
#adddosageファイル完成
#ここまでコメントアウト
####################################################################################


for (i in 1:length(rownames(anuelist))) {
  if(sum(!is.na(anuelist[i,])) == 0) {
  } else {
    for (j in 1:length(colnames(anuelist))) {
      anue_j <- anuelist[i, j]
      if (is.na(anue_j) == 0) {
        if (anue_j == 3) {
          chr_j <- which(snp.chromosome == paste0("chr", sprintf("%02d", j)))
          dosage[i, chr_j] <- dosage_3somy[i, chr_j]
        } else if (anue_j == 5) {
          chr_j <- which(snp.chromosome == paste0("chr", sprintf("%02d", j)))
          dosage[i, chr_j] <- dosage_5somy[i, chr_j]
        } else {
          stop("stop!!!")
        }
      } else {
      }
    }
  }
}

####################################################################################
#2025dataに置換
dosage <- dosage[-which(rownames(dosage) %in% rownames(adddosage)),]
colnames(adddosage) <- colnames(dosage)
dosage <- rbind(dosage, adddosage)
dosage <- dosage[order(rownames(dosage)),]
#置換完了
#ここまでコメントアウト
####################################################################################


#歪み除去のため組み替え近傍マーカー除こう

geno_ancestor <- read.csv("R/dev/lazygas_2025.11.17/genotype.csv",
                          row.names = 1) # Your file
anc_marker <- rownames(geno_ancestor)
anc_ds <- geno_ancestor[,"X20TP41.6"]
anc_SNP <- data.frame(marker = anc_marker,
                      chr = sub("_.*", "", anc_marker),
                      pos = as.numeric(sub(".*_", "", anc_marker)),
                      dos = anc_ds)

prog_SNP <- data.frame(marker = paste0(snp.chromosome, "_", snp.position),
                       chr = snp.chromosome,
                       pos = snp.position)

#データの前処理
#####################################
#TP102集団のSNP近傍の親個体のSNPを検索
chr_order <- paste0("chr", as.character(formatC(1:12, width = 2, flag = "0")))
SNPref <- NULL
for (i in 1:12) {
  anc_i <- filter(anc_SNP, chr == chr_order[i])
  prog_i <- filter(prog_SNP, chr == chr_order[i])

  SNPref_j <- NULL
  for (j in 1:length(rownames(prog_i))) {
    ref_j <-  which.min(abs(anc_i[, "pos"] - prog_i[, "pos"][j]))

    if(is.vector(ref_j) == T){
      ref_j <- ref_j[1]
    } else{
      ref_j <- ref_j
    }
    ref_j <- anc_i$pos[ref_j]
    SNPref_j <- c(SNPref_j, ref_j)
  }
  SNPref <- c(SNPref, SNPref_j)
}
prog_SNP <- cbind(prog_SNP, SNPref)
prog_SNP$refmarker <- paste0(prog_SNP$chr, "_", as.character(SNPref))
#近傍のSNP検索終了
#####################################

#####################################
#近傍マーカーのdosage検索
ref_dos <- NULL
for (h in 1:length(rownames(prog_SNP))) {
  marker_h <- prog_SNP$refmarker[h]
  ref_h <- anc_SNP[anc_SNP$marker == marker_h, "dos"]
  ref_dos <- c(ref_dos, ref_h)
}
prog_SNP$refdosage <- ref_dos


#breaking pointの検索
##############################
df_forBP <- prog_SNP
#近傍のbreaking pointを推定、距離を算出
BPvec <- NULL
for (i in 1:12) {
  chr_i <- filter(df_forBP, chr == chr_order[i])
  phasing_BP <- which((chr_i$refdosage[-length(rownames(chr_i))] - chr_i$refdosage[-1])　!= 0)
  breakingpoint <- (chr_i[phasing_BP, ]$pos + chr_i[(phasing_BP + 1), ]$pos)/2
  breakingpoint <- c(chr_i$pos[1], breakingpoint, chr_i$pos[length(rownames(chr_i))])
  for (j in 1:length(rownames(chr_i))) {
    closeBP<- breakingpoint[which.min(abs(chr_i$pos[j] - breakingpoint))]
    BPvec <- c(BPvec, closeBP)
  }
}
df_forBP <- cbind(df_forBP, BPvec)
df_forBP$BPdis <- abs(df_forBP$pos - df_forBP$BPvec)
#breaking pointの推定と計算完了
##############################
#breaking pointの検索
################################################################################

##############################
#除外領域の指定
deldis <- 400000
del <- paste0(as.character(deldis/1000), "kb")
removemarker <- which(df_forBP$BPdis <= deldis)
#除外領域の指定完了
##############################
#諸々完了
################################################################################

#組み替え点近傍除外
genotype <- genotype[, -removemarker]
snp.id <- snp.id[-removemarker]
snp.chromosome <- snp.chromosome[-removemarker]
snp.position <- snp.position[-removemarker]
snp.allele <-  snp.allele[-removemarker]
haplotype <- haplotype[,, -removemarker]
dosage <- dosage[, -removemarker]
#table((dosage))
#2025不要ならコメントアウト
##############################
sample.id <- rownames(dosage)
dosage <- as.matrix(dosage)
##############################
result <- apply(dosage, 2, function(x) {
  table(factor(x, levels = c(0, 12, 15, 20, 24, 30, 36, 40, 45, 48, 60)))
})

result <- as.data.frame(t(result))
rownames(result) <- paste0(snp.chromosome, "_", snp.position)
result <- result / sum(result[1,])
a <- apply(result, 1, max)
result$tf <- a > 0.9 & a < 1

dosage[, which(result$tf == T)] <- NA






aneuploid_ind <- apply(anuelist, 1, function(x){!(sum(na.omit(x)) == 0)})
#sample.id %in% rownames(anuelist)
aneuploid_ind <- c(aneuploid_ind, F)
names(aneuploid_ind)[189] <- "sample_060"
aneuploid_ind <- aneuploid_ind[order(names(aneuploid_ind))]
#which(aneuploid_ind == T)

# sample.id <- sample.id[-which(aneuploid_ind == T)]
# genotype <- genotype[-which(aneuploid_ind == T),]
# haplotype <- haplotype[,-which(aneuploid_ind == T), ]
# dosage <- dosage[-which(aneuploid_ind == T), ]




temp_dir <- tempdir()
create_gds <-  list(genotype = dosage,
                    sample.id = sample.id,
                    snp.id = snp.id,
                    snp.rs.id = NULL,
                    snp.chromosome = snp.chromosome,
                    snp.position = snp.position,
                    snp.allele = snp.allele,
                    haplotype = haplotype,
                    dosage = dosage)



sample_gds2 <- tempfile("sample", temp_dir, ".gds")
lg <- buildLazyGas(gds_fn = sample_gds2,
                   create_gds = create_gds)



# save.image(file = "R/dev/lg.rdata")

#
# create_gds$sample.id
#
# lg[[sample.id]]
#j <- 3
###############################################################################################################################
recalc_df <- NULL
peak_df <- NULL
peak_c_df <- NULL
for (j in 1:length(phe_vec)) {


  pheno <- phe_72[, j + 2]
  #pheno <- pheno[-which(aneuploid_ind == T)]
  pheno <- data.frame(sample.id, pheno)
  colnames(pheno) <- c("id", "pheno")

  lg <- assignPheno(object = lg,
                    pheno = pheno,
                    rename = phe_vec[j])

  pheno <- getPheno(object = lg)
  for(i in seq_along(pheno$pheno_names)){
    p <- plotPheno(object = lg, pheno = i, xlab = pheno$pheno_names[i])
    print(p)
  }


  g <- getGenoPerMarker(object = lg, geno_format = "dosage")
  print(g)

  g <- as.integer(dosage[,1])

  #add = ifelse(g == 0, -1, ifelse(g == 60, 1, 0))
  #add = ifelse(g == 15, -1, ifelse(g == 45, 1, ifelse(g == 20, -2/3, ifelse(g == 40, 2/3, ifelse(g == 12, -6/5, ifelse(g == 24, -2/5, ifelse(g == 36, 2/5, ifelse(g == 48, 6/5, 0))))))))
  #add = ifelse(g == 0, -2, ifelse(g == 12, -6/5, ifelse(g == 15, -1, ifelse(g == 20, -2/3, ifelse(g == 24, -2/5, ifelse(g == 30, 0, ifelse(g == 36, 2/5, ifelse(g == 40, 2/3, ifelse(g == 45, 1, ifelse(g == 48, 6/5, 0))))))))))
  #dom = ifelse(g == 15, -2, ifelse(g == 45, 2, ifelse(g == 20, -4/3, ifelse(g == 40, 4/3, ifelse(g == 12, -4/5, ifelse(g == 24, -4/5, ifelse(g == 36, 4/5, ifelse(g == 48, 12/5, 0))))))))
  #dom = ifelse(g == 12, -6/5, ifelse(g == 15, -1, ifelse(g == 20, -2/3, ifelse(g == 24, -2/5, ifelse(g == 36, 2/5, ifelse(g == 40, 2/3, ifelse(g == 45, 1, ifelse(g == 48, 6/5, 0))))))))
  conv_fun <- function(g) {
    add <- ifelse(g == 0, -1, ifelse(g == 60, 1, 0))
    dom <- as.numeric(g %in% c(12, 15, 20, 24, 30, 36, 40, 45, 48))
    dose <- ifelse(g == 12, -6/5, ifelse(g == 15, -1, ifelse(g == 20, -2/3, ifelse(g == 24, -2/5, ifelse(g == 36, 2/5, ifelse(g == 40, 2/3, ifelse(g == 45, 1, ifelse(g == 48, 6/5, 0))))))))
    out <- data.frame(add = add, dom = dom, dose = dose)
    return(out)
  }

  # For the model matrix above, the formula can be the following.
  formula <- "add + dom + dose + Aneuploid"

  Aneuploid <- aneuploid_ind

  null_formula <- "Aneuploid"


  scanAssoc(object = lg,
            formula = formula,
            null_formula = null_formula,
            fixed_effect = data.frame(Aneuploid = sample(c(0, 1), size = length(sample.id), replace = T)),
            conv_fun = conv_fun,
            geno_format = "dosage",
            kruskal = NULL)

  #ここでsignif
  pheno <- getPheno(object = lg)
  for(i in seq_along(pheno$pheno_names)){
    if (fdr_boarder == "0.05") {
      p <- plotManhattan(object = lg, pheno = i)
    } else {
      p <- plotManhattan(object = lg, pheno = i, signif = as.numeric(fdr_boarder))
    }
    p <- p + labs(title = pheno$pheno_names[i])
    print(p)
  }


  #行った気がする


  scan <- lazyData(object = lg, dataset = "scan", pheno = pheno$pheno_names[1])
  scan[which.max(scan$negLog10P), ]


  callPeakBlock(object = lg, signif = as.numeric(fdr_boarder), threshold = 0.8)

  for(i in seq_along(pheno$pheno_names)){
    p <- plotPeaks(object = lg, pheno = i, recalc = FALSE)
    p <- p + labs(title = pheno$pheno_names[i])
    print(p)
  }
  #変更注意
  #groupthrashみたいなんで重回帰の閾値
  recalcAssoc(object = lg, n_threads = 10, grouping_threshold = 0.05)
  #recalcAssoc(object = lg, n_threads = 10, grouping_threshold = 0.1)

  for(i in seq_along(pheno$pheno_names)){
    p <- plotPeaks(object = lg, pheno = i, recalc = TRUE)
    p <- p + labs(title = pheno$pheno_names[i])
    print(p)
  }

  ## List up candidate genes
  # gff_fn <- "~/Desktop/gdsdata_2025Apr/tp_GBSCleanR_2025Apr_fromFuruta/gff/nb_combined_all.gff"
  # gff <- rtracklayer::import.gff(gff_fn)
  #
  # ann <- read.csv("/Users/okapi/Desktop/gdsdata_2025Apr/tp_GBSCleanR_2025Apr_fromFuruta/ortholog/orthopair_summary_annotated.csv")
  # msa <- read.csv("/Users/okapi/Desktop/gdsdata_2025Apr/tp_GBSCleanR_2025Apr_fromFuruta/ortholog/orthopair_summary_msa.csv")
  # ann <- subset(ann, select = c(NB, WK21, NB_original_gene_id, NB_chr, NB_start, WK21_chr, WK21_start,
  #                               NB_AA_len:WK21_PANTHER, NB_Pfam, WK21_Pfam))
  # msa <- subset(msa, select = c(NB, WK21,  NB_WK21_AA_Del_Length:NB_WK21_AA_Mis_Number))
  # ann <- left_join(x = ann, y = msa, by = c("NB", "WK21"))
  # names(ann)[3] <- "Gene_ID"


  # listCandidate(object = lg, gff = gff,
  #             ann = ann, recalc = TRUE)

  # #TF注意
  # for(i in seq_along(pheno$pheno_names)){
  #   pdf(file = paste0(phe_vec_2[j], "_", fdr_boarder,"_boxplot.pdf"), onefile = T)
  #   out <- haploPlot(object = lg, pheno = i, recalc =F)
  #   print(out)
  #   dev.off()
  # }


  #############################################################################################################################
  # peakcall <- lazyData(object = lg, dataset = "peakcall", pheno = pheno$pheno_names[1])
  # write.csv(peakcall, paste0("allpeak_TP72_", phe_vec_2[j], "_", fdr_boarder, ".csv"))
  recalc <- lazyData(object = lg, dataset = "recalc", pheno = pheno$pheno_names[1])
  write.csv(recalc, paste0("TP72_", phe_vec_2[j], "_", fdr_boarder, ".csv"))
  # candidate <- lazyData(object = lg, dataset = "candidate", pheno = pheno$pheno_names[1])
  #write.csv(candidate, paste0("TP72_", phe_vec_2[j], "Candidate_", fdr_boarder, ".csv"))
  peakscan <- lazyData(object = lg, dataset = "scan", pheno = pheno$pheno_names[1])
  makeInteractiveSummary(object = lg, pheno = pheno$pheno_names[1],
                         out_fn = paste0("TP72_", phe_vec_2[j], "_", fdr_boarder, ".html"),
                         what = c("scan", "scan_png", "peakcall", "recalc", "groups", "preakcall_haplo", "recalc_haplo"))
  #############################################################################################################################








  if (is.data.frame(recalc) == T) {
    recalc <- left_join(recalc, peakscan, by = c("Chr", "Pos"))
    for (peak in 1:length(unique(recalc$peak_ID))) {
      peak_i <- filter(recalc, peak_ID == unique(recalc$peak_ID)[peak])
      peak_SNP <- filter(peak_i, dist2peak == 0)
      peak_SNP$QTL <- paste0("q72", phe_vec_2[j], sub("chr", "",peak_i$peak_Chr[1]), "_", sprintf("%02d",round(peak_i$peak_Pos[1] / 1000000)))
      peak_SNP$SNP <- paste0(peak_i$peak_Chr[1], "_", peak_i$peak_Pos[1])
      peak_SNP$start <- peak_i$Pos[1]
      peak_SNP$end <- peak_i$Pos[length(rownames(peak_i))]
      peak_SNP$pheno <- phe_vec_2[j]
      peak_df <- rbind(peak_df, peak_SNP)
    }
  }




  # peakscan <- lazyData(object = lg, dataset = "scan", pheno = pheno$pheno_names[1])
  # if (is.data.frame(peakcall) == T) {
  #   for (peak in 1:peakcall$peak_ID[length(rownames(peakcall))]) {
  #     peak_c_i <- filter(peakcall, peak_ID == peak)
  #     peak_c_SNP <- filter(peak_c_i, dist2peak == 0)
  #     peak_c_SNP$QTL <- paste0("q72", phe_vec_2[j], sub("chr", "",peak_c_i$peak_Chr[1]), "_", sprintf("%05d",floor(peak_c_i$peak_Pos[1] / 1000)))
  #     peak_c_SNP$SNP <- paste0(peak_c_i$peak_Chr[1], "_", peak_c_i$peak_Pos[1])
  #     peak_c_SNP$start <- peak_c_i$Pos[1]
  #     peak_c_SNP$end <- peak_c_i$Pos[length(rownames(peak_c_i))]
  #     peak_c_SNP$pheno <- pheno$pheno_names[1]
  #     scan_i <- filter(peakscan, Chr == peak_c_SNP$peak_Chr, Pos == peak_c_SNP$peak_Pos)[, 4:11]
  #     peak_c_SNP <- cbind(peak_c_SNP, scan_i)
  #     peak_c_df <- rbind(peak_c_df, peak_c_SNP)
  #   }
  # }

}

QTL_list_recalc <- data.frame(Phenotype = peak_df$pheno,
                              Population = rep("21TP72", length(rownames(peak_df))),
                              QTL = peak_df$QTL,
                              Chr = peak_df$peak_Chr,
                              peak = peak_df$SNP,
                              Start = peak_df$start,
                              End = peak_df$end,
                              P.model = peak_df$P.model,
                              P.add = peak_df$P.add,
                              P.dom = peak_df$P.dom,
                              P.dose = peak_df$P.dose,
                              Coef.add = peak_df$Coef.add,
                              Coef.dom = peak_df$Coef.dom,
                              Coef.dose = peak_df$Coef.dose,
                              PEV = peak_df$PVE)

yobi <- QTL_list_recalc

write.csv(QTL_list_recalc, paste0("QTLlit_TP72_recalc_", fdr_boarder,".csv"))


#2025不要ならコメントアウト
##############################
anuelist_new <- anuelist
rownames(anuelist25) <- sub("21TP072_", "sample_", rownames(anuelist25))
anuelist_new <- anuelist_new[-which(rownames(anuelist) %in% rownames(anuelist25)),]
anuelist_new <- rbind(anuelist_new, anuelist25)
anuelist_new <- anuelist_new[order(rownames(anuelist_new)),]
anuelist <- anuelist_new
##############################

# outdir <- "./output" # Your directory
# dir.create(outdir, showWarnings = FALSE)

setwd("/Users/okapi/Desktop/LazyGas物置/")
outdir <- "./boxplot" # Your directory
dir.create(outdir, showWarnings = FALSE)
phename_df <- cbind(colnames(phe_72)[-c(1, 2)], phe_vec, phe_vec_2)
for (i in 1:length(rownames(QTL_list_recalc))) {
  box_df <- cbind(phe_72[,phename_df[which(phename_df[,3] == QTL_list_recalc$Phenotype[i]), 1]],
                  dosage[ ,which(paste0(snp.chromosome, "_", snp.position) == (QTL_list_recalc$peak[i]))])
  colnames(box_df) <- c("pheno", "plex")
  box_df <- as.data.frame(box_df)

  anue_ind <- rownames(anuelist)[which(!is.na(anuelist[,QTL_list_recalc$Chr[i]]))]
  box_df<- box_df[-which(rownames(box_df) %in% anue_ind),]

  box_df$plex <- as.numeric(box_df$plex / 15)
  box_df <- rbind(box_df, cbind(plex = 0:4, pheno = NA))
  box_df$plex <- factor(box_df$plex, levels=c(0:4))
  p <- ggplot(box_df,aes(x = plex, y = pheno,
                         fill = plex)) +
    geom_boxplot(width=0.8,
                 alpha=0.7,
                 lwd=2) +
    xlab(NULL) +
    scale_fill_manual(breaks = c(0:4),
                      values = c("#fde725", "#5ec962", "#21918c", "#3b528b", "#440154")) +
    ylab(phename_df[which(phename_df[,3] == QTL_list_recalc$Phenotype[i]), 2]) +
    xlab("plex") +
    guides(fill="none") +
    geom_jitter(color="black", size=2, alpha=0.9, width = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 7.5, fill = "red") +
    labs(title = QTL_list_recalc$QTL[i]) +
    theme(axis.text.x = element_text(size = 25, color = "black"),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          title = element_text(size = 20),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 5))

  print(p)
  ggsave(filename = file.path(outdir,
                              paste0("21TP72", QTL_list_recalc$Phenotype[i], "_", QTL_list_recalc$QTL[i], "_", "TP72", "_boxplot.pdf")),
         plot = p,
         width = 7.5, height = 7.5)
}



setwd("/Users/okapi/Desktop/LazyGas物置/")
keywords_list <- read.csv("/Users/okapi/Desktop/gdsdata_2025Apr/candidate_keywords.csv")
outdir <- "./candidate" # Your directory
dir.create(outdir, showWarnings = FALSE)
ann_df <- as.data.frame(ann)
ann_df$gene_number <- 1:length(rownames(ann_df))
for (i in unique(QTL_list_recalc$Phenotype)) {
  candidate_list <- NULL
  QTL_i <- filter(QTL_list_recalc, Phenotype == i)
  for (j in 1:length(rownames(QTL_i))) {
    QTL_j <- QTL_i[j,]
    start_j <- QTL_j$Start - 400000
    end_j <- QTL_j$End + 400000
    candidate_NB <- filter(ann_df, NB_chr == QTL_j$Chr)
    candidate_NB <- filter(candidate_NB, NB_start >= start_j)
    candidate_NB <- filter(candidate_NB, NB_start <= end_j)
    WKnumber <- filter(candidate_NB, WK21_chr == QTL_j$Chr)
    start_WK <- WKnumber[1,]
    end_WK <- WKnumber[length(WKnumber),]
    candidate_WK <- filter(ann_df, WK21_chr == QTL_j$Chr)
    candidate_WK <- filter(candidate_WK, gene_number >= start_WK$gene_number)
    candidate_WK <- filter(candidate_WK, gene_number <= end_WK$gene_number)
    candidate_number <- unique(c(candidate_NB$gene_number, candidate_WK$gene_number))
    candidate_j <- filter(ann_df, gene_number %in% candidate_number)
    candidate_j$QTL <- QTL_j$QTL
    candidate_j$population <- QTL_j$Population
    candidate_list <- rbind(candidate_list, candidate_j)
  }

  key_j <- filter(keywords_list, phenotype == i)
  key_j <- as.vector(key_j)[-c(1, 2)]

  if (length(key_j[-which(key_j %in% "")]) > 0) {
    key_j <- key_j[-which(key_j %in% "")]
  }

  key_cand <- NULL
  for (k in 1:length(key_j)) {
    keywords <- key_j[k]

    key_cand_1 <- candidate_list[grepl(keywords, candidate_list$RAP_Note,
                                       ignore.case = T), ]
    key_cand_2 <- candidate_list[grepl(keywords, candidate_list$RAPDB_geneName,
                                       ignore.case = T), ]
    key_cand_3 <- candidate_list[grepl(keywords, candidate_list$Oryzabase_geneName,
                                       ignore.case = T), ]
    key_cand_4 <- candidate_list[grepl(keywords, candidate_list$CGSNL_geneName,
                                       ignore.case = T), ]
    key_cand_k <- rbind(key_cand_1, key_cand_2, key_cand_3, key_cand_4)

    key_cand <- rbind(key_cand, key_cand_k)


    # #ここで一回output
    # write.csv(key_cand,
    #           file = file.path(outdir, paste0("TP72_",i, "_CandList.csv")))
  }
  key_cand_unique <- unique(key_cand)


  cand_k <- list()
  for (k in unique(key_cand_unique$QTL)) {
    key_cand_k <- filter(key_cand_unique, QTL == k)
    cand_l <- NULL

    for (l in 1:length(rownames(key_cand_k))) {
      key_cand_l <- subset(key_cand_k, select = c(RAPDB_geneSymbol, CGSNL_geneSymbol, Gene_ID, NB_chr, NB_start, WK21_chr, WK21_start, NB_WK21_AA_Del_Length, NB_WK21_AA_Ins_Length, NB_WK21_AA_Mis_Number))
      colnames(key_cand_l) <- c("gene", "symbol", "ID", "NBchr", "NBpos", "WKchr", "WKpos", "del", "ins", "mis")
      key_cand_l <- key_cand_l[l,]
      candvec_l <- as.vector(key_cand_l)
      cand_l <- c(cand_l, key_cand_l)
    }
    cand_l$QTL <- k
    cand_l <- c(QTL = cand_l$QTL, cand_l[names(cand_l) != "QTL"])
    cand_k[[which(unique(key_cand_unique$QTL) == k)]] <- as.data.frame(t(cand_l))

  }
  df <- bind_rows(cand_k)
  df$QTL <- as.character(df$QTL)

  QTL_list_recalc <- left_join(x = QTL_list_recalc, y = df, by = "QTL")


  # write.csv(candidate_list,
  #           file = file.path(outdir, paste0("TP72_",i, "_CandList.csv")))
}

Candlist <- QTL_list_recalc[,-c(1:15)]

qtl_n <- Candlist[1,]
vec_n <- unlist(qtl_n, use.names = FALSE)
mat_n <- matrix(vec_n, nrow = 10)
mat_n <- unique(as.data.frame(t(mat_n)))
mat_n <- as.vector(t(mat_n))
cand_mat <- as.data.frame(t(as.data.frame(mat_n)))

for (n in 2:length(rownames(Candlist))) {
  qtl_n <- Candlist[n,]
  vec_n <- unlist(qtl_n, use.names = FALSE)
  if (!is.null(vec_n)) {
    mat_n <- matrix(vec_n, nrow = 10)
    mat_n <- unique(as.data.frame(t(mat_n)))
    mat_n <- as.vector(t(mat_n))
    mat_n <- as.data.frame(t(as.data.frame(mat_n)))
  } else {
    mat_n <- as.data.frame(matrix(rep(NA, 10), ncol = 10))
  }
  cand_mat <- pad_and_rbind(list(cand_mat, mat_n))
}

QTL_list_recalc <- cbind(QTL_list_recalc[,c(1:15)], cand_mat)
write.csv(QTL_list_recalc, "QTLCandidate_List.csv")



phe_in <- phe_72
ds_in <- create_gds$dosage
colnames(ds_in) <- paste0(create_gds$snp.chromosome, "_",
                          create_gds$snp.position)

setwd("/Users/okapi/Desktop/LazyGas物置/")
outdir <- "./CMplot"
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

###########################################
# Loop regression analysis over phenotypes
p_list <- NULL
for(i in seq_len(ncol(phe_in[, -(1:2)]))){
  phe_i <- phe_in[, i + 2]

  ###############################################
  # Loop regression over markers
  p_values <- apply(ds_in, 2, function(g){

    if(length(unique(na.omit(g))) == 1){
      return(rep(NA, 4))
    }

    #add = ifelse(g == 0, -1, ifelse(g == 60, 1, 0))
    #add = ifelse(g == 15, -1, ifelse(g == 45, 1, ifelse(g == 20, -2/3, ifelse(g == 40, 2/3, ifelse(g == 12, -6/5, ifelse(g == 24, -2/5, ifelse(g == 36, 2/5, ifelse(g == 48, 6/5, 0))))))))
    #add = ifelse(g == 0, -2, ifelse(g == 12, -6/5, ifelse(g == 15, -1, ifelse(g == 20, -2/3, ifelse(g == 24, -2/5, ifelse(g == 30, 0, ifelse(g == 36, 2/5, ifelse(g == 40, 2/3, ifelse(g == 45, 1, ifelse(g == 48, 6/5, 0))))))))))
    #dom = ifelse(g == 15, -2, ifelse(g == 45, 2, ifelse(g == 20, -4/3, ifelse(g == 40, 4/3, ifelse(g == 12, -4/5, ifelse(g == 24, -4/5, ifelse(g == 36, 4/5, ifelse(g == 48, 12/5, 0))))))))
    #dom = ifelse(g == 12, -6/5, ifelse(g == 15, -1, ifelse(g == 20, -1, ifelse(g == 24, -2/5, ifelse(g == 36, 2/5, ifelse(g == 40, 2/3, ifelse(g == 45, 1, ifelse(g == 48, 6/5, 0))))))))
    tmpdf <- data.frame(phe = phe_i,
                        add = ifelse(g == 0, -1, ifelse(g == 60, 1, 0)),
                        dom = as.numeric(g %in% c(12, 15, 20, 24, 30, 36, 40, 45, 48)),
                        dose = ifelse(g == 12, -6/5, ifelse(g == 15, -1, ifelse(g == 20, -2/3, ifelse(g == 24, -2/5, ifelse(g == 36, 2/5, ifelse(g == 40, 2/3, ifelse(g == 45, 1, ifelse(g == 48, 6/5, 0))))))))
    )
    tmpdf <- na.omit(tmpdf)  # 欠損のある行を除去
    if(nrow(tmpdf) == 0){
      return(rep(NA, 4))  # 全部欠損ならスキップ
    }
    res <- lm(phe ~ add + dom + dose, data = tmpdf)
    # 以下略

    s <- summary(res)
    f <- s$fstatistic
    if(is.numeric(f)){
      p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
    } else {
      return(rep(NA, 4))
    }
    coef <- s$coefficients
    coef <- data.frame(t(coef[, grepl("Pr", colnames(coef))]))
    out <- data.frame(all = NA, add = NA, dom = NA, dose = NA)
    if(!is.null(p)){
      out$all <- p
    }
    if(!is.null(coef$add)){
      out$add <- coef$add
    }
    if(!is.null(coef$dom)){
      out$dom <- coef$dom
    }
    if(!is.null(coef$dose)){
      out$dose <- coef$dose
    }
    return(out)
  })

  p_values <- do.call("rbind", p_values)
  colnames(p_values) <- c("P.all", "P.add", "P.dom", "P.dose")

  df <- p_values
  id <- rownames(df)
  chr <- sub("_.*", "", id)
  pos <- as.numeric(sub(".*_", "", id))
  q_val <- p.adjust(df$P.all, "fdr")

  p_values <- data.frame(Marker = id,
                         Chr = chr,
                         Position = pos)
  q_val <- as.data.frame(q_val)
  p_values <- cbind(p_values, df, q_val)
  p_values$"FDR<0.10" <- q_val < 0.1
  p_list <- c(p_list, list(p_values))
}

names(p_list) <- colnames(phe_in[, -(1:2)])
for(i in names(p_list)){

  write.csv(p_list[[i]], file = file.path(outdir,
                                          paste0(i, "_glmResult.csv")))
}
setwd("/Users/okapi/Desktop/LazyGas物置/")

for(i in seq_along(p_list)){
  df <- p_list[[i]]
  SNP <- rownames(df)
  CHR <- sub("_.*", "", SNP)
  BP <- as.numeric(sub(".*_", "", SNP))
  a <- as.data.frame(p_list[i])
  P <- a[,4]
  df <- data.frame(SNP = SNP, CHR = CHR, BP = BP, P = P)

  GWAS_result <- as.data.frame(p_list[i])
  GWAS_result <- as.data.frame(p_list[[i]])[,4:8]
  colnames(GWAS_result) <- c("p_val", "p.add", "p.dom", "p.dose", "q_val")

  qfilt <- dplyr::filter(GWAS_result, q_val < 0.1)
  if (nrow(qfilt) == 0) {
    fdr0.10 <- NA
  } else {
    top <- qfilt[which.max(qfilt$p_val), ]
    fdr0.10 <- top$p_val * 0.10 / top$q_val
  }

  fname <- paste0(names(p_list)[i], "_Manhattan")

  if (is.na(fdr0.10)) {
    CMplot(df,
           type = "h",
           plot.type = "m",
           file.name = fname,
           file = "pdf",
           dpi = 300,
           file.output = TRUE,
           verbose = TRUE,
           width = 15,
           height = 8)
  } else {
    CMplot(df,
           type = "h",
           plot.type = "m",
           threshold = fdr0.10,
           threshold.lty = 1,
           threshold.lwd = 1,
           threshold.col = "red",
           amplify = F,
           file.name = fname,
           file = "pdf",
           dpi = 300,
           file.output = TRUE,
           verbose = TRUE,
           width = 15,
           height = 8)
  }
}

