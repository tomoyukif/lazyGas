library(lazyGas)
library(GBScleanR)
library(ggplot2)
library(GenovisR)
library(SNPRelate)
library(SeqArray)
library(Biostrings)
library(dplyr)

#closeGDS(gds)

anuelist <- read.csv("R/dev/LazyGas_Anueploid_forFurutaSensei/anueploidlist.csv", row.names = 1)
anuelist <- anuelist[-c(60, 98),]
fdr_boarder <- "0.10"
peak_list_allpheno <- NULL

ID <- read.csv("R/dev/LazyGas_Anueploid_forFurutaSensei/Book2.csv")

phe_data <- read.csv("R/dev/LazyGas_Anueploid_forFurutaSensei/Phenotype_緊急事態.csv")
phe_72 <- phe_data[phe_data$Population_ID == "21TP72",]
phe_72$Serial_N <- formatC(phe_72$Serial_N,width=3,flag="0")
phe_72 <- phe_72[-c(60, 98),]
phe_72 <- phe_72[, -c(12, 13, 17)]
phe_vec <- c("Days to Heading", "Shattering", "Awn Length", "Culm Length", "Panicle Length", "Panicle Number",
             "Grain Number", "Spikelet Namber", "Seed Fertility","Grain Area Size", "Grain Length", "Grain Width")
phe_vec_2 <- c("Dh", "Sh", "Al", "Cl", "Pl", "Pn", "Gn", "Sn", "Sf", "Ga", "Gl", "Gw")



dosage_3somy <- read.csv("R/dev/LazyGas_Anueploid_forFurutaSensei/TP72genotype_ploidy3.csv", row.names = 1)
dosage_3somy <- dosage_3somy * 20
dosage_3somy <- dosage_3somy[-c(60, 98, 191, 192),]
colnames(dosage_3somy) <- 1:length(colnames(dosage_3somy))
dosage_3somy <- as.matrix(dosage_3somy)

dosage_5somy <- read.csv("R/dev/LazyGas_Anueploid_forFurutaSensei/TP72genotype_ploidy5.csv", row.names = 1)
dosage_5somy <- dosage_5somy * 12
dosage_5somy <- dosage_5somy[-c(60, 98, 191, 192),]
colnames(dosage_5somy) <- 1:length(colnames(dosage_5somy))
dosage_5somy <- as.matrix(dosage_5somy)


gds <- loadGDS("R/dev/LazyGas_Anueploid_forFurutaSensei/disallow_missing_genotype.gds",
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


#歪み除去のため組み替え近傍マーカー除こう

geno_ancestor <- read.csv("R/dev/LazyGas_Anueploid_forFurutaSensei/20TP41-6.csv", row.names = 1) # Your file
anc_marker <- rownames(geno_ancestor)
anc_ds <- geno_ancestor[,1]
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
deldis <- 300000
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


































temp_dir <- tempdir()
create_gds <-  list(genotype = dosage,
                    sample.id = sample.id,
                    snp.id = snp.id,
                    snp.rs.id = NULL,
                    snp.chromosome = snp.chromosome,
                    snp.position = snp.position,
                    snp.allele = snp.allele,haplotype=haplotype,
                    dosage = dosage)



sample_gds2 <- tempfile("sample", temp_dir, ".gds")
lg <- buildLazyGas(gds_fn = sample_gds2,
                   create_gds = create_gds)


#
# create_gds$sample.id
#
# lg[[sample.id]]

###############################################################################################################################
recalc_df <- NULL
peak_df <- NULL
peak_c_df <- NULL
for (j in 1:length(phe_vec)) {

  pheno <- phe_72[, j + 2]
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

  conv_fun <- function(g) {
    add <- ifelse(g == 0, -1, ifelse(g == 60, 1, 0))
    dom <- as.numeric(g %in% c(12, 15, 20, 24, 30, 36, 40, 45, 48, 60))
    dose <- ifelse(g == 15, -1, ifelse(g == 45, 1, ifelse(g == 20, -2/3, ifelse(g == 40, 2/3, ifelse(g == 12, -6/5, ifelse(g == 24, -2/5, ifelse(g == 36, 2/5, ifelse(g == 48, 6/5, 0))))))))
    out <- data.frame(add = add, dom = dom, dose = dose)
    return(out)
  }

  # For the model matrix above, the formula can be the following.
  formula <- "add + dom + dose"

  scanAssoc(object = lg,
            formula = formula,
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
  recalcAssoc(object = lg, n_threads = 10)

  for(i in seq_along(pheno$pheno_names)){
    p <- plotPeaks(object = lg, pheno = i, recalc = TRUE)
    p <- p + labs(title = pheno$pheno_names[i])
    print(p)
  }

  #TF注意
  for(i in seq_along(pheno$pheno_names)){
    pdf(file = paste0(phe_vec_2[j], "_", fdr_boarder,"_boxplot.pdf"), onefile = T)
    out <- haploPlot(object = lg, pheno = i, recalc =F)
    print(out)
    dev.off()
  }


  #############################################################################################################################
  peakcall <- lazyData(object = lg, dataset = "peakcall", pheno = pheno$pheno_names[1])
  write.csv(peakcall, paste0("allpeak_TP72_", phe_vec_2[j], "_", fdr_boarder, ".csv"))
  recalc <- lazyData(object = lg, dataset = "recalc", pheno = pheno$pheno_names[1])
  write.csv(recalc, paste0("TP72_", phe_vec_2[j], "_", fdr_boarder, ".csv"))
  makeInteractiveSummary(object = lg, pheno = pheno$pheno_names[1], out_fn = paste0("TP72_", phe_vec_2[j], "_", fdr_boarder, ".html"))


  #############################################################################################################################
  if (is.data.frame(recalc) == T) {
    for (peak in 1:recalc$peak_ID[length(rownames(recalc))]) {
      peak_i <- filter(recalc, peak_ID == peak)
      peak_SNP <- filter(peak_i, Dist2peak == 0)
      peak_SNP$QTL <- paste0("q72", phe_vec_2[j], sub("chr", "",peak_i$peak_Chr[1]), "_", sprintf("%05d",floor(peak_i$peak_Pos[1] / 1000)))
      peak_SNP$SNP <- paste0(peak_i$peak_Chr[1], "_", peak_i$peak_Pos[1])
      peak_SNP$start <- peak_i$Pos[1]
      peak_SNP$end <- peak_i$Pos[length(rownames(peak_i))]
      peak_SNP$pheno <- pheno$pheno_names[1]
      peak_df <- rbind(peak_df, peak_SNP)
    }
  }

  peakscan <- lazyData(object = lg, dataset = "scan", pheno = pheno$pheno_names[1])
  if (is.data.frame(peakcall) == T) {
    for (peak in 1:peakcall$peak_ID[length(rownames(peakcall))]) {
      peak_c_i <- filter(peakcall, peak_ID == peak)
      peak_c_SNP <- filter(peak_c_i, dist2peak == 0)
      peak_c_SNP$QTL <- paste0("q72", phe_vec_2[j], sub("chr", "",peak_c_i$peak_Chr[1]), "_", sprintf("%05d",floor(peak_c_i$peak_Pos[1] / 1000)))
      peak_c_SNP$SNP <- paste0(peak_c_i$peak_Chr[1], "_", peak_c_i$peak_Pos[1])
      peak_c_SNP$start <- peak_c_i$Pos[1]
      peak_c_SNP$end <- peak_c_i$Pos[length(rownames(peak_c_i))]
      peak_c_SNP$pheno <- pheno$pheno_names[1]
      scan_i <- filter(peakscan, Chr == peak_c_SNP$peak_Chr, Pos == peak_c_SNP$peak_Pos)[, 4:11]
      peak_c_SNP <- cbind(peak_c_SNP, scan_i)
      peak_c_df <- rbind(peak_c_df, peak_c_SNP)
    }
  }

}

QTL_list_recalc <- data.frame(Phenotype = peak_df$pheno,
                              Population = rep("21TP72", length(rownames(peak_df))),
                              QTL = peak_df$QTL,
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
                              PEV = peak_df$PVE,
                              q_val = peak_df$FDR)

QTL_NA_vec <- NULL
for (i in 1:length(rownames(QTL_list_recalc))) {
  QTL_i <- QTL_list_recalc$peak[i]
  QTL_NA <- filter(na_df, SNP == QTL_i)
  QTL_NA_vec <- c(QTL_NA_vec, QTL_NA$NA_ratio)
}
QTL_list_recalc$NA_ratio <- QTL_NA_vec
write.csv(QTL_list_recalc, paste0("QTLlit_TP72_recalc_", fdr_boarder,".csv"))


QTL_list_peakcall <- data.frame(Phenotype = peak_c_df$pheno,
                                Population = rep("21TP72", length(rownames(peak_c_df))),
                                QTL = peak_c_df$QTL,
                                peak = peak_c_df$SNP,
                                Start = peak_c_df$start,
                                End = peak_c_df$end,
                                P.model = peak_c_df$P.model,
                                P.add = peak_c_df$P.add,
                                P.dom = peak_c_df$P.dom,
                                P.dose = peak_c_df$P.dose,
                                Coef.add = peak_c_df$Coef.add,
                                Coef.dom = peak_c_df$Coef.dom,
                                Coef.dose = peak_c_df$Coef.dose,
                                PEV = peak_c_df$PVE,
                                q_val = peak_c_df$FDR)

QTL_NA_vec <- NULL
for (i in 1:length(rownames(QTL_list_peakcall))) {
  QTL_i <- QTL_list_peakcall$peak[i]
  QTL_NA <- filter(na_df, SNP == QTL_i)
  QTL_NA_vec <- c(QTL_NA_vec, QTL_NA$NA_ratio)
}
QTL_list_peakcall$NA_ratio <- QTL_NA_vec
write.csv(QTL_list_peakcall, paste0("QTLlit_TP72_allpeak_", fdr_boarder,".csv"))




