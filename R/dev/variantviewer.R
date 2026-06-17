scan_df <- read.csv("2025_list_variants/Rout/scan_out.csv")
sel_df <- read.csv("docs/2025/GWAS 2025 selected genes 20250430.csv")

library(rtracklayer)
gff_fn <- "~/01_wd/genomeData/rice/cultivar_sativa/nb_combined/version_2023/gff/nb_combined_all_annotated.gff"
gff <- import.gff(gff_fn)

library(ggplot2)
library(plotly)
library(htmltools)
library(reactable)

out_dir <- "2025_list_variants/Rout/htmlviewer"
dir.create(out_dir)
dataset <- unique(scan_df$dataset)
scan_df$GWAS.file <- gsub("/", "_", scan_df$GWAS.file)
scan_df$GWAS.file <- gsub(" ", "_", scan_df$GWAS.file)
for(i in seq_along(dataset)){
  scan_df_i <- scan_df[scan_df$dataset == dataset[i], ]
  scan_df_i$snp_id <- factor(scan_df_i$snp_id, sort(unique(scan_df_i$snp_id)))
  scan_df_i$snp_id <- as.numeric(scan_df_i$snp_id)
  base_info <- scan_df_i[1, 1:8]
  ann_i <- strsplit(scan_df_i$ann, "\\|")
  ann_i <- do.call("rbind", ann_i)
  ann_i <- data.frame(Gene_ID = scan_df_i$Gene.ID[1],
                      ID = scan_df_i$snp_id,
                      P_value = scan_df_i$p_value,
                      negLog10p = scan_df_i$neglog10p,
                      ann_i[, c(1, 3, 7, 13, 10, 14, 11)])
  names(ann_i)[-(1:4)] <- c("Effective_Allele", "Effect", "Transcript_ID",
                            "Position_in_CDS", "Change_in_DNA",
                            "Position_in_AA", "Change_in_AA")
  cds <- gff[gff$type == "CDS" & gff$gene_id == base_info$Gene.ID]
  cds_df <- data.frame(Transcript_ID = unlist(cds$Parent), cds_start = start(cds), cds_end = end(cds))
  base_info$gene_chr <- as.character(seqnames(gff[gff$ID == base_info$Gene.ID]))
  base_info$gene_start <- start(gff[gff$ID == base_info$Gene.ID])
  base_info$gene_end <- end(gff[gff$ID == base_info$Gene.ID])
  base_info$strand <- as.character(strand(gff[gff$ID == base_info$Gene.ID]))

  plot_df <- data.frame(Chr = scan_df_i$V4,
                        Start = scan_df_i$V5,
                        var_len = nchar(sub(",.+", "", scan_df_i$V6)),
                        Allele = scan_df_i$V6,
                        ann_i)
  plot_df$End <- plot_df$Start + plot_df$var_len - 1
  inter_gene <- grepl(base_info$Gene.ID, plot_df$Transcript_ID)
  if(sum(inter_gene) > 0){
    check <- duplicated(rbind(plot_df[!inter_gene, c("ID", "Effect")],
                              plot_df[inter_gene, c("ID", "Effect")]))
    check <- check[seq(sum(!inter_gene) + 1, length.out = sum(inter_gene))]
    inter_gene[inter_gene] <- inter_gene[inter_gene] & check
    plot_df <- plot_df[!inter_gene, ]
  }

  scale <- 5
  plot_df$Transcript_ID <- factor(plot_df$Transcript_ID, sort(unique(plot_df$Transcript_ID)))
  plot_df$ymin <- as.numeric(plot_df$Transcript_ID)
  plot_df$ymin <- max(plot_df$ymin) - plot_df$ymin + scale * 1.2
  plot_df$ymax <- plot_df$ymin + 0.8
  plot_df$y_pos <- plot_df$ymin + (plot_df$ymax - plot_df$ymin) / 2
  y_lab <- rev(sort(unique(plot_df$y_pos)))
  plot_df$x_pos <- plot_df$Start + (plot_df$End - plot_df$Start) / 2
  plot_df$Effect <- factor(plot_df$Effect, c("HIGH", "MODERATE", "LOW", "MODIFIER"))

  cds_df$Transcript_ID <- factor(cds_df$Transcript_ID, levels(plot_df$Transcript_ID))
  cds_df$ymin <- as.numeric(cds_df$Transcript_ID)
  cds_df$ymin <- max(cds_df$ymin) - cds_df$ymin + scale * 1.2
  cds_df$ymax <- cds_df$ymin + 0.8

  plot_man <- subset(plot_df, select = c(ID, Chr, Start, End,
                                         Allele, P_value, negLog10p,
                                         y_pos, x_pos,
                                         Position_in_CDS:Change_in_AA))
  plot_man$scaled_score <- plot_man$negLog10p / max(plot_man$negLog10p, na.rm = TRUE) * scale
  plot_man <- unique(plot_man)

  p <- ggplot() +
    geom_rect(data = cds_df,
              aes(ymin = ymin, ymax = ymax, xmin = cds_start, xmax = cds_end),
              fill = "gray99") +
    geom_point(data = plot_df,
               aes(y = y_pos, x = x_pos, color = Effect,
                   ID = ID,
                   Start = Start,
                   End = End,
                   Allele = Allele,
                   Position_in_CDS = Position_in_CDS,
                   Change_in_DNA = Change_in_DNA,
                   Position_in_AA = Position_in_AA,
                   Change_in_AA = Change_in_AA,
                   P_value = P_value,
                   negLog10p = negLog10p), size = 0.8) +
    geom_point(data = plot_man,
               aes(y = scaled_score, x = x_pos,
                   ID = ID,
                   Start = Start,
                   End = End,
                   Allele = Allele,
                   Position_in_CDS = Position_in_CDS,
                   Change_in_DNA = Change_in_DNA,
                   Position_in_AA = Position_in_AA,
                   Change_in_AA = Change_in_AA,
                   P_value = P_value,
                   negLog10p = negLog10p), size = 0.8) +
    scale_color_manual(breaks = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
                       values = c("magenta", "green", "blue", "gray30"),
                       labels = c("HIGH", "MODERATE", "LOW", "MODIFIER")) +
    scale_y_continuous(breaks = y_lab, labels = levels(plot_df$Transcript_ID)) +
    ylab("") +
    xlab("Physical position (bp)")

  plot_height <- length(levels(plot_df$Transcript_ID)) * 50 + 100
  if(plot_height < 300){
    plot_height <- 300
  }

  tag_list <- tagList(div(h1(base_info$GWAS.file),
                          style = "text-align:center"))
  tag_list <- tagList(tag_list,
                      div(h1(base_info$Gene.ID), style = "text-align:center"))
  tag_list <- tagList(tag_list,
                      div(h1(paste0(base_info$gene_chr, ": ",
                                    base_info$gene_start, "..",
                                    base_info$gene_end, ": ",
                                    base_info$strand)),
                          style = "text-align:center"))
  tag_list <- tagList(tag_list,
                      div(h1(base_info$Annotation), style = "text-align:center"))
  tag_list <- tagList(tag_list,
                      div(h1(base_info$X), style = "text-align:center"))
  tag_list <- tagList(tag_list,
                      div(h1(base_info$Expression.pattern), style = "text-align:center"))
  tag_list <- tagList(tag_list,
                      div(ggplotly(p, height = plot_height), style = "margin:auto;width:95vw"))

  geno_df <- subset(scan_df_i, select = -(c(GWAS.file:V3, ann:dataset)))
  names(geno_df)[1:4] <- c("Chr", "Pos", "Allele", "ID")
  geno_df <- unique(geno_df)
  geno_df$p_value <- signif(geno_df$p_value, digits = 3)
  geno_df$neglog10p <- signif(geno_df$neglog10p, digits = 3)
  names(geno_df)[grep("p_value", names(geno_df))] <- "P_value"
  names(geno_df)[grep("neglog10p", names(geno_df))] <- "negLog10p"
  table <- reactable(data = geno_df, sortable = TRUE,
                     resizable = TRUE, filterable = TRUE,
                     searchable = TRUE, showPageSizeOptions = TRUE,
                     wrap = TRUE, striped = TRUE, onClick = JS())
  tag_list <- tagList(tag_list,
                      div(h1("Variant list"), style = "text-align:center"),
                      div(table, style = "margin:auto;width:95vw;"))
  save_html(tag_list,
            file = file.path(out_dir,
                             paste(paste(base_info$GWAS.file, base_info$Gene.ID, sep = "_"),
                                   "html", sep = ".")))
  plot_df <- subset(plot_df, select = c(Chr, Start, End, Allele,
                                        Gene_ID,
                                        P_value:Change_in_AA))
  write.csv(plot_df, file.path(out_dir,
                               paste(paste(base_info$GWAS.file, base_info$Gene.ID, sep = "_"),
                                     "csv", sep = ".")), row.names = FALSE)
}


