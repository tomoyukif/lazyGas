library(ggplot2)
library(rtracklayer)
library(vcfR)
library(dplyr)
library(GBScleanR)
library(rrBLUP)

buildLazyQTL <- function(geno,
                         pheno,
                         sampleID_geno = NULL,
                         sampleID_pheno = NULL,
                         marker_chr,
                         marker_pos,
                         geno_fmt = "dosage",
                         geno_levels = 0:2){
  # Validation of genotype and phenotype data
  if(is.null(sampleID_pheno)){
    sampleID_pheno <- pheno[, grepl("ID|id", colnames(pheno))]
    pheno <- subset(pheno, select = !grepl("ID|id", colnames(pheno)))
  }
  if(is.null(sampleID_geno)){
    sampleID_geno <- rownames(geno)
  }
  nopheno <- sampleID_geno[!sampleID_geno %in% sampleID_pheno]
  nogeno <- sampleID_pheno[!sampleID_pheno %in% sampleID_geno]
  if(length(nopheno) > 0){
    message("The following samples have no phenotype info: \n",
            paste(nopheno, collapse = " "))
  }
  if(length(nogeno) > 0){
    message("The following samples have no genotype info: \n",
            paste(nogeno, collapse = " "))
  }

  # Validate genotype format
  uniq_geno <- na.omit(unique(as.vector(unlist(geno))))
  if(!all(uniq_geno %in% geno_levels)){
    stop("The input geno contains the following levels: \n",
         paste(sort(uniq_geno), collapse = " "),
         "\n But you specified geno_levels as: \n",
         paste(sort(geno_levels), collapse = " "),
         call. = FALSE)
  }
  if(any(is.na(suppressWarnings(as.numeric(uniq_geno))))){
    stop("The input geno contains the level(s) that cannot be coerced",
         "into numeric value(s). \n",
         'Set geno_fmt = "haplotype" if your geno indicates haplotype info.',
         call. = FALSE)
  }

  # Reorder phenotype data to match it with genotype data
  pheno <- subset(pheno, subset = !sampleID_pheno %in% nogeno)
  geno <- subset(geno, subset = !sampleID_geno %in% nopheno)
  colnames(geno) <- NULL
  sampleID_geno <- sampleID_geno[!sampleID_geno %in% nopheno]
  hitid <- match(sampleID_geno, sampleID_pheno)
  pheno <- pheno[hitid, ]

  # Calculate MAF
  if(geno_fmt == "dosage"){
    maf <- apply(geno, 2, function(x){
      tabx <- table(factor(x, geno_levels))
      tabx <- c(sum(tabx * rev(geno_levels)),
                sum(tabx * geno_levels))
      return(min(tabx)/sum(tabx))
    })
    maf <- unname(maf)

  } else {
    maf <- apply(geno, 2, function(x){
      tabx <- table(factor(x, geno_levels))
      return(min(tabx)/sum(tabx))
    })
    maf <- unname(maf)
  }

  out <- list(geno = geno,
              pheno = pheno,
              pheno_names = colnames(pheno),
              sample_id = sampleID_geno,
              marker_chr = marker_chr,
              marker_pos = marker_pos,
              maf = maf,
              geno_fmt = geno_fmt,
              geno_levels = geno_levels)
  class(out) <- c(class(out), "lazyQTL")
  gc();gc()
  return(out)
}


buildLazyGWAS <- function(gds,
                          pheno,
                          sampleID_pheno = NULL){
  # Validation of genotype and phenotype data
  if(is.null(sampleID_pheno)){
    sampleID_pheno <- pheno[, grepl("ID|id", colnames(pheno))]
    pheno <- subset(pheno, select = !grepl("ID|id", colnames(pheno)))
  }
  sampleID_geno <- getSamID(gds)
  nopheno <- sampleID_geno[!sampleID_geno %in% sampleID_pheno]
  nogeno <- sampleID_pheno[!sampleID_pheno %in% sampleID_geno]
  if(length(nopheno) > 0){
    message("The following samples have no phenotype info: \n",
            paste(nopheno, collapse = " "))
  }
  if(length(nogeno) > 0){
    message("The following samples have no genotype info: \n",
            paste(nogeno, collapse = " "))
  }

  # Prepare genotype data
  geno <- seqGetData(gds, "$dosage")
  geno <- geno[, validMar(gds)]
  geno <- -1 * (geno - 1)
  rownames(geno) <- getSamID(gds)

  # Reorder phenotype data to match it with genotype data
  pheno <- subset(pheno, subset = !sampleID_pheno %in% nogeno)
  geno <- subset(geno, subset = !sampleID_geno %in% nopheno)
  colnames(geno) <- NULL
  sampleID_geno <- sampleID_geno[!sampleID_geno %in% nopheno]
  hitid <- match(sampleID_geno, sampleID_pheno)
  pheno <- pheno[hitid, ]

  # Calculate a genomic relationship matrix
  message("Calculate the genomic relationship matrix.")
  grm <- A.mat(geno, return.imputed = TRUE, min.MAF = 0, max.missing = 1)
  geno <- data.frame(MarkerName = getMarID(gds),
                     Chr = getChromosome(gds),
                     Pos = getPosition(gds),
                     t(grm$imputed))
  rownames(grm$A) <- paste0("X", rownames(grm$A))
  colnames(grm$A) <- paste0("X", colnames(grm$A))

  # Calculate MAF
  id <- getSamID(gds)
  gds <- setSamFilter(gds, id = id[!id %in% sampleID_geno])
  seqSetFilter(gds, variant.sel = validMar(gds), sample.sel = validSam(gds))
  maf <- seqAlleleFreq(gds, minor = TRUE)

  out <- list(geno = geno,
              grm = grm$A,
              pheno = pheno,
              pheno_names = colnames(pheno),
              sample_id = sampleID_geno,
              maf = maf)
  class(out) <- c(class(out), "lazyGWAS")
  gc();gc()
  return(out)
}


print.lazyQTL <- function(x){
  message("Genotype")
  if(nrow(x$geno) >= 5){
    r_i <- seq_len(5)
  } else {
    r_i <- seq_len(nrow(x$geno))
  }
  if(ncol(x$geno) >= 5){
    c_i <- seq_len(5)
  } else {
    c_i <- seq_len(ncol(x$geno))
  }
  print(x$geno[r_i, c_i])
  message("Phenotype")
  print(head(x$pheno))
  message("Sample ID")
  print(head(x$sample_id))
  message("Marker position")
  print(rbind(head(x$marker_chr), head(x$marker_pos)))
  message("Minor allele frequency")
  print(head(x$maf))
  message("Genotype format")
  print(head(x$geno_fmt))
  message("Genotype levels")
  print(head(x$geno_levels))
}

print.lazyGWAS <- function(x){
  message("Genotype")
  if(nrow(x$geno) >= 7){
    r_i <- seq(3, 7)
  } else {
    r_i <- seq_len(nrow(x$geno))
  }
  if(ncol(x$geno) >= 7){
    c_i <- seq(3, 7)
  } else {
    c_i <- seq_len(ncol(x$geno))
  }
  print(x$geno[r_i, c_i])
  message("Phenotype")
  print(head(x$pheno))
  message("Sample ID")
  print(head(x$sample_id))
  message("Marker position")
  print(rbind(head(x$geno[, 1]), head(x$geno[, 2])))
  message("Minor allele frequency")
  print(head(x$maf))
}

scanQTL <- function(x, ...){
  UseMethod("scanQTL", x)
}

scanQTL.lazyQTL <- function(x,
                            out_fn,
                            formula = "",
                            makeDF_FUN = NULL){
  if(!inherits(x, "lazyQTL")){
    stop("The input should be the lazyQTL object",
         call. = FALSE)
  }

  message("Going to analyse the following phenotypes: \n",
          paste(x$pheno_names, collapse = ", "))

  for(i in seq_along(x$pheno_names)){
    p_values <- apply(x$geno, 2, function(g){
      if(length(unique(na.omit(g))) == 1){
        return(rep(NA, 6))
      }
      if(x$geno_fmt == "dosage"){
        df <- makeDF(g = g,
                     phe = x$pheno[, i],
                     makeDF_FUN = makeDF_FUN,
                     formula = formula)

        if(length(na.omit(unique(df$phe))) == 2){
          out <- doGLM(df$df, df$fml)

        } else {
          out <- doLM(df$df, df$fml)
        }

      } else if(x$geno_fmt == "haplotype"){
        if(is.null(makeDF_FUN)){
          df <- data.frame(phe = x$pheno[, i], group = g)
          fml <- "phe ~ group"

        } else {
          df <-  data.frame(phe = x$pheno[, i], makeDF_FUN(g))
          if(ncol(df) > 2){
            stop("Haplotype association test accepts only one ",
                 "explanatory variable.\n",
                 "Please revise the makeDF_FUN.",
                 call. = FALSE)
          }
        }
        out <- doKruskal(df, fmt)
      }

      return(out)
    })

    if(x$geno_fmt == "haplotype"){
      p_values <- data.frame(P.fit = p_values)

    } else {
      if(is.list(p_values)){
        p_values <- data.frame(do.call("rbind", p_values))
      } else {
        p_values <- data.frame(t(p_values))
      }
    }

    p_values <- data.frame(Chr = x$marker_chr,
                           Pos = x$marker_pos,
                           MAF = x$maf,
                           p_values)
    colnames(p_values)[ncol(p_values)] <- "%var"
    p_values$FDR <- p.adjust(p_values$P.fit, "fdr")
    p_values$P.fit <- -log10(p_values$P.fit)
    write.csv(p_values,
              paste0(out_fn, x$pheno_names[i], "_scanQTL.csv"),
              row.names = FALSE)
    gc();gc()
  }
  out <- list(geno = x$geno,
              pheno_names = x$pheno_names,
              pvalues_fn = paste0(out_fn, x$pheno_names, "_scanQTL.csv"))
  class(out) <- c(class(out), "QTLscan")
  invisible(out)
}

makeDF <- function(g, phe, makeDF_FUN, formula){
  g <- as.numeric(g)

  if(is.null(makeDF_FUN)){
    df <- data.frame(add = g)
    fml <- formula("phe ~ add")

  } else {
    df <- makeDF_FUN(g)
    if(formula == ""){
      stop("Provide formula if you specified makeDF_FUN",
           call. = FALSE)
    }
    fml <- formula(paste0("phe ~ ", formula))
  }

  return(list(df = data.frame(phe = phe, df), fml = fml))
}

doGLM <- function(df, fml){
  res <- try(glm(formula = fml, data = df, family = binomial))
  if(inherits(res, "try-error")){ return(NA) }

  p <- pchisq(res$null.deviance - res$deviance,
              res$df.null - res$df.residual,
              lower.tail = FALSE)
  pervar <- 1 - res$deviance / res$null.deviance
  s <- summary(res)
  coef <- s$coefficients

  if(nrow(coef) == 1){ return(NA) }

  att <- attributes(res$terms)
  return(makeOut(p, coef, pervar, att$term.labels))
}

doLM <- function(df, fml){
  res <- try(lm(formula = fml, data = df))
  if(inherits(res, "try-error")){ return(NA) }
  s <- summary(res)
  f <- s$fstatistic
  pervar <- s$r.squared
  coef <- s$coefficients
  if(is.numeric(f)){
    p <- pf(f[1], f[2], f[3], lower.tail = FALSE)

  } else {
    return(NA)
  }
  att <- attributes(res$terms)
  return(makeOut(p, coef, pervar, att$term.labels))
}

doKruskal <- function(df, fml){
  res <- try(kruskal.test(formula = formula(fml), data = df))
  if(inherits(res, "try-error")){ return(NA) }
  return(res$p.value)
}

makeOut <- function(p, coef, pervar, terms){
  labs <- c("P.fit", paste("P", terms, sep= "."),
            paste("Coef", terms, sep = "."),
            "%var")
  out <- rep(NA, length(labs))
  out[1] <- p
  coef_row <- rownames(coef)
  for(j in seq_along(terms)){
    if(terms[j] %in% coef_row){
      out[1 + j] <- coef[terms[j], 4]
      out[1 + j + length(terms)] <- coef[terms[j], 1]
    }
  }
  out[2 + length(terms)*2] <- pervar
  names(out) <- labs
  return(out)
}

print.QTLscan <- function(x){
  message("Genotype")
  if(nrow(x$geno) >= 5){
    r_i <- seq_len(5)
  } else {
    r_i <- seq_len(nrow(x$geno))
  }
  if(ncol(x$geno) >= 5){
    c_i <- seq_len(5)
  } else {
    c_i <- seq_len(ncol(x$geno))
  }
  print(x$geno[r_i, c_i])
  message("File names storing regression results")
  print(cbind(x$pheno_names, x$pvalues_fn))
}

scanGWAS <- function(x, ...){
  UseMethod("scanGWAS", x)
}

scanGWAS.lazyGWAS <- function(x,
                              out_fn,
                              nPC = 2,
                              n_core = 1){
  for(i in seq_along(x$pheno_names)){
    pheno_df <- data.frame(id = paste0("X", x$sample_id),
                           pheno = x$pheno[, i])
    colnames(pheno_df)[2] <- x$pheno_names[i]
    gwas <- GWAS(pheno = pheno_df, geno = x$geno,
                 n.PC = nPC, K = x$grm, plot = FALSE, n.core = n_core)

    names(gwas)[4] <- "P.fit"
    gwas$P.fit[gwas$P.fit == 0] <- NA
    gwas$FDR <- p.adjust(10^-gwas$P.fit, method = "fdr")
    gwas$MAF <- x$maf
    gwas <- subset(gwas, select = c(Chr, Pos, MAF, P.fit, FDR))
    write.csv(gwas,
              paste0(out_fn, x$pheno_names[i], "_scanGWAS.csv"),
              row.names = FALSE)
    gc();gc()
  }

  out <- list(geno = x$geno,
              pheno_names = x$pheno_names,
              pvalues_fn = paste0(out_fn, x$pheno_names, "_scanGWAS.csv"))
  class(out) <- c(class(out), "GWASscan")
  invisible(out)
}

print.GWASscan <- function(x){
  message("Genotype")
  if(nrow(x$geno) >= 7){
    r_i <- seq(3, 7)
  } else {
    r_i <- seq_len(nrow(x$geno))
  }
  if(ncol(x$geno) >= 7){
    c_i <- seq(3, 7)
  } else {
    c_i <- seq_len(ncol(x$geno))
  }
  print(x$geno[r_i, c_i])
  message("File names storing regression results")
  print(cbind(x$pheno_names, x$pvalues_fn))
}

plotManhattan <- function(x, ...){
  UseMethod("plotManhattan", x)
}

plotManhattan.QTLscan <- function(x,
                                  pheno = NULL,
                                  chr = NULL,
                                  start = NULL,
                                  end = NULL,
                                  signif = NULL,
                                  out_fn = "",
                                  out_fmt = "pdf"){
  if(!is.null(pheno)){
    if(is.numeric(pheno)){
      phe_index <- pheno
    } else if(is.logical(pheno)){
      phe_index <- which(pheno)
    } else if(is.character(pheno)){
      phe_index <- which(x$pheno_names %in% pheno)
    }
  } else {
    phe_index <- seq_along(x$pvalues_fn)
  }
  for(i in phe_index){
    pvalues <- read.csv(x$pvalues_fn[i])
    tmp_fn <- paste0(out_fn, x$pheno_names[i], "_plotManhattan.", out_fmt)
    p <- plotManhattan.data.frame(pvalues, chr, start, end, signif, tmp_fn, out_fmt)
    print(p)
  }
}

plotManhattan.GWASscan <- function(x,
                                   pheno = NULL,
                                   chr = NULL,
                                   start = NULL,
                                   end = NULL,
                                   signif = NULL,
                                   out_fn = "",
                                   out_fmt = "pdf"){
  if(!is.null(pheno)){
    if(is.numeric(pheno)){
      phe_index <- pheno
    } else if(is.logical(pheno)){
      phe_index <- which(pheno)
    } else if(is.character(pheno)){
      phe_index <- which(x$pheno_names %in% pheno)
    }
  } else {
    phe_index <- seq_along(x$pvalues_fn)
  }
  for(i in phe_index){
    pvalues <- read.csv(x$pvalues_fn[i])
    tmp_fn <- paste0(out_fn, x$pheno_names[i], "_plotManhattan.", out_fmt)
    p <- plotManhattan.data.frame(pvalues, chr, start, end, signif, tmp_fn, out_fmt)
    invisible(p)
  }
}

plotManhattan.character <- function(x,
                                    chr = NULL,
                                    start = NULL,
                                    end = NULL,
                                    signif = NULL,
                                    out_fn = "",
                                    out_fmt = "pdf"){
  pvalues <- read.csv(x)
  p <- plotManhattan.data.frame(pvalues, chr, start, end, signif, out_fn, out_fmt)
  invisible(p)
}

plotManhattan.data.frame <- function(x,
                                     chr = NULL,
                                     start = NULL,
                                     end = NULL,
                                     signif = NULL,
                                     out_fn = "",
                                     out_fmt = "pdf"){
  if(is.null(signif)){
    signif <- rep(TRUE, nrow(x))
  }
  if(!is.null(chr)){
    signif <- signif[x$Chr == chr]
    x <- subset(x, subset = Chr == chr)
  }
  if(!is.null(start)){
    signif <- signif[x$Pos >= start]
    x <- subset(x, subset = Pos >= start)
  }
  if(!is.null(end)){
    signif <- signif[x$Pos <= end]
    x <- subset(x, subset = Pos <= end)
  }

  if(!is.logical(signif)){
    if(is.character(signif)){
      signif <- eval(parse(text = signif))
    } else {
      stop("signif should be a string or a vector of logical values.",
           call. = FALSE)
    }
  }

  signif[is.na(signif)] <- FALSE
  x_signif <- subset(x, subset = signif)
  x <- subset(x, subset = !signif)
  p <- ggplot() +
    geom_point(data = x,
               mapping = aes(x = Pos, y = P.fit),
               color = "darkgray",
               size = 1,
               shape = 20)
  if(nrow(x_signif) != 0){
    p <- p + geom_point(data = x_signif,
                        mapping = aes(x = Pos, y = P.fit),
                        color = "magenta", size = 1, shape = 18)
  }
  p <- p + facet_wrap(~ Chr,
                      nrow = 1,
                      scales = "free_x",
                      strip.position = "bottom") +
    ylab("-log10(P)") +
    xlab("Chromosome") +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          strip.text.x = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          axis.line.x.bottom = element_line(colour = "black"),
          panel.spacing.x = unit(0.2, "lines"),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "gray90"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.placement = "outside",
          strip.background = element_rect(fill = "white", colour = "white"))
  if(out_fmt == "pdf"){
    pdf(out_fn)
    print(p)
    dev.off()
  } else if(out_fmt == "png"){
    png(out_fn, res = 300, width = 480*4, height = 480*4)
    print(p)
    dev.off()
  }
  gc();gc()
  invisible(p)
}

callPeakBlock <- function(x, ...){
  UseMethod("callPeakBlock", x)
}

callPeakBlock.data.frame <- function(x,
                                     pvalues,
                                     signif,
                                     out_fn,
                                     rsquare = 0.6){
  geno <- matrix(as.numeric(as.matrix(x)), nrow(x), ncol(x))
  callPeakBlock.matrix(geno, pvalues, signif, out_fn, rsquare)
}

callPeakBlock.QTLscan <- function(x,
                                  signif,
                                  out_fn,
                                  rsquare = 0.6){
  geno <- matrix(as.numeric(as.matrix(x$geno)), nrow(geno), ncol(geno))
  for(i in seq_along(x$pvalues_fn)){
    pvalues <- read.csv(x$pvalues_fn[i])
    tmp_fn <- paste0(out_fn, x$pheno_names[i], "_peakBlock.csv")
    callPeakBlock.matrix(geno, pvalues, signif, tmp_fn, rsquare)
  }
  out <- list(peakblock_fn = paste0(out_fn, x$pheno_names, "_peakBlock.csv"),
              pheno_names = x$pheno_names)
  class(out) <- c(class(out), "peakCall")
  attributes(out) <- c(attributes(out), scan = "QTL")
  invisible(out)
}

callPeakBlock.GWASscan <- function(x,
                                   signif,
                                   out_fn,
                                   rsquare = 0.6){
  geno <- t(as.matrix(subset(x$geno, select = -c(MarkerName:Pos))))
  geno <- matrix(as.numeric(geno), nrow(geno), ncol(geno))
  for(i in seq_along(x$pvalues_fn)){
    pvalues <- read.csv(x$pvalues_fn[i])
    tmp_fn <- paste0(out_fn, x$pheno_names[i], "_peakBlock.csv")
    callPeakBlock.matrix(geno, pvalues, signif, tmp_fn, rsquare)
  }
  out <- list(peakblock_fn = paste0(out_fn, x$pheno_names, "_peakBlock.csv"),
              pheno_names = x$pheno_names)
  class(out) <- c(class(out), "peakCall")
  attributes(out) <- c(attributes(out), scan = "GWAS")
  invisible(out)
}

callPeakBlock.matrix <- function(x,
                                 pvalues,
                                 signif,
                                 out_fn,
                                 rsquare = 0.6){
  if(!is.numeric(x)){
    stop("The input x should be the numeric matrix.",
         call. = FALSE)
  }
  if(nrow(pvalues) != ncol(x)){
    stop("The number of markers does not match between geno and pvalues.",
         call. = FALSE)
  }

  if(!is.logical(signif)){
    if(is.character(signif)){
      signif <- eval(parse(text = sub("x\\$", "pvalues$", signif)))
    } else {
      stop("signif should be a string or a vector of logical values.",
           call. = FALSE)
    }
  }

  signif[is.na(signif)] <- FALSE
  pvalues <- subset(pvalues, subset = signif)
  rownames(pvalues) <- seq_len(nrow(pvalues))
  x <- subset(x, select = signif)
  x <- x + 1

  snpgds_fn <- tempfile(pattern = "lazyPeakcall", fileext = ".snpgds")
  snpgdsCreateGeno(snpgds_fn,
                   x,
                   seq_len(nrow(x)),
                   rownames(pvalues),
                   rownames(pvalues),
                   pvalues$Chr,
                   pvalues$Pos,
                   snpfirstdim = FALSE)
  snpgds <- snpgdsOpen(snpgds_fn)

  colnames(pvalues) <- sub("X.var", "%var", colnames(pvalues))

  write.table(t(c("Peak_id", "Peak_chr", "Peak_pos", "Peak_pval", "Dist2peak",
                  colnames(pvalues))),
              out_fn,
              sep = ",", quote = TRUE, row.names = FALSE, col.names = FALSE)

  peak_id <- 0
  while(TRUE){
    if(nrow(pvalues) == 0){
      break
    }
    peak_id <- peak_id + 1
    peak_index <- which.max(pvalues$P.fit)
    chr <- pvalues$Chr[peak_index]
    pos <- pvalues$Pos[peak_index]

    target_variants <- rownames(pvalues)[pvalues$Chr == chr]
    if(length(target_variants) == 1){
      peak_block <- target_variants

    } else {
      ld <- snpgdsLDMat(snpgds,
                        snp.id = target_variants,
                        slide = 0, num.thread = 1, method = "composite")
      peak_pos <- which(ld$snp.id == rownames(pvalues)[peak_index])
      peak_ld  <- ld$LD[, peak_pos]^2
      ld_block <- peak_ld >= rsquare
      ld_block[is.na(ld_block)] <- FALSE
      peak_block <- target_variants[ld_block]
    }
    peakblock_pvalues <- subset(pvalues,
                                subset = rownames(pvalues) %in% peak_block)
    write.table(data.frame(peak_id = peak_id,
                           peak_chr = chr,
                           peak_pos = pos,
                           peak_pval = max(peakblock_pvalues$P.fit),
                           dist2peak = peakblock_pvalues$Pos - pos,
                           peakblock_pvalues),
                out_fn,
                sep = ",", quote = TRUE, row.names = FALSE, col.names = FALSE,
                append = TRUE)
    x <- subset(x, select = !rownames(pvalues) %in% peak_block)
    pvalues <- subset(pvalues, subset = !rownames(pvalues) %in% peak_block)
  }
  snpgdsClose(snpgds)
  unlink(snpgds_fn, force = TRUE)
  gc();gc()
}

print.peakCall <- function(x){
  message("Peakcall file names")
  print(cbind(x$pheno_names, x$peakblock_fn))
  message("Scan type")
  print(attributes(x)$scan)
}

listCandidate <- function(x, ...){
  UseMethod("listCandidate", x)
}

listCandidate.peakCall <- function(x,
                                   annotation_fn,
                                   gff_fn,
                                   snpeff_fn = NULL,
                                   out_fn){
  if(attributes(x)$scan == "QTL"){
    scan <- "QTL"
  } else if(attributes(x)$scan == "GWAS"){
    scan <- "GWAS"
  }

  if(!is.null(snpeff_fn)){
    snpeff_fn <- read.vcfR(snpeff_fn, checkFile = FALSE, verbose = FALSE)
  }

  if(grepl("\\.csv$", basename(annotation_fn))){
    annotation_fn <- read.csv(annotation_fn)
  } else if(grepl("\\.tsv$", basename(annotation_fn))){
    annotation_fn <- read.table(annotation_fn, sep = "\t", header = TRUE)
  }

  for(i in seq_along(x$peakblock_fn)){
    peakblock <- read.csv(x$peakblock_fn[i])
    tmp_fn <- paste0(out_fn, x$pheno_names[i], "_candidateList.csv")
    listCandidate.data.frame(peakblock,
                             annotation_fn,
                             gff_fn,
                             snpeff_fn,
                             tmp_fn,
                             scan = scan)
  }
  out <- paste0(out_fn, x$pheno_names, "_candidateList.csv")
  names(out) <- x$pheno_names
  invisible(out)
}

listCandidate.data.frame <- function(x,
                                     annotation_fn,
                                     gff_fn,
                                     snpeff_fn = NULL,
                                     out_fn,
                                     scan = "QTL"){
  out <- NULL
  if(nrow(x) != 0){
    gff <- import.gff(gff_fn)
    gff <- gff[gff$type %in% c("transcript", "mRNA")]

    if(!is.data.frame(annotation_fn)){
      if(grepl("\\.csv$", basename(annotation_fn))){
        ann <- read.csv(annotation_fn)
      } else if(grepl("\\.tsv$", basename(annotation_fn))){
        ann <- read.table(annotation_fn, sep = "\t", header = TRUE)
      }
    }

    if(!is.null(snpeff_fn)){
      if(inherits(snpeff_fn, "vcfR")){
        snpeff <- snpeff_fn

      } else {
        snpeff <- read.vcfR(snpeff_fn, checkFile = FALSE, verbose = FALSE)
      }

    } else {
      snpeff <- NULL
    }

    for(i_peak in unique(x$Peak_id)){
      i_peakblock <- x[x$Peak_id == i_peak, ]

      if(scan == "QTL"){
        tmp <- data.frame(Peak_id = i_peak,
                          getQTLcandidate(gff, i_peakblock, snpeff))

      } else if(scan == "GWAS"){
        tmp <- data.frame(Peak_id = i_peak,
                          getGWAScandidate(gff, i_peakblock, snpeff))
      }
      out <- rbind(out, tmp)
    }
    if(!is.null(out)){
      out <- left_join(out, ann, by = c("GeneID", "TxID"))
    }
  }
  write.csv(out, out_fn, row.names = FALSE)
}

getQTLcandidate <- function(gff, peakblock, snpeff){
  peak_gff <- GRanges(seqnames = peakblock$Peak_chr[1],
                      ranges = IRanges(start = peakblock$Pos[1],
                                       end = tail(peakblock$Pos, 1)))
  hit <- gff[queryHits(findOverlaps(gff, peak_gff))]
  hit <- hit[order(start(hit))]
  hit$fromPeak <- start(hit) - peakblock$Peak_pos[1]
  nearest_p <- sapply(start(hit), function(x){
    return(peakblock$P.fit[which.min(abs(peakblock$Pos - x))])
  })
  out <- data.frame(P.fit = nearest_p,
                    Dist2peak = hit$fromPeak,
                    Gene_chr = as.character(seqnames(hit)),
                    Gene_start = start(hit),
                    GeneID = unlist(hit$Parent),
                    TxID = hit$ID)

  if(!is.null(snpeff)){
    snpeff_out <- addSnpEff(snpeff, peakblock, scan = "QTL")
    out <- left_join(out, snpeff_out, by ="TxID")
  }
  return(out)
}

getGWAScandidate <- function(gff, peakblock, snpeff){
  snpeff_out <- addSnpEff(snpeff, peakblock, scan = "GWAS")
  hit <- gff[match(snpeff_out$TxID, gff$ID)]
  hit$fromPeak <- start(hit) - peakblock$Peak_pos[1]
  out <- left_join(data.frame(Dist2peak = hit$fromPeak,
                              Gene_chr = as.character(seqnames(hit)),
                              Gene_start = start(hit),
                              GeneID = unlist(hit$Parent),
                              TxID = hit$ID),
                   snpeff_out, by = "TxID")
  return(out)
}

snpeff2df <- function(info, chr, pos, var_chr, var_pos){
  hit <- info[which(chr == var_chr & pos == var_pos)]
  if(length(hit) == 0){ return(NA) }
  hit <- sub(";LOF.*", "", hit)
  ann <- suppressWarnings(strsplit(sub(".*ANN=", "", hit), ",")[[1]])
  ann[grepl("\\|$", ann)] <- paste(ann[grepl("\\|$", ann)], " ")
  ann <- do.call("rbind", strsplit(ann, "\\|"))
  ann <- subset(ann, subset = ann[, 2] != "intergenic_region",
                select = c(2:4, 7))
  return(ann)
}

addSnpEff <- function(snpeff, peakblock, scan){
  fix <- getFIX(snpeff)
  chr <- fix[, 1]
  pos <- as.numeric(fix[, 2])
  info <- getINFO(snpeff)
  peak_ann <- NULL

  if(scan == "QTL"){
    target_chr <- peakblock$Chr[1] == chr
    target_pos <- pos >= min(peakblock$Pos) & pos <= max(peakblock$Pos)
    target_pos <- pos[target_chr & target_pos]
    for(i in seq_along(target_pos)){
      peak_ann <- rbind(peak_ann,
                        data.frame(snpeff2df(info, chr, pos,
                                             var_chr = peakblock$Chr[1],
                                             var_pos = target_pos[i])))
    }
    colnames(peak_ann) <- c("type", "impact", "GeneId", "TxID")
    genewise <- tapply(peak_ann$impact, peak_ann$TxId, function(x){
      as.vector(table(factor(x, c("HIGH", "MODERATE", "LOW", "MODIFIER"))))
    })
    genewise <- data.frame(TxID = names(genewise),
                           do.call("rbind", genewise))
    colnames(genewise)[-1] <- c("HIGH", "MODERATE", "LOW", "MODIFIER")

  } else if(scan == "GWAS"){
    for(i in seq_len(nrow(peakblock))){
      peak_ann <- rbind(peak_ann,
                        data.frame(peakblock$P.fit[i], peakblock$MAF[i],
                                   snpeff2df(info, chr, pos,
                                             var_chr = peakblock$Chr[i],
                                             var_pos = peakblock$Pos[i])))
    }
    colnames(peak_ann) <- c("P.fit", "MAF", "Type",
                            "Impact", "GeneId", "TxID")
    impacts <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
    genewise <- tapply(seq_along(peak_ann$Impact), peak_ann$TxID, function(i){
      data.frame(Max_Pval = max(peak_ann$P.fit[i]),
                 Min_Pval = min(peak_ann$P.fit[i]),
                 Maf_Max = peak_ann$MAF[i][which.max(peak_ann$P.fit[i])],
                 Maf_Min = peak_ann$MAF[i][which.min(peak_ann$P.fit[i])],
                 t(as.vector(table(factor(peak_ann$Impact[i], impacts)))))
    })
    genewise <- data.frame(TxID = names(genewise),
                           do.call("rbind", genewise))
    colnames(genewise)[-(1:5)] <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  }
  return(genewise)
}
