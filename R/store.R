################################################################################
# Companion storage for lazyGas results (Parquet default, SQLite optional, GDS legacy)
#' @importFrom jsonlite fromJSON write_json
#' @importFrom arrow write_parquet read_parquet
NULL

#' @importFrom gdsfmt exist.gdsn
.store_gds_fn <- function(object) {
  if (!is.null(object$filename) && nzchar(object$filename)) {
    return(normalizePath(object$filename, winslash = "/", mustWork = FALSE))
  }
  if (!is.null(object@gds.fn) && nzchar(object@gds.fn)) {
    return(normalizePath(object@gds.fn, winslash = "/", mustWork = FALSE))
  }
  stop("Cannot determine the GDS file path from the LazyGas object.", call. = FALSE)
}

.store_companion_path <- function(gds_fn, store_mode) {
  gds_fn <- normalizePath(gds_fn, winslash = "/", mustWork = FALSE)
  stem <- sub("\\.gds$", "", gds_fn, ignore.case = TRUE)
  if (identical(stem, gds_fn)) {
    stem <- gds_fn
  }
  if (store_mode == "sqlite") {
    return(paste0(stem, ".lazygas.sqlite"))
  }
  paste0(stem, ".lazygas")
}

.store_init <- function(object, gds_fn, store = "parquet", overwrite = FALSE) {
  store <- match.arg(store, c("parquet", "sqlite", "gds", "auto"))
  gds_fn <- normalizePath(gds_fn, winslash = "/", mustWork = FALSE)

  if (store == "auto") {
    if (exist.gdsn(node = object$root, path = "lazygas")) {
      store <- "gds"
    } else {
      store <- "parquet"
    }
  }

  if (store == "gds") {
    object@store <- list(
      mode = "gds",
      gds_fn = gds_fn,
      path = NULL,
      peak_buffer = new.env(parent = emptyenv())
    )
    return(object)
  }

  path <- .store_companion_path(gds_fn = gds_fn, store_mode = store)
  if (overwrite && store == "parquet" && dir.exists(path)) {
    unlink(path, recursive = TRUE)
  }
  if (overwrite && store == "sqlite" && file.exists(path)) {
    unlink(path)
  }

  if (store == "parquet") {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "scan"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "peakcall"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "recalc"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "candidate"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "snpeff"), recursive = TRUE, showWarnings = FALSE)
  }

  object@store <- list(
    mode = store,
    gds_fn = gds_fn,
    path = path,
    peak_buffer = new.env(parent = emptyenv())
  )

  if (store != "gds" && exist.gdsn(node = object$root, path = "lazygas")) {
    companion_empty <- !dir.exists(path) ||
      length(list.files(path, all.files = TRUE, no.. = TRUE)) == 0L
    if (store == "sqlite") {
      companion_empty <- !file.exists(path)
    }
    if (companion_empty) {
      .store_migrate_from_gds(object)
    } else {
      message(
        "Companion store already exists; keeping it. ",
        "Use importLazyGasResults() to copy legacy GDS lazygas/ data."
      )
    }
  }

  object
}

.store_mode <- function(object) {
  if (length(object@store) == 0L || is.null(object@store$mode)) {
    return("gds")
  }
  object@store$mode
}

.store_is_gds <- function(object) {
  identical(.store_mode(object), "gds")
}

.store_path <- function(object) {
  object@store$path
}

.store_meta_path <- function(object) {
  file.path(.store_path(object), "meta.json")
}

.store_scan_path <- function(object, pheno_name) {
  root <- .store_path(object)
  if (.store_mode(object) == "sqlite") {
    return(paste0("scan_", .store_safe_name(pheno_name)))
  }
  file.path(root, "scan", paste0(.store_safe_name(pheno_name), ".parquet"))
}

.store_table_path <- function(object, section, kind, pheno_name) {
  root <- .store_path(object)
  if (.store_mode(object) == "sqlite") {
    return(paste0(section, "_", kind, "_", .store_safe_name(pheno_name)))
  }
  file.path(root, section, paste0(kind, "_", .store_safe_name(pheno_name), ".parquet"))
}

.store_safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

.store_read_meta <- function(object) {
  if (.store_is_gds(object)) {
    return(.store_read_meta_gds(object))
  }
  meta_path <- .store_meta_path(object)
  if (!file.exists(meta_path)) {
    return(list())
  }
  jsonlite::fromJSON(txt = meta_path, simplifyVector = TRUE)
}

.store_write_meta <- function(object, meta) {
  if (.store_is_gds(object)) {
    return(invisible(NULL))
  }
  meta_path <- .store_meta_path(object)
  jsonlite::write_json(meta, path = meta_path, auto_unbox = TRUE, pretty = TRUE)
  invisible(NULL)
}

.store_update_meta <- function(object, ...) {
  meta <- .store_read_meta(object)
  upd <- list(...)
  for (nm in names(upd)) {
    meta[[nm]] <- upd[[nm]]
  }
  .store_write_meta(object, meta)
}

.store_read_meta_gds <- function(object) {
  if (!exist.gdsn(node = object$root, path = "lazygas/scan")) {
    return(list())
  }
  nodes <- c("kruskal", "formula", "null_formula", "conv_fun", "geno_format")
  out <- list()
  for (n in nodes) {
    path <- paste0("lazygas/scan/", n)
    if (exist.gdsn(node = object$root, path = path)) {
      out[[n]] <- .get_data_gds_only(object, path)
    }
  }
  if (exist.gdsn(node = object$root, path = "lazygas/scan/fixed_effect")) {
    out$fixed_effect <- .get_data_gds_only(object, "lazygas/scan/fixed_effect")
    out$fixed_effect_col_names <- get.attr.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = "lazygas/scan/fixed_effect")
    )$col_names
  }
  if (exist.gdsn(node = object$root, path = "lazygas/peakcall/peaks")) {
    att <- gdsfmt::get.attr.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = "lazygas/peakcall/peaks")
    )
    out$peakcall_signif <- att$signif
    out$peakcall_threshold <- att$threshold
    out$peakcall_peaks_col_names <- att$col_names
    out$peakcall_blocks_col_names <- gdsfmt::get.attr.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = "lazygas/peakcall/blocks")
    )$col_names
  }
  out
}

.store_write_parquet <- function(df, path) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required for Parquet storage.", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(df, sink = path, compression = "zstd")
  invisible(path)
}

.store_read_parquet <- function(path, columns = NULL) {
  if (!file.exists(path)) {
    return(NULL)
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required for Parquet storage.", call. = FALSE)
  }
  tab <- as.data.frame(arrow::read_parquet(file = path))
  if (!is.null(columns)) {
    keep <- intersect(columns, names(tab))
    tab <- tab[, keep, drop = FALSE]
  }
  tab
}

.store_sqlite_con <- function(object, create = TRUE) {
  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("RSQLite", quietly = TRUE)) {
    stop("Packages 'DBI' and 'RSQLite' are required for SQLite storage.",
         call. = FALSE)
  }
  db_path <- .store_path(object)
  DBI::dbConnect(RSQLite::SQLite(), dbname = db_path, create = create)
}

.store_write_table <- function(object, df, table_id) {
  if (.store_mode(object) == "parquet") {
    .store_write_parquet(df, table_id)
  } else {
    con <- .store_sqlite_con(object)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    DBI::dbWriteTable(con, name = table_id, value = df, overwrite = TRUE)
  }
  invisible(NULL)
}

.store_read_table <- function(object, table_id, columns = NULL) {
  if (.store_mode(object) == "parquet") {
    return(.store_read_parquet(table_id, columns = columns))
  }
  if (!file.exists(.store_path(object))) {
    return(NULL)
  }
  con <- .store_sqlite_con(object, create = FALSE)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  if (!table_id %in% DBI::dbListTables(con)) {
    return(NULL)
  }
  if (is.null(columns)) {
    return(DBI::dbReadTable(con, name = table_id))
  }
  cols <- paste(columns, collapse = ", ")
  DBI::dbGetQuery(con, paste("SELECT", cols, "FROM", DBI::dbQuoteIdentifier(con, table_id)))
}

.store_table_exists <- function(object, table_id) {
  if (.store_is_gds(object)) {
    return(FALSE)
  }
  if (.store_mode(object) == "parquet") {
    return(file.exists(table_id))
  }
  if (!file.exists(.store_path(object))) {
    return(FALSE)
  }
  con <- .store_sqlite_con(object, create = FALSE)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  table_id %in% DBI::dbListTables(con)
}

.store_check_scan_exists <- function(object) {
  if (.store_is_gds(object)) {
    return(exist.gdsn(node = object$root, path = "lazygas/scan"))
  }
  meta <- .store_read_meta(object)
  length(meta$scan_phenotypes) > 0L ||
    dir.exists(file.path(.store_path(object), "scan"))
}

.store_scan_pheno_exists <- function(object, pheno_name) {
  if (.store_is_gds(object)) {
    return(exist.gdsn(node = object$root,
                      path = paste0("lazygas/scan/", pheno_name)))
  }
  .store_table_exists(object, .store_scan_path(object, pheno_name))
}

.store_list_scan_phenos <- function(object) {
  if (.store_is_gds(object)) {
    if (!exist.gdsn(node = object$root, path = "lazygas/scan")) {
      return(character())
    }
    nodes <- gdsfmt::ls.gdsn(node = gdsfmt::index.gdsn(node = object, path = "lazygas/scan"))
    meta_nodes <- c("kruskal", "formula", "null_formula", "conv_fun",
                    "geno_format", "fixed_effect")
    return(setdiff(nodes, meta_nodes))
  }
  meta <- .store_read_meta(object)
  if (length(meta$scan_phenotypes) > 0L) {
    return(unlist(meta$scan_phenotypes, use.names = FALSE))
  }
  scan_dir <- file.path(.store_path(object), "scan")
  if (!dir.exists(scan_dir)) {
    return(character())
  }
  files <- list.files(scan_dir, pattern = "\\.parquet$", full.names = FALSE)
  sub("\\.parquet$", "", files)
}

.store_write_scan <- function(object, pheno_name, mat, colnames_stats) {
  if (.store_is_gds(object)) {
  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = pheno_name,
               val = mat,
               storage = "float",
               valdim = dim(mat),
               attr = list(colnames = colnames_stats))
    return(invisible(NULL))
  }

  df <- data.frame(mat, check.names = FALSE)
  names(df) <- colnames_stats
  path <- .store_scan_path(object, pheno_name)
  .store_write_table(object, df, path)
  phenos <- unique(c(.store_list_scan_phenos(object), pheno_name))
  meta <- .store_read_meta(object)
  meta$scan_phenotypes <- phenos
  meta$scan_colnames[[pheno_name]] <- colnames_stats
  .store_write_meta(object, meta)
  invisible(NULL)
}

.store_read_scan_matrix <- function(object, pheno_name, columns = NULL) {
  if (.store_is_gds(object)) {
    scan_node <- paste0("lazygas/scan/", pheno_name)
    col_names <- gdsfmt::get.attr.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = scan_node)
    )$colnames
    if (is.null(columns)) {
      mat <- gdsfmt::read.gdsn(
        node = gdsfmt::index.gdsn(node = object$root, path = scan_node)
      )
    } else {
      dim <- gdsfmt::objdesp.gdsn(
        node = gdsfmt::index.gdsn(node = object$root, path = scan_node)
      )
      sel <- list(rep(TRUE, dim$dim[1]), col_names %in% columns)
      mat <- gdsfmt::readex.gdsn(
        node = gdsfmt::index.gdsn(node = object$root, path = scan_node),
        sel = sel
      )
      col_names <- col_names[sel[[2]]]
    }
    return(list(matrix = mat, colnames = col_names))
  }

  path <- .store_scan_path(object, pheno_name)
  if (!.store_table_exists(object, path)) {
    return(NULL)
  }
  df <- .store_read_table(object, path, columns = columns)
  list(matrix = as.matrix(df), colnames = names(df))
}

.store_write_scan_meta <- function(object,
                                   kruskal,
                                   formula,
                                   null_formula,
                                   fixed_effect,
                                   conv_fun,
                                   geno_format) {
  if (.store_is_gds(object)) {
    .store_additional_info_gds(object, kruskal, formula, null_formula,
                             fixed_effect, conv_fun, geno_format)
    return(invisible(NULL))
  }

  if (!is.character(conv_fun)) {
    conv_fun <- deparse(conv_fun)
  }
  if (!is.character(formula)) {
    formula <- deparse(formula)
  }
  if (is.null(null_formula)) {
    null_formula <- "NULL"
  }
  fe <- fixed_effect
  fe_col <- NULL
  if (is.null(fixed_effect)) {
    fe <- matrix(NA_real_, nrow = 1, ncol = 1)
  } else {
    fe_col <- names(fixed_effect)
    fe <- as.matrix(fixed_effect)
  }

  meta <- .store_read_meta(object)
  meta$kruskal <- as.character(kruskal)
  meta$formula <- as.character(formula)
  meta$null_formula <- as.character(null_formula)
  meta$conv_fun <- as.character(conv_fun)
  meta$geno_format <- as.character(geno_format)
  meta$fixed_effect <- fe
  meta$fixed_effect_col_names <- fe_col
  .store_write_meta(object, meta)
  invisible(NULL)
}

.store_read_scan_scalar <- function(object, name) {
  if (.store_is_gds(object)) {
    path <- paste0("lazygas/scan/", name)
    if (!exist.gdsn(node = object$root, path = path)) {
      return(NULL)
    }
    return(.get_data_gds_only(object, path))
  }
  meta <- .store_read_meta(object)
  meta[[name]]
}

.store_read_fixed_effect <- function(object) {
  if (.store_is_gds(object)) {
    fe <- .get_data_gds_only(object, "lazygas/scan/fixed_effect")
    if (length(fe) == 1L) {
      return(NULL)
    }
    fe <- as.data.frame(fe)
    col_names <- gdsfmt::get.attr.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = "lazygas/scan/fixed_effect")
    )
    colnames(fe) <- col_names
    return(fe)
  }
  meta <- .store_read_meta(object)
  fe <- meta$fixed_effect
  if (is.null(fe) || (length(fe) == 1L && is.na(fe[[1]]))) {
    return(NULL)
  }
  fe <- as.data.frame(fe)
  if (!is.null(meta$fixed_effect_col_names)) {
    colnames(fe) <- meta$fixed_effect_col_names
  }
  fe
}

.store_ensure_scan_folder <- function(object) {
  if (.store_is_gds(object)) {
    if (!exist.gdsn(node = object$root, path = "lazygas")) {
      .create_gdsn(root_node = object$root, target_node = "",
                   new_node = "lazygas", is_folder = TRUE)
    }
    if (!exist.gdsn(node = object$root, path = "lazygas/scan")) {
      .create_gdsn(root_node = object$root, target_node = "lazygas",
                   new_node = "scan", is_folder = TRUE)
    }
  }
  invisible(NULL)
}

.store_set_peakcall_params <- function(object, signif, threshold) {
  if (.store_is_gds(object)) {
    return(invisible(NULL))
  }
  .store_update_meta(object,
                     peakcall_signif = signif,
                     peakcall_threshold = threshold,
                     peakcall_peaks_col_names = c("peak_ID", "peak_variant_ID"),
                     peakcall_blocks_col_names = c(
                       "peak_ID", "variant_ID", "dist2peak", "LD2peak",
                       "chr_start", "chr_end"
                     ))
}

.store_peak_buffer_key <- function(section, pheno_name, kind) {
  paste(section, pheno_name, kind, sep = "|")
}

.store_peak_buffer_init <- function(object, section, pheno_name) {
  env <- object@store$peak_buffer
  assign(.store_peak_buffer_key(section, pheno_name, "peaks"),
         data.frame(peak_ID = integer(), peak_variant_ID = character()),
         envir = env)
  assign(.store_peak_buffer_key(section, pheno_name, "blocks"),
         data.frame(
           peak_ID = integer(),
           variant_ID = character(),
           dist2peak = numeric(),
           LD2peak = numeric(),
           chr_start = integer(),
           chr_end = integer()
         ),
         envir = env)
}

.store_peak_buffer_append <- function(object, section, pheno_name,
                                      peak_id, peak_info, peak_block) {
  env <- object@store$peak_buffer
  peaks <- get(.store_peak_buffer_key(section, pheno_name, "peaks"), envir = env)
  blocks <- get(.store_peak_buffer_key(section, pheno_name, "blocks"), envir = env)

  peaks <- rbind(
    peaks,
    data.frame(
      peak_ID = peak_id,
      peak_variant_ID = peak_info$variantID,
      stringsAsFactors = FALSE
    )
  )

  n_var <- length(peak_block$ids)
  blocks <- rbind(
    blocks,
    data.frame(
      peak_ID = rep(peak_id, n_var),
      variant_ID = peak_block$ids,
      dist2peak = peak_block$dist,
      LD2peak = peak_block$ld,
      chr_start = rep(peak_block$chr_start, n_var),
      chr_end = rep(peak_block$chr_end, n_var),
      stringsAsFactors = FALSE
    )
  )

  assign(.store_peak_buffer_key(section, pheno_name, "peaks"), peaks, envir = env)
  assign(.store_peak_buffer_key(section, pheno_name, "blocks"), blocks, envir = env)
}

.store_peak_buffer_flush <- function(object, section, pheno_name) {
  if (.store_is_gds(object)) {
    return(invisible(NULL))
  }
  env <- object@store$peak_buffer
  pk <- .store_peak_buffer_key(section, pheno_name, "peaks")
  bk <- .store_peak_buffer_key(section, pheno_name, "blocks")
  if (!exists(pk, envir = env, inherits = FALSE)) {
    return(invisible(NULL))
  }
  peaks <- get(pk, envir = env)
  blocks <- get(bk, envir = env)
  .store_write_peaks_df(object, section, pheno_name, "peaks", peaks)
  .store_write_peaks_df(object, section, pheno_name, "blocks", blocks)
  invisible(NULL)
}

.store_write_peaks_df <- function(object, section, pheno_name, kind, df) {
  if (.store_is_gds(object)) {
    return(invisible(NULL))
  }
  path <- .store_table_path(object, section, kind, pheno_name)
  if (nrow(df) == 0L) {
    meta <- .store_read_meta(object)
    col_key <- paste0(section, "_", kind, "_col_names")
    if (!is.null(meta[[col_key]])) {
      df <- as.data.frame(
        setNames(
          rep(list(if (kind == "peaks") integer() else numeric()), length(meta[[col_key]])),
          meta[[col_key]]
        )
      )
    }
  }
  .store_write_table(object, df, path)
  invisible(NULL)
}

.store_write_empty_peak_tables <- function(object, section, pheno_name) {
  if (.store_is_gds(object)) {
    return(invisible(NULL))
  }
  meta <- .store_read_meta(object)
  if (section == "peakcall") {
    pcols <- meta$peakcall_peaks_col_names
    bcols <- meta$peakcall_blocks_col_names
  } else {
    pcols <- c("peak_ID", "peak_variant_ID")
    bcols <- c("peak_ID", "variant_ID", "dist2peak", "LD2peak", "chr_start", "chr_end")
  }
  .store_write_peaks_df(
    object, section, pheno_name, "peaks",
    as.data.frame(setNames(rep(list(integer()), length(pcols)), pcols))
  )
  .store_write_peaks_df(
    object, section, pheno_name, "blocks",
    as.data.frame(setNames(rep(list(numeric()), length(bcols)), bcols))
  )
}

.store_read_peaks_df <- function(object, section, pheno_name, kind) {
  if (.store_is_gds(object)) {
    return(NULL)
  }
  path <- .store_table_path(object, section, kind, pheno_name)
  .store_read_table(object, path)
}

.store_section_exists <- function(object, section) {
  if (.store_is_gds(object)) {
    return(exist.gdsn(node = object$root, path = paste0("lazygas/", section)))
  }
  dir.exists(file.path(.store_path(object), section))
}

.store_pheno_table_exists <- function(object, section, kind, pheno_name) {
  if (.store_is_gds(object)) {
    return(exist.gdsn(
      node = object$root,
      path = paste0("lazygas/", section, "/", kind, "/", pheno_name)
    ))
  }
  .store_table_exists(object, .store_table_path(object, section, kind, pheno_name))
}

.store_write_recalc <- function(object, pheno_name, peaks, blocks, groups) {
  if (.store_is_gds(object)) {
    return(invisible(NULL))
  }
  meta <- .store_read_meta(object)
  .store_write_peaks_df(object, "recalc", pheno_name, "peaks", peaks)
  .store_write_peaks_df(object, "recalc", pheno_name, "blocks", blocks)
  path <- .store_table_path(object, "recalc", "groups", pheno_name)
  .store_write_table(object, groups, path)
  meta$recalc_groups_col_names <- names(groups)
  .store_write_meta(object, meta)
  invisible(NULL)
}

.store_write_recalc_empty <- function(object, pheno_name) {
  if (.store_is_gds(object)) {
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas/recalc/peaks",
                 new_node = pheno_name,
                 storage = "uint32",
                 replace = TRUE)
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas/recalc/blocks",
                 new_node = pheno_name,
                 storage = "double",
                 replace = TRUE)
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas/recalc/groups",
                 new_node = pheno_name,
                 storage = "double",
                 replace = TRUE)
    return(invisible(NULL))
  }
  .store_write_empty_peak_tables(object, "recalc", pheno_name)
  empty_groups <- data.frame(
    peak_ID = integer(),
    variant_ID = character(),
    stringsAsFactors = FALSE
  )
  path <- .store_table_path(object, "recalc", "groups", pheno_name)
  .store_write_table(object, empty_groups, path)
  invisible(NULL)
}

.store_save_recalc_results <- function(object, pheno_name, out1, out2, out3) {
  if (.store_is_gds(object)) {
    out1_gdsn <- .create_gdsn(
      root_node = object$root,
      target_node = "lazygas/recalc/peaks",
      new_node = pheno_name,
      val = as.matrix(out1),
      storage = "uint32",
      replace = TRUE
    )
    out2_gdsn <- .create_gdsn(
      root_node = object$root,
      target_node = "lazygas/recalc/blocks",
      new_node = pheno_name,
      val = as.matrix(out2),
      storage = "double",
      replace = TRUE
    )
    out3_gdsn <- .create_gdsn(
      root_node = object$root,
      target_node = "lazygas/recalc/groups",
      new_node = pheno_name,
      val = as.matrix(out3),
      storage = "double",
      replace = TRUE
    )
    gdsfmt::put.attr.gdsn(node = out1_gdsn, name = "col_names", val = names(out1))
    gdsfmt::put.attr.gdsn(node = out2_gdsn, name = "col_names", val = names(out2))
    gdsfmt::put.attr.gdsn(node = out3_gdsn, name = "col_names", val = names(out3))
    return(invisible(NULL))
  }
  .store_write_recalc(object, pheno_name, out1, out2, out3)
}

.store_dataset_exists <- function(object, dataset, pheno_name = NULL) {
  if (.store_is_gds(object)) {
    if (dataset == "scan") {
      if (is.null(pheno_name)) {
        return(exist.gdsn(node = object$root, path = "lazygas/scan"))
      }
      return(exist.gdsn(node = object$root,
                        path = paste0("lazygas/scan/", pheno_name)))
    }
    if (dataset == "groups") {
      return(exist.gdsn(node = object$root,
                        path = paste0("lazygas/recalc/groups/", pheno_name)))
    }
    if (dataset %in% c("candidate", "snpeff")) {
      return(exist.gdsn(node = object$root,
                        path = paste0("lazygas/", dataset, "/", pheno_name)))
    }
    return(exist.gdsn(
      node = object$root,
      path = paste0("lazygas/", dataset, "/peaks/", pheno_name)
    ))
  }
  if (dataset == "scan") {
    if (is.null(pheno_name)) {
      return(.store_check_scan_exists(object))
    }
    return(.store_scan_pheno_exists(object, pheno_name))
  }
  if (dataset == "groups") {
    return(.store_pheno_table_exists(object, "recalc", "groups", pheno_name))
  }
  if (dataset == "candidate") {
    return(.store_pheno_table_exists(object, "candidate", "candidate", pheno_name))
  }
  if (dataset == "snpeff") {
    return(.store_pheno_table_exists(object, "snpeff", "snpeff", pheno_name))
  }
  .store_pheno_table_exists(object, dataset, "peaks", pheno_name)
}

#' Import lazyGas results from the legacy GDS `lazygas/` subtree into the companion store.
#' @param object A LazyGas object using Parquet or SQLite storage.
#' @export
importLazyGasResults <- function(object) {
  if (.store_is_gds(object)) {
    stop("object already uses GDS storage (lazygas_store = \"gds\").", call. = FALSE)
  }
  if (!exist.gdsn(node = object$root, path = "lazygas")) {
    message("No lazygas/ node in the GDS file.")
    return(invisible(object))
  }
  .store_migrate_from_gds(object)
  invisible(object)
}

.store_write_candidate <- function(object, pheno_name, candidate, snpeff) {
  if (.store_is_gds(object)) {
    if (is.null(candidate)) {
      .create_gdsn(root_node = object$root,
                   target_node = "lazygas/candidate",
                   new_node = pheno_name,
                   storage = "string",
                   replace = TRUE)
    } else {
      .create_gdsn(root_node = object$root,
                   target_node = "lazygas/candidate",
                   new_node = pheno_name,
                   val = as.matrix(candidate),
                   storage = "string",
                   replace = TRUE)
      gdsfmt::put.attr.gdsn(
        node = gdsfmt::index.gdsn(node = object$root, path = "lazygas/candidate"),
        name = "col_names",
        val = colnames(candidate)
      )
    }
    if (is.null(snpeff) || (is.data.frame(snpeff) && nrow(snpeff) == 0L)) {
      .create_gdsn(root_node = object$root,
                   target_node = "lazygas/snpeff",
                   new_node = pheno_name,
                   storage = "string",
                   replace = TRUE)
    } else {
      .create_gdsn(root_node = object$root,
                   target_node = "lazygas/snpeff",
                   new_node = pheno_name,
                   val = as.matrix(snpeff),
                   storage = "string",
                   replace = TRUE)
      gdsfmt::put.attr.gdsn(
        node = gdsfmt::index.gdsn(node = object$root, path = "lazygas/snpeff"),
        name = "col_names",
        val = colnames(snpeff)
      )
    }
    return(invisible(NULL))
  }
  meta <- .store_read_meta(object)
  path_c <- .store_table_path(object, "candidate", "candidate", pheno_name)
  path_s <- .store_table_path(object, "snpeff", "snpeff", pheno_name)
  if (is.null(candidate)) {
    .store_write_table(object, data.frame(), path_c)
  } else {
    .store_write_table(object, as.data.frame(candidate), path_c)
    meta$candidate_col_names <- colnames(candidate)
  }
  if (is.null(snpeff) || (is.data.frame(snpeff) && nrow(snpeff) == 0L)) {
    .store_write_table(object, data.frame(), path_s)
  } else {
    .store_write_table(object, as.data.frame(snpeff), path_s)
    meta$snpeff_col_names <- colnames(snpeff)
  }
  .store_write_meta(object, meta)
  invisible(NULL)
}

.store_read_matrix_dataset <- function(object, dataset, pheno_name) {
  if (.store_is_gds(object)) {
    return(NULL)
  }
  if (dataset == "groups") {
    path <- .store_table_path(object, "recalc", "groups", pheno_name)
    df <- .store_read_table(object, path)
    if (is.null(df) || nrow(df) == 0L) {
      return(NULL)
    }
    meta <- .store_read_meta(object)
    colnames(df) <- meta$recalc_groups_col_names
    return(df)
  }
  if (dataset == "candidate") {
    path <- .store_table_path(object, "candidate", "candidate", pheno_name)
    df <- .store_read_table(object, path)
    if (is.null(df)) {
      return(NULL)
    }
    meta <- .store_read_meta(object)
    if (!is.null(meta$candidate_col_names)) {
      colnames(df) <- meta$candidate_col_names
    }
    return(df)
  }
  if (dataset == "snpeff") {
    path <- .store_table_path(object, "snpeff", "snpeff", pheno_name)
    df <- .store_read_table(object, path)
    if (is.null(df)) {
      return(NULL)
    }
    meta <- .store_read_meta(object)
    if (!is.null(meta$snpeff_col_names)) {
      colnames(df) <- meta$snpeff_col_names
    }
    return(df)
  }
  NULL
}

.store_peak_meta <- function(object, section, pheno_name) {
  if (.store_is_gds(object)) {
    if (section == "recalc") {
      peaks_gdsn <- "lazygas/recalc/peaks"
      blocks_gdsn <- "lazygas/recalc/blocks"
      peaks_col <- gdsfmt::get.attr.gdsn(
        node = gdsfmt::index.gdsn(
          node = object$root,
          path = paste0(peaks_gdsn, "/", pheno_name)
        )
      )
      blocks_col <- gdsfmt::get.attr.gdsn(
        node = gdsfmt::index.gdsn(
          node = object$root,
          path = paste0(blocks_gdsn, "/", pheno_name)
        )
      )
    } else {
      peaks_col <- gdsfmt::get.attr.gdsn(
        node = gdsfmt::index.gdsn(node = object$root, path = "lazygas/peakcall/peaks")
      )
      blocks_col <- gdsfmt::get.attr.gdsn(
        node = gdsfmt::index.gdsn(node = object$root, path = "lazygas/peakcall/blocks")
      )
    }
    return(list(
      peaks_col_names = peaks_col$col_names,
      blocks_col_names = blocks_col$col_names,
      groups_col_names = NULL,
      signif = peaks_col$signif,
      threshold = peaks_col$threshold
    ))
  }
  meta <- .store_read_meta(object)
  list(
    peaks_col_names = meta$peakcall_peaks_col_names,
    blocks_col_names = meta$peakcall_blocks_col_names,
    groups_col_names = meta$recalc_groups_col_names,
    signif = meta$peakcall_signif,
    threshold = meta$peakcall_threshold
  )
}

.store_migrate_from_gds <- function(object) {
  if (!exist.gdsn(node = object$root, path = "lazygas")) {
    return(invisible(NULL))
  }
  if (.store_is_gds(object)) {
    return(invisible(NULL))
  }
  message("Migrating lazygas/ results from GDS to companion storage...")
  meta <- .store_read_meta_gds(object)
  .store_write_meta(object, meta)

  if (exist.gdsn(node = object$root, path = "lazygas/scan")) {
    phenos <- .store_list_scan_phenos_gds(object)
    for (pheno in phenos) {
      scan_node <- paste0("lazygas/scan/", pheno)
      col_names <- gdsfmt::get.attr.gdsn(
        node = gdsfmt::index.gdsn(node = object$root, path = scan_node)
      )$colnames
      mat <- gdsfmt::read.gdsn(
        node = gdsfmt::index.gdsn(node = object$root, path = scan_node)
      )
      .store_write_scan(object, pheno, mat, col_names)
    }
  }

  for (section in c("peakcall", "recalc")) {
    for (kind in c("peaks", "blocks", "groups")) {
      base <- paste0("lazygas/", section, "/", kind)
      if (!exist.gdsn(node = object$root, path = base)) {
        next
      }
      phenos <- gdsfmt::ls.gdsn(
        node = gdsfmt::index.gdsn(node = object, path = base)
      )
      for (pheno in phenos) {
        dfs <- .store_gds_peak_kind_to_df(
          object = object,
          section = section,
          kind = kind,
          pheno_name = pheno
        )
        if (is.null(dfs)) {
          next
        }
        .store_write_peaks_df(object, section, pheno, kind, dfs)
      }
    }
  }

  for (kind in c("candidate", "snpeff")) {
    base <- paste0("lazygas/", kind)
    if (!exist.gdsn(node = object$root, path = base)) {
      next
    }
    phenos <- gdsfmt::ls.gdsn(
      node = gdsfmt::index.gdsn(node = object, path = base)
    )
    att <- gdsfmt::get.attr.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = base)
    )
    for (pheno in phenos) {
      raw <- .get_data_gds_only(object, paste0(base, "/", pheno))
      if (length(raw) == 0L) {
        next
      }
      df <- as.data.frame(raw, stringsAsFactors = FALSE)
      if (length(att$col_names)) {
        colnames(df) <- att$col_names
      }
      if (kind == "candidate") {
        .store_write_candidate(object, pheno, df, NULL)
      } else {
        .store_write_candidate(object, pheno, NULL, df)
      }
    }
  }
  invisible(NULL)
}

.store_gds_peak_kind_to_df <- function(object, section, kind, pheno_name) {
  node <- paste0("lazygas/", section, "/", kind, "/", pheno_name)
  if (!exist.gdsn(node = object$root, path = node)) {
    return(NULL)
  }
  raw <- .get_data_gds_only(object, node)
  if (length(raw) == 0L) {
    return(NULL)
  }
  if (!is.matrix(raw)) {
    raw <- matrix(raw, nrow = 1)
  }
  if (section == "recalc") {
    df <- data.frame(raw, stringsAsFactors = FALSE)
  } else {
    df <- data.frame(t(raw), stringsAsFactors = FALSE)
  }
  att_path <- if (section == "recalc") {
    node
  } else {
    paste0("lazygas/", section, "/", kind)
  }
  att <- gdsfmt::get.attr.gdsn(
    node = gdsfmt::index.gdsn(node = object$root, path = att_path)
  )
  if (length(att$col_names)) {
    colnames(df) <- att$col_names
  }
  df
}

.store_list_scan_phenos_gds <- function(object) {
  nodes <- gdsfmt::ls.gdsn(node = gdsfmt::index.gdsn(node = object, path = "lazygas/scan"))
  meta_nodes <- c("kruskal", "formula", "null_formula", "conv_fun",
                  "geno_format", "fixed_effect")
  setdiff(nodes, meta_nodes)
}

.get_data_gds_only <- function(object, node, sel = NULL) {
  if (is.null(sel)) {
    gdsfmt::read.gdsn(node = gdsfmt::index.gdsn(node = object$root, path = node))
  } else {
    gdsfmt::readex.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = node),
      sel = sel
    )
  }
}

.store_additional_info_gds <- function(object,
                                       kruskal,
                                       formula,
                                       null_formula,
                                       fixed_effect,
                                       conv_fun,
                                       geno_format) {
  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = "kruskal",
               val = kruskal,
               storage = "string32",
               valdim = dim(kruskal))

  if (!is.character(conv_fun)) {
    conv_fun <- deparse(conv_fun)
  }
  if (!is.character(formula)) {
    formula <- deparse(formula)
  }

  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = "formula",
               val = formula,
               storage = "string32",
               valdim = dim(formula))

  if (is.null(null_formula)) {
    null_formula <- "NULL"
  }

  if (is.null(fixed_effect)) {
    fixed_effect <- matrix(data = NA, nrow = 1, ncol = 1)
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas/scan",
                 new_node = "fixed_effect",
                 val = fixed_effect,
                 storage = "float",
                 valdim = dim(fixed_effect))
  } else {
    fixed_effect_attr <- names(fixed_effect)
    fixed_effect <- as.matrix(fixed_effect)
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas/scan",
                 new_node = "fixed_effect",
                 val = fixed_effect,
                 storage = "float",
                 valdim = dim(fixed_effect),
                 attr = list(col_names = fixed_effect_attr))
  }

  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = "null_formula",
               val = null_formula,
               storage = "string32",
               valdim = dim(null_formula))

  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = "conv_fun",
               val = conv_fun,
               storage = "string32",
               valdim = dim(deparse(conv_fun)))

  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = "geno_format",
               val = geno_format,
               storage = "string32",
               valdim = dim(geno_format))
}

.store_close <- function(object) {
  if (.store_mode(object) == "sqlite" && !is.null(object@store$db_con)) {
    DBI::dbDisconnect(object@store$db_con)
    object@store$db_con <- NULL
  }
  invisible(object)
}
