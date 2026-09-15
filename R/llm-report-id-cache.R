# Phase 2 E5–E7: official ID → name caches (not bundled; download once)

#' Cache root for GO / KEGG / Pfam / InterPro label maps
#' @keywords internal
.lazygas_id_cache_dir <- function(cache_dir = NULL) {
  if (!is.null(cache_dir) && nzchar(as.character(cache_dir)[1L])) {
    d <- path.expand(as.character(cache_dir)[1L])
  } else {
    d <- file.path(tools::R_user_dir("lazyGas", which = "cache"), "id_maps")
  }
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Download a URL to a local file (text or binary)
#' @keywords internal
.lazygas_download_file <- function(url, dest, timeout = 300, mode = "wb") {
  dest <- path.expand(dest)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  ok <- tryCatch(
    {
      if (requireNamespace("curl", quietly = TRUE)) {
        curl::curl_download(url, destfile = dest, quiet = TRUE, mode = mode)
        TRUE
      } else {
        utils::download.file(
          url, destfile = dest, quiet = TRUE, mode = mode,
          timeout = timeout
        )
        TRUE
      }
    },
    error = function(e) FALSE
  )
  isTRUE(ok) && file.exists(dest) && file.info(dest)$size > 0
}

#' Parse go-basic.obo into ID → name named character
#' @keywords internal
.parse_go_obo_id_name <- function(path) {
  lines <- readLines(path, warn = FALSE)
  ids <- character()
  names_v <- character()
  cur_id <- NA_character_
  cur_name <- NA_character_
  in_term <- FALSE
  flush <- function() {
    if (isTRUE(in_term) && !is.na(cur_id) && nzchar(cur_id) &&
        !is.na(cur_name) && nzchar(cur_name)) {
      ids <<- c(ids, cur_id)
      names_v <<- c(names_v, cur_name)
    }
  }
  for (ln in lines) {
    if (identical(ln, "[Term]")) {
      flush()
      in_term <- TRUE
      cur_id <- NA_character_
      cur_name <- NA_character_
      next
    }
    if (startsWith(ln, "[")) {
      flush()
      in_term <- FALSE
      next
    }
    if (!isTRUE(in_term)) next
    if (startsWith(ln, "id: ")) {
      cur_id <- sub("^id:\\s*", "", ln)
      cur_id <- .normalize_go_id(cur_id)
    } else if (startsWith(ln, "name: ")) {
      cur_name <- sub("^name:\\s*", "", ln)
    }
  }
  flush()
  if (!length(ids)) {
    return(character())
  }
  # first wins
  keep <- !duplicated(ids)
  setNames(names_v[keep], ids[keep])
}

#' @keywords internal
.normalize_go_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^GO[:_ ]*", "GO:", x)
  # pad digits
  m <- regexec("^GO:([0-9]+)$", x)
  hit <- regmatches(x, m)
  vapply(seq_along(x), function(i) {
    h <- hit[[i]]
    if (length(h) >= 2L) {
      paste0("GO:", sprintf("%07d", as.integer(h[2L])))
    } else {
      x[i]
    }
  }, character(1L))
}

#' Load GO ID→name map (cache + official OBO)
#'
#' @param cache_dir Optional cache directory.
#' @param obo_path Optional local OBO path (tests / offline).
#' @param download If \code{TRUE} and cache missing, fetch go-basic.obo.
#' @return Named character vector (names = GO IDs).
#' @keywords internal
.go_term_name_map <- function(cache_dir = NULL,
                              obo_path = NULL,
                              download = TRUE) {
  cache <- .lazygas_id_cache_dir(cache_dir)
  rds <- file.path(cache, "go_id_name.rds")
  if (file.exists(rds) && is.null(obo_path)) {
    return(readRDS(rds))
  }
  obo <- obo_path
  if (is.null(obo) || !file.exists(obo)) {
    obo <- file.path(cache, "go-basic.obo")
  }
  if (!file.exists(obo)) {
    if (!isTRUE(download)) {
      return(character())
    }
    url <- "http://current.geneontology.org/ontology/go-basic.obo"
    if (!.lazygas_download_file(url, obo, mode = "wb")) {
      warning("Failed to download GO OBO; unresolved GO IDs will be dropped.",
              call. = FALSE)
      return(character())
    }
  }
  map <- .parse_go_obo_id_name(obo)
  saveRDS(map, rds)
  map
}

#' Parse KEGG REST list (tab-separated id\\tname)
#' @keywords internal
.parse_kegg_list_text <- function(text) {
  text <- as.character(text)
  if (!length(text)) {
    return(character())
  }
  if (length(text) == 1L && grepl("\n", text, fixed = TRUE)) {
    text <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  }
  text <- text[nzchar(trimws(text))]
  ids <- character()
  nms <- character()
  for (ln in text) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 2L) next
    id <- trimws(parts[1L])
    nm <- trimws(parts[2L])
    # strip path: prefix e.g. path:map00500 → map00500
    id <- sub("^(path|ko|ec):", "", id)
    if (!nzchar(id) || !nzchar(nm)) next
    ids <- c(ids, id)
    nms <- c(nms, nm)
  }
  if (!length(ids)) {
    return(character())
  }
  keep <- !duplicated(ids)
  setNames(nms[keep], ids[keep])
}

#' Load KEGG list map for pathway / ko / enzyme
#' @keywords internal
.kegg_list_name_map <- function(kind = c("pathway", "ko", "enzyme"),
                                org = NULL,
                                cache_dir = NULL,
                                text = NULL,
                                download = TRUE) {
  kind <- match.arg(kind)
  cache <- .lazygas_id_cache_dir(cache_dir)
  tag <- if (!is.null(org) && nzchar(org) && identical(kind, "pathway")) {
    paste0("pathway_", org)
  } else {
    kind
  }
  rds <- file.path(cache, paste0("kegg_", tag, ".rds"))
  if (file.exists(rds) && is.null(text)) {
    return(readRDS(rds))
  }
  if (!is.null(text)) {
    map <- .parse_kegg_list_text(text)
    saveRDS(map, rds)
    return(map)
  }
  if (!isTRUE(download)) {
    return(character())
  }
  url <- if (identical(kind, "pathway") && !is.null(org) && nzchar(org)) {
    paste0("https://rest.kegg.jp/list/pathway/", org)
  } else if (identical(kind, "pathway")) {
    "https://rest.kegg.jp/list/pathway"
  } else if (identical(kind, "ko")) {
    "https://rest.kegg.jp/list/ko"
  } else {
    "https://rest.kegg.jp/list/enzyme"
  }
  raw <- tryCatch(
    {
      if (requireNamespace("curl", quietly = TRUE)) {
        rawToChar(curl::curl_fetch_memory(url)$content)
      } else {
        paste(readLines(url, warn = FALSE), collapse = "\n")
      }
    },
    error = function(e) NULL
  )
  if (is.null(raw) || !nzchar(raw)) {
    warning("Failed to download KEGG list '", kind, "'; unresolved IDs dropped.",
            call. = FALSE)
    return(character())
  }
  map <- .parse_kegg_list_text(raw)
  saveRDS(map, rds)
  map
}

#' Parse Pfam-A.clans.tsv (or minimal 2-col) accession → description
#' @keywords internal
.parse_pfam_clans_tsv <- function(path) {
  # Pfam-A.clans.tsv.gz columns: accession, id, clan, clan_id, description...
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  ids <- character()
  nms <- character()
  for (ln in lines) {
    if (!nzchar(ln) || startsWith(ln, "#")) next
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 2L) next
    acc <- trimws(parts[1L])
    # description often last column; prefer col 5 if present else col 2
    desc <- if (length(parts) >= 5L) trimws(parts[5L]) else trimws(parts[2L])
    if (!nzchar(acc) || !grepl("^PF", acc)) next
    if (!nzchar(desc)) next
    ids <- c(ids, acc)
    nms <- c(nms, desc)
  }
  if (!length(ids)) {
    return(character())
  }
  keep <- !duplicated(ids)
  setNames(nms[keep], ids[keep])
}

#' Parse InterPro entry.list
#' @keywords internal
.parse_interpro_entry_list <- function(path) {
  lines <- readLines(path, warn = FALSE)
  ids <- character()
  nms <- character()
  for (ln in lines) {
    if (!nzchar(ln) || startsWith(ln, "ENTRY") || startsWith(ln, "#")) next
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 2L) next
    acc <- trimws(parts[1L])
    nm <- trimws(parts[length(parts)])
    if (!grepl("^IPR", acc) || !nzchar(nm)) next
    ids <- c(ids, acc)
    nms <- c(nms, nm)
  }
  if (!length(ids)) {
    return(character())
  }
  keep <- !duplicated(ids)
  setNames(nms[keep], ids[keep])
}

#' @keywords internal
.pfam_name_map <- function(cache_dir = NULL,
                           tsv_path = NULL,
                           download = TRUE) {
  cache <- .lazygas_id_cache_dir(cache_dir)
  rds <- file.path(cache, "pfam_id_name.rds")
  if (file.exists(rds) && is.null(tsv_path)) {
    return(readRDS(rds))
  }
  path <- tsv_path
  if (is.null(path) || !file.exists(path)) {
    path <- file.path(cache, "Pfam-A.clans.tsv.gz")
  }
  if (!file.exists(path)) {
    if (!isTRUE(download)) {
      return(character())
    }
    url <- "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"
    if (!.lazygas_download_file(url, path, mode = "wb")) {
      warning("Failed to download Pfam clans TSV; unresolved Pfam IDs dropped.",
              call. = FALSE)
      return(character())
    }
  }
  map <- .parse_pfam_clans_tsv(path)
  saveRDS(map, rds)
  map
}

#' @keywords internal
.interpro_name_map <- function(cache_dir = NULL,
                               list_path = NULL,
                               download = TRUE) {
  cache <- .lazygas_id_cache_dir(cache_dir)
  rds <- file.path(cache, "interpro_id_name.rds")
  if (file.exists(rds) && is.null(list_path)) {
    return(readRDS(rds))
  }
  path <- list_path
  if (is.null(path) || !file.exists(path)) {
    path <- file.path(cache, "interpro_entry.list")
  }
  if (!file.exists(path)) {
    if (!isTRUE(download)) {
      return(character())
    }
    url <- "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list"
    if (!.lazygas_download_file(url, path, mode = "wb")) {
      warning("Failed to download InterPro entry.list; unresolved IPR IDs dropped.",
              call. = FALSE)
      return(character())
    }
  }
  map <- .parse_interpro_entry_list(path)
  saveRDS(map, rds)
  map
}

#' Resolve GO tokens to names only (drop unresolved IDs)
#' @keywords internal
.e5_resolve_term_names <- function(tokens,
                                   cache_dir = NULL,
                                   obo_path = NULL,
                                   download = TRUE) {
  tokens <- unique(as.character(tokens))
  tokens <- tokens[nzchar(trimws(tokens))]
  if (!length(tokens)) {
    return(character())
  }
  map <- .go_term_name_map(cache_dir = cache_dir, obo_path = obo_path, download = download)
  out <- character()
  for (tok in tokens) {
    t <- trimws(tok)
    if (grepl("^GO[:_0-9]", t, ignore.case = TRUE)) {
      id <- .normalize_go_id(t)
      nm <- unname(map[id])
      if (!is.na(nm) && nzchar(nm)) {
        out <- c(out, nm)
      }
      # unresolved ID → drop (never pass ID as name)
    } else {
      out <- c(out, t)
    }
  }
  unique(out)
}

#' Classify a KEGG token into kind + bare id
#' @keywords internal
.e6_classify_token <- function(tok) {
  t <- trimws(as.character(tok)[1L])
  if (!nzchar(t)) {
    return(NULL)
  }
  # gene entry org:digits — skip
  if (grepl("^[a-z]{3}:[0-9]+$", t)) {
    return(NULL)
  }
  if (grepl("^K[0-9]{5}$", t, ignore.case = TRUE)) {
    return(list(id = toupper(t), kind = "ko"))
  }
  if (grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9n-]+$", t)) {
    return(list(id = t, kind = "ec"))
  }
  # map00500 / osa00500 / path:map00500
  t2 <- sub("^path:", "", t)
  if (grepl("^(map|[a-z]{3})[0-9]{5}$", t2)) {
    return(list(id = t2, kind = "pathway", org = if (grepl("^[a-z]{3}", t2)) {
      substr(t2, 1L, 3L)
    } else {
      NULL
    }))
  }
  # free name
  list(id = NULL, kind = NULL, name = t)
}

#' Resolve KEGG tokens to name+kind entries
#' @keywords internal
.e6_resolve_entries <- function(tokens,
                                cache_dir = NULL,
                                maps = NULL,
                                download = TRUE) {
  tokens <- unique(as.character(tokens))
  tokens <- tokens[nzchar(trimws(tokens))]
  if (!length(tokens)) {
    return(list())
  }
  if (is.null(maps)) {
    maps <- list(
      pathway = .kegg_list_name_map("pathway", cache_dir = cache_dir, download = download),
      ko = .kegg_list_name_map("ko", cache_dir = cache_dir, download = download),
      enzyme = .kegg_list_name_map("enzyme", cache_dir = cache_dir, download = download)
    )
  }
  org_maps <- list()
  entries <- list()
  seen <- character()
  for (tok in tokens) {
    cl <- .e6_classify_token(tok)
    if (is.null(cl)) next
    if (!is.null(cl$name)) {
      key <- paste0("name:", cl$name)
      if (key %in% seen) next
      seen <- c(seen, key)
      entries[[length(entries) + 1L]] <- list(name = cl$name, kind = NULL)
      next
    }
    kind <- cl$kind
    id <- cl$id
    mp <- if (identical(kind, "pathway")) {
      if (!is.null(cl$org) && nzchar(cl$org)) {
        if (is.null(org_maps[[cl$org]])) {
          org_maps[[cl$org]] <- .kegg_list_name_map(
            "pathway", org = cl$org, cache_dir = cache_dir, download = download
          )
        }
        org_maps[[cl$org]]
      } else {
        maps$pathway
      }
    } else if (identical(kind, "ko")) {
      maps$ko
    } else {
      maps$enzyme
    }
    nm <- unname(mp[id])
    if (is.na(nm) || !nzchar(nm)) {
      # try mapXXXX from osaXXXX
      if (identical(kind, "pathway") && grepl("^[a-z]{3}", id)) {
        mid <- paste0("map", substr(id, 4L, nchar(id)))
        nm <- unname(maps$pathway[mid])
      }
    }
    if (is.na(nm) || !nzchar(nm)) next
    key <- paste(kind, nm, sep = ":")
    if (key %in% seen) next
    seen <- c(seen, key)
    entries[[length(entries) + 1L]] <- list(name = nm, kind = kind)
  }
  entries
}

#' Resolve domain tokens to name+source
#' @keywords internal
.e7_resolve_domains <- function(tokens,
                                cache_dir = NULL,
                                pfam_path = NULL,
                                interpro_path = NULL,
                                download = TRUE) {
  tokens <- unique(as.character(tokens))
  tokens <- tokens[nzchar(trimws(tokens))]
  if (!length(tokens)) {
    return(list())
  }
  pfam <- .pfam_name_map(cache_dir = cache_dir, tsv_path = pfam_path, download = download)
  ipr <- .interpro_name_map(
    cache_dir = cache_dir, list_path = interpro_path, download = download
  )
  out <- list()
  seen <- character()
  for (tok in tokens) {
    t <- trimws(tok)
    source <- "other"
    nm <- NULL
    if (grepl("^PF[0-9]+", t, ignore.case = TRUE)) {
      acc <- toupper(sub("\\..*$", "", t))
      acc <- sub("^(PF[0-9]+).*", "\\1", acc)
      # pad PF00001 style
      m <- regexec("^PF([0-9]+)$", acc)
      h <- regmatches(acc, m)[[1L]]
      if (length(h) >= 2L) {
        acc <- paste0("PF", sprintf("%05d", as.integer(h[2L])))
      }
      nm <- unname(pfam[acc])
      source <- "pfam"
      if (is.na(nm) || !nzchar(nm)) next
    } else if (grepl("^IPR[0-9]+", t, ignore.case = TRUE)) {
      acc <- toupper(sub("\\..*$", "", t))
      m <- regexec("^IPR([0-9]+)$", acc)
      h <- regmatches(acc, m)[[1L]]
      if (length(h) >= 2L) {
        acc <- paste0("IPR", sprintf("%06d", as.integer(h[2L])))
      }
      nm <- unname(ipr[acc])
      source <- "interpro"
      if (is.na(nm) || !nzchar(nm)) next
    } else {
      # name-only
      nm <- t
      source <- "other"
    }
    key <- paste(source, nm, sep = ":")
    if (key %in% seen) next
    seen <- c(seen, key)
    out[[length(out) + 1L]] <- list(name = nm, source = source)
  }
  out
}
