################################################################################
#' AI-lazyGas configuration helpers (Phase 1)
#'
#' Resolution priority for settings: **function arguments > YAML config >
#' environment variables > built-in defaults**.
#'
#' @name ai-config
#' @keywords internal
NULL

`%||%` <- function(x, y) if (is.null(x)) y else x

.PHASE1_DEFAULT_WEIGHTS <- c(
  annotation = 0.30,
  snpeff = 0.35,
  gwas = 0.25,
  expression = 0.10,
  finemap = 0.50
)

.PHASE1_DEFAULT_SOURCES <- c("annotation", "snpeff", "gwas", "expression")

.EVIDENCE_SOURCE_CHOICES <- c(
  "annotation", "snpeff", "gwas", "expression", "literature", "ortholog",
  "finemap"
)

#' Load AI-lazyGas YAML config
#'
#' @param which One of \code{"models"}, \code{"data"}, or \code{"both"}.
#' @return Named list (\code{models} / \code{data} elements when \code{both}).
#' @export
loadAiLazyGasConfig <- function(which = c("both", "models", "data")) {
  which <- match.arg(which)
  models <- if (which %in% c("both", "models")) {
    .ai_config_read_yaml("ai-lazygas-models.yaml")
  } else {
    NULL
  }
  data <- if (which %in% c("both", "data")) {
    .ai_config_read_yaml("ai-lazygas-data.yaml")
  } else {
    NULL
  }
  if (identical(which, "models")) {
    return(models)
  }
  if (identical(which, "data")) {
    return(data)
  }
  list(models = models, data = data)
}

#' Phase 1 data path resolver
#'
#' @param key Key under \code{phase1} in \code{ai-lazygas-data.yaml}
#'   (e.g. \code{"snpeff_gds"}, \code{"gff"}, \code{"ann"}).
#' @param env_var Optional environment variable override.
#' @param default Fallback path when unset.
#' @param must_exist Stop if the resolved path does not exist.
#' @return Character path (may be empty when unset and \code{must_exist = FALSE}).
#' @export
phase1DataPath <- function(key,
                           env_var = NULL,
                           default = "",
                           must_exist = FALSE) {
  key <- as.character(key)[[1L]]
  if (!is.null(env_var) && nzchar(env_var)) {
    env <- Sys.getenv(env_var, unset = "")
    if (nzchar(env)) {
      path <- normalizePath(env, winslash = "/", mustWork = FALSE)
      if (isTRUE(must_exist) && !file.exists(path)) {
        stop("Path from ", env_var, " not found: ", path, call. = FALSE)
      }
      return(path)
    }
  }
  cfg <- .ai_config_phase1_data()
  val <- cfg[[key]]
  if (!is.null(val) && nzchar(as.character(val)[[1L]])) {
    path <- normalizePath(as.character(val)[[1L]], winslash = "/", mustWork = FALSE)
  } else if (nzchar(default)) {
    path <- normalizePath(default, winslash = "/", mustWork = FALSE)
  } else {
    path <- ""
  }
  if (isTRUE(must_exist) && (!nzchar(path) || !file.exists(path))) {
    stop("Phase 1 path '", key, "' not found: ", path, call. = FALSE)
  }
  path
}

#' Directory for AI reports (\code{lazygas/ai_report/})
#'
#' @param object Optional \code{LazyGas} object. Reports are stored under
#'   \code{<dirname(store)>/lazygas/ai_report/}.
#' @param out_dir Explicit directory (overrides \code{object}).
#' @return Character path.
#' @export
aiReportDir <- function(object = NULL, out_dir = NULL) {
  if (!is.null(out_dir) && nzchar(out_dir)) {
    return(normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  }
  if (!is.null(object) && inherits(object, "LazyGas")) {
    root <- dirname(.store_path(object))
    return(file.path(root, "lazygas", "ai_report"))
  }
  normalizePath(
    file.path("lazygas", "ai_report"),
    winslash = "/",
    mustWork = FALSE
  )
}

#' Default Phase 1 ranking weights
#' @return Named numeric vector.
#' @export
phase1DefaultWeights <- function() {
  cfg <- .ai_config_models()
  w <- cfg$weights
  if (is.list(w) && length(w)) {
    out <- unlist(w)
    out <- as.numeric(out)
    names(out) <- names(w)
    keep <- intersect(
      c("annotation", "snpeff", "gwas", "expression", "finemap",
        "literature", "ortholog"),
      names(out)
    )
    if (length(keep)) {
      return(out[keep])
    }
  }
  .PHASE1_DEFAULT_WEIGHTS
}

#' Default Phase 1 evidence sources
#' @return Character vector.
#' @export
phase1DefaultSources <- function() {
  .PHASE1_DEFAULT_SOURCES
}

# ---- internals -------------------------------------------------------------

.ai_config_find <- function(filename) {
  pkg <- tryCatch(
    system.file("config", filename, package = "lazyGas"),
    error = function(e) ""
  )
  if (nzchar(pkg) && file.exists(pkg)) {
    return(pkg)
  }
  # Development tree / source install
  candidates <- c(
    file.path(getwd(), "inst", "config", filename),
    file.path(dirname(getwd()), "inst", "config", filename)
  )
  for (p in candidates) {
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/"))
    }
  }
  ""
}

.ai_config_read_yaml <- function(filename) {
  path <- .ai_config_find(filename)
  if (!nzchar(path) || !file.exists(path)) {
    return(list())
  }
  if (requireNamespace("yaml", quietly = TRUE)) {
    return(yaml::read_yaml(path))
  }
  .ai_config_parse_simple_yaml(path)
}

.ai_config_parse_simple_yaml <- function(path) {
  # Minimal indented-key parser for Phase 1 configs (no nested lists beyond maps).
  lines <- readLines(path, warn = FALSE)
  lines <- sub("#.*$", "", lines)
  root <- list()
  stack <- list(list(indent = -1L, target = "root"))
  env <- new.env(parent = emptyenv())
  env$root <- root

  set_at <- function(path_keys, value) {
    cur <- env$root
    if (!length(path_keys)) {
      return(invisible(NULL))
    }
    assign_recursive <- function(node, keys, val) {
      k <- keys[[1L]]
      if (length(keys) == 1L) {
        node[[k]] <- val
        return(node)
      }
      if (is.null(node[[k]]) || !is.list(node[[k]])) {
        node[[k]] <- list()
      }
      node[[k]] <- assign_recursive(node[[k]], keys[-1L], val)
      node
    }
    env$root <- assign_recursive(env$root, path_keys, value)
  }

  key_stack <- character()
  indent_stack <- -1L

  for (line in lines) {
    if (!nzchar(trimws(line))) {
      next
    }
    indent <- nchar(sub("^(\\s*).*$", "\\1", line))
    content <- trimws(line)
    if (!grepl(":", content, fixed = TRUE)) {
      next
    }
    key <- trimws(sub(":.*$", "", content))
    val <- trimws(sub("^[^:]+:\\s*", "", content))
    val <- gsub("^['\"]|['\"]$", "", val)

    while (length(indent_stack) && indent <= indent_stack[[length(indent_stack)]]) {
      indent_stack <- indent_stack[-length(indent_stack)]
      key_stack <- key_stack[-length(key_stack)]
    }

    if (!nzchar(val)) {
      key_stack <- c(key_stack, key)
      indent_stack <- c(indent_stack, indent)
      set_at(key_stack, list())
    } else {
      num <- suppressWarnings(as.numeric(val))
      parsed <- if (!is.na(num) && grepl("^-?[0-9.]+$", val)) {
        num
      } else if (identical(tolower(val), "true")) {
        TRUE
      } else if (identical(tolower(val), "false")) {
        FALSE
      } else {
        val
      }
      set_at(c(key_stack, key), parsed)
    }
  }
  env$root
}

.ai_config_models <- function() {
  .ai_config_read_yaml("ai-lazygas-models.yaml")
}

.ai_config_phase1_data <- function() {
  cfg <- .ai_config_read_yaml("ai-lazygas-data.yaml")
  if (is.null(cfg$phase1)) {
    return(list())
  }
  cfg$phase1
}

.ai_config_llm_settings <- function(model = NULL, base_url = NULL, timeout = NULL) {
  cfg <- .ai_config_models()
  def <- cfg$default %||% list()
  agents <- cfg$agents %||% list()
  main <- agents$main %||% list()

  resolved_model <- if (!is.null(model) && nzchar(model)) {
    model
  } else if (!is.null(main$model) && nzchar(main$model)) {
    as.character(main$model)
  } else {
    Sys.getenv("LAZYGAS_LLM_MODEL", unset = "gemma4-31b-64k:latest")
  }

  resolved_url <- if (!is.null(base_url) && nzchar(base_url)) {
    base_url
  } else if (!is.null(def$base_url) && nzchar(def$base_url)) {
    as.character(def$base_url)
  } else {
    Sys.getenv("LAZYGAS_LLM_URL", unset = "http://127.0.0.1:11435")
  }

  resolved_timeout <- if (!is.null(timeout) && is.finite(as.numeric(timeout)[1L])) {
    as.numeric(timeout)[1L]
  } else if (!is.null(def$timeout) && is.finite(as.numeric(def$timeout)[1L])) {
    as.numeric(def$timeout)[1L]
  } else {
    600
  }

  list(
    model = resolved_model,
    base_url = sub("/+$", "", resolved_url),
    timeout = resolved_timeout,
    ranking = agents$ranking %||% list()
  )
}

.ai_config_ranking_flags <- function() {
  cfg <- .ai_config_models()
  ranking <- cfg$agents$ranking %||% list()
  list(
    use_keyword = isTRUE(ranking$use_keyword %||% TRUE),
    use_semantic = FALSE,
    use_llm_relevance = isTRUE(ranking$use_llm_relevance %||% TRUE),
    llm_model = as.character(
      ranking$llm_model %||%
        (cfg$agents$main$model %||% "gemma4-31b-64k:latest")
    )
  )
}
