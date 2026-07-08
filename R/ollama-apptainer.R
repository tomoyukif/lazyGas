################################################################################
#' Default directory for lazyGas Ollama models and wrapper scripts
#'
#' Uses \code{LAZYGAS_OLLAMA_HOME} when set; otherwise
#' \code{tools::R_user_dir("lazyGas", "cache")/ollama}.
#'
#' @return Character scalar (directory path).
#' @export
lazyGasOllamaHome <- function() {
  .lazygas_ollama_home()
}

#' @keywords internal
.lazygas_ollama_home <- function() {
  env <- Sys.getenv("LAZYGAS_OLLAMA_HOME", unset = "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  if (requireNamespace("tools", quietly = TRUE)) {
    return(file.path(tools::R_user_dir("lazyGas", which = "cache"), "ollama"))
  }
  file.path(path.expand("~/.cache/lazygas"), "ollama")
}

#' @keywords internal
.apptainer_bin <- function() {
  env <- Sys.getenv("LAZYGAS_APPTAINER_BIN", unset = "")
  if (nzchar(env) && file.exists(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  bin <- Sys.which("apptainer")
  if (nzchar(bin)) {
    return(bin)
  }
  Sys.which("singularity")
}

#' @keywords internal
.lazygas_ollama_sif_default <- function() {
  env <- Sys.getenv("LAZYGAS_OLLAMA_SIF", unset = "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  file.path(.lazygas_ollama_home(), "ollama.sif")
}

#' @keywords internal
.lazygas_ollama_wrapper_default <- function() {
  env <- Sys.getenv("LAZYGAS_OLLAMA_CLI", unset = "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  opt <- getOption("lazygas.ollama.cli", default = "")
  if (nzchar(opt)) {
    return(normalizePath(opt, winslash = "/", mustWork = FALSE))
  }
  file.path(.lazygas_ollama_home(), "bin", "lazygas-ollama")
}

#' @keywords internal
.lazygas_ollama_host_from_url <- function(base_url) {
  url <- sub("^https?://", "", .llm_base_url(base_url))
  if (grepl(":", url)) {
    return(url)
  }
  paste0(url, ":11434")
}

#' @keywords internal
.write_lazygas_ollama_wrapper <- function(wrapper_path,
                                            apptainer_bin,
                                            sif_path,
                                            home_dir,
                                            base_url) {
  template <- system.file("apptainer", "lazygas-ollama.in", package = "lazyGas")
  if (!nzchar(template) || !file.exists(template)) {
    dev <- file.path(getwd(), "inst", "apptainer", "lazygas-ollama.in")
    if (file.exists(dev)) {
      template <- dev
    } else {
      stop("Wrapper template inst/apptainer/lazygas-ollama.in not found.", call. = FALSE)
    }
  }
  txt <- readLines(template, warn = FALSE)
  txt <- gsub("__APPTAINER_BIN__", apptainer_bin, txt, fixed = TRUE)
  txt <- gsub("__SIF_PATH__", sif_path, txt, fixed = TRUE)
  txt <- gsub("__OLLAMA_HOME__", home_dir, txt, fixed = TRUE)
  txt <- gsub(
    "__OLLAMA_HOST__",
    .lazygas_ollama_host_from_url(base_url),
    txt,
    fixed = TRUE
  )
  writeLines(txt, wrapper_path)
  Sys.chmod(wrapper_path, mode = "0755")
  invisible(wrapper_path)
}

#' Build the official lazyGas Ollama Apptainer image (\code{.sif})
#'
#' Pulls \code{docker://ollama/ollama:latest} by default. Requires
#' \code{apptainer} or \code{singularity} on \code{PATH}.
#'
#' @param dest Output \code{.sif} path (default: under [lazyGasOllamaHome()]).
#' @param method \code{"pull"} (recommended) or \code{"build"} from
#'   \code{inst/apptainer/ollama.def}.
#' @param apptainer_bin Path to apptainer/singularity binary.
#'
#' @return Invisibly, the path to the created \code{.sif} file.
#' @export
buildLazyGasOllamaSif <- function(dest = NULL,
                                method = c("pull", "build"),
                                apptainer_bin = NULL) {
  method <- match.arg(method)
  apptainer_bin <- apptainer_bin %||% .apptainer_bin()
  if (!nzchar(apptainer_bin)) {
    stop(
      "Apptainer/Singularity not found. Install apptainer or set LAZYGAS_APPTAINER_BIN.",
      call. = FALSE
    )
  }

  dest <- dest %||% .lazygas_ollama_sif_default()
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  if (method == "pull") {
    message("Pulling docker://ollama/ollama:latest -> ", dest)
    exit_code <- system2(
      apptainer_bin,
      args = c("pull", dest, "docker://ollama/ollama:latest")
    )
  } else {
    def <- system.file("apptainer", "ollama.def", package = "lazyGas")
    if (!nzchar(def) || !file.exists(def)) {
      def <- file.path(getwd(), "inst", "apptainer", "ollama.def")
    }
    if (!file.exists(def)) {
      stop("Apptainer definition file not found.", call. = FALSE)
    }
    message("Building ", dest, " from ", def)
    exit_code <- system2(apptainer_bin, args = c("build", dest, def))
  }

  if (!is.null(exit_code) && exit_code != 0L) {
    stop("Failed to create Ollama SIF (exit code ", exit_code, ").", call. = FALSE)
  }
  if (!file.exists(dest)) {
    stop("Ollama SIF was not created: ", dest, call. = FALSE)
  }

  message("Ollama SIF ready: ", dest)
  invisible(normalizePath(dest, winslash = "/", mustWork = FALSE))
}

#' Configure lazyGas to use Ollama inside Apptainer (official setup)
#'
#' Installs a wrapper script (\code{lazygas-ollama}) and sets environment variables
#' so [startOllamaServer()], [pullOllamaModel()], and the Shiny explorer use the
#' containerised Ollama CLI instead of a host installation.
#'
#' @param sif Path to the Ollama \code{.sif} image. Built via [buildLazyGasOllamaSif()]
#'   when missing and \code{build_sif_if_missing = TRUE}.
#' @param home_dir Host directory bound to \code{/root/.ollama} (models/cache).
#' @param apptainer_bin Path to \code{apptainer} or \code{singularity}.
#' @param base_url HTTP URL for the Ollama API (default \code{http://127.0.0.1:11434}).
#' @param model Default model name (\code{LAZYGAS_LLM_MODEL}).
#' @param build_sif_if_missing Run [buildLazyGasOllamaSif()] when \code{sif} is absent.
#'
#' @return Invisibly, a list with \code{wrapper}, \code{sif}, \code{home_dir},
#'   and \code{base_url}.
#' @export
#'
#' @seealso [buildLazyGasOllamaSif()], [startLocalLLM()], [describeOllamaStatus()]
configureLazyGasOllama <- function(sif = NULL,
                                    home_dir = NULL,
                                    apptainer_bin = NULL,
                                    base_url = "http://127.0.0.1:11434",
                                    model = "llama3.2:3b",
                                    build_sif_if_missing = TRUE) {
  apptainer_bin <- apptainer_bin %||% .apptainer_bin()
  if (!nzchar(apptainer_bin)) {
    stop(
      "Apptainer/Singularity not found. Install apptainer (https://apptainer.org).",
      call. = FALSE
    )
  }

  home_dir <- home_dir %||% .lazygas_ollama_home()
  dir.create(home_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(home_dir, "bin"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(home_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

  sif <- sif %||% .lazygas_ollama_sif_default()
  if (!file.exists(sif)) {
    if (!isTRUE(build_sif_if_missing)) {
      stop(
        "Ollama SIF not found: ", sif,
        ". Run buildLazyGasOllamaSif() first.",
        call. = FALSE
      )
    }
    sif <- buildLazyGasOllamaSif(dest = sif, apptainer_bin = apptainer_bin)
  }

  wrapper <- file.path(home_dir, "bin", "lazygas-ollama")
  .write_lazygas_ollama_wrapper(
    wrapper_path = wrapper,
    apptainer_bin = apptainer_bin,
    sif_path = sif,
    home_dir = home_dir,
    base_url = base_url
  )

  Sys.setenv(
    LAZYGAS_APPTAINER_BIN = apptainer_bin,
    LAZYGAS_OLLAMA_SIF = sif,
    LAZYGAS_OLLAMA_HOME = home_dir,
    LAZYGAS_OLLAMA_CLI = wrapper,
    LAZYGAS_LLM_URL = base_url,
    LAZYGAS_LLM_MODEL = model
  )
  options(
    lazygas.ollama.cli = wrapper,
    lazygas.ollama.runtime = "apptainer"
  )

  message("lazyGas Ollama (Apptainer) configured.")
  message("  SIF: ", sif)
  message("  Home: ", home_dir)
  message("  CLI wrapper: ", wrapper)

  invisible(list(
    wrapper = wrapper,
    sif = sif,
    home_dir = home_dir,
    base_url = .llm_base_url(base_url),
    model = model
  ))
}

#' @keywords internal
.lazygas_ollama_ensure_configured <- function() {
  if (nzchar(.ollama_cli_path())) {
    return(invisible(TRUE))
  }
  if (!nzchar(.apptainer_bin())) {
    return(invisible(FALSE))
  }
  sif <- .lazygas_ollama_sif_default()
  if (file.exists(sif)) {
    configureLazyGasOllama(sif = sif, build_sif_if_missing = FALSE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

#' @keywords internal
.lazygas_ollama_using_apptainer <- function() {
  identical(getOption("lazygas.ollama.runtime", ""), "apptainer") ||
    nzchar(Sys.getenv("LAZYGAS_OLLAMA_SIF", unset = ""))
}
