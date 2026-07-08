################################################################################
#' Locate the Ollama command-line binary
#'
#' @return Character scalar (empty if not found).
#' @keywords internal
.ollama_cli_path <- function() {
  wrapper <- .lazygas_ollama_wrapper_default()
  if (nzchar(wrapper) && file.exists(wrapper)) {
    return(normalizePath(wrapper, winslash = "/", mustWork = FALSE))
  }
  env_cli <- Sys.getenv("LAZYGAS_OLLAMA_CLI", unset = "")
  if (nzchar(env_cli) && file.exists(env_cli)) {
    return(normalizePath(env_cli, winslash = "/", mustWork = FALSE))
  }
  Sys.which("ollama")
}

#' @keywords internal
.llm_ollama_state_env <- function() {
  if (!exists(".lazygas_ollama_state", envir = .GlobalEnv, inherits = FALSE)) {
    assign(".lazygas_ollama_state", new.env(parent = emptyenv()), envir = .GlobalEnv)
  }
  get(".lazygas_ollama_state", envir = .GlobalEnv)
}

#' @keywords internal
.llm_ollama_state_key <- function(base_url) {
  sub("/+$", "", .llm_base_url(base_url))
}

#' @keywords internal
.llm_ollama_get_state <- function(base_url) {
  env <- .llm_ollama_state_env()
  key <- .llm_ollama_state_key(base_url)
  if (!exists(key, envir = env, inherits = FALSE)) {
    return(NULL)
  }
  get(key, envir = env)
}

#' @keywords internal
.llm_ollama_set_state <- function(base_url, state) {
  env <- .llm_ollama_state_env()
  key <- .llm_ollama_state_key(base_url)
  assign(key, state, envir = env)
}

#' @keywords internal
.llm_ollama_clear_state <- function(base_url) {
  env <- .llm_ollama_state_env()
  key <- .llm_ollama_state_key(base_url)
  if (exists(key, envir = env, inherits = FALSE)) {
    rm(list = key, envir = env)
  }
}

#' @keywords internal
.llm_wait_for_server <- function(base_url, wait_seconds = 30, interval = 0.5) {
  deadline <- Sys.time() + wait_seconds
  while (Sys.time() < deadline) {
    if (llmHealthCheck(base_url = base_url, timeout = 2)) {
      return(TRUE)
    }
    Sys.sleep(interval)
  }
  FALSE
}

#' List models available on a local Ollama server
#'
#' @param base_url Ollama base URL.
#' @return Character vector of model names.
#' @export
listOllamaModels <- function(base_url = NULL) {
  base_url <- .llm_base_url(base_url)
  if (!llmHealthCheck(base_url = base_url, timeout = 5)) {
    stop("Ollama server is not reachable at ", base_url, call. = FALSE)
  }
  resp <- .llm_http_get(path = "/api/tags", base_url = base_url, timeout = 10)
  if (is.null(resp) || is.null(resp$models) || !length(resp$models)) {
    return(character())
  }
  vapply(resp$models, function(m) as.character(m$name), character(1L))
}

#' @keywords internal
.ollama_model_available <- function(model, base_url) {
  installed <- listOllamaModels(base_url = base_url)
  if (!length(installed)) {
    return(FALSE)
  }
  model <- sub(":latest$", "", model)
  any(vapply(installed, function(m) {
    m <- sub(":latest$", "", m)
    identical(m, model) || startsWith(m, paste0(model, ":"))
  }, logical(1L)))
}

#' Start the Ollama server as a background process
#'
#' Spawns \code{ollama serve} in a separate process (unless a server is already
#' reachable at \code{base_url}). The process PID is recorded so that
#' [stopOllamaServer()] can stop only servers started from this R session.
#'
#' Requires the \code{ollama} CLI on \code{PATH}, or an Apptainer wrapper from
#' [configureLazyGasOllama()] (official lazyGas setup).
#'
#' @param base_url Ollama base URL.
#' @param wait_seconds Seconds to wait for the HTTP API to become ready.
#' @param log_file Optional log file path (temporary file by default).
#'
#' @return Invisibly, a list with \code{pid}, \code{base_url}, \code{log_file},
#'   and \code{already_running}.
#' @export
startOllamaServer <- function(base_url = NULL,
                              wait_seconds = 45,
                              log_file = NULL) {
  base_url <- .llm_base_url(base_url)

  .lazygas_ollama_ensure_configured()

  if (llmHealthCheck(base_url = base_url, timeout = 2)) {
    message("Ollama is already running at ", base_url)
    return(invisible(list(
      pid = NA_integer_,
      base_url = base_url,
      log_file = NA_character_,
      already_running = TRUE
    )))
  }

  ollama_bin <- .ollama_cli_path()
  if (!nzchar(ollama_bin)) {
    stop(
      "Ollama CLI not found. Official lazyGas setup: install Apptainer and run ",
      "configureLazyGasOllama() (see vignette). Alternatively install Ollama from ",
      "https://ollama.com and ensure 'ollama' is on PATH.",
      call. = FALSE
    )
  }

  if (is.null(log_file) || !nzchar(log_file)) {
    if (.lazygas_ollama_using_apptainer()) {
      log_file <- file.path(.lazygas_ollama_home(), "logs", "ollama_serve.log")
    } else {
      log_file <- tempfile(pattern = "lazygas_ollama_", fileext = ".log")
    }
  }
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

  state <- .llm_ollama_spawn(ollama_bin = ollama_bin, log_file = log_file)
  state$base_url <- base_url
  state$started_at <- Sys.time()
  .llm_ollama_set_state(base_url, state)

  if (!.llm_wait_for_server(base_url = base_url, wait_seconds = wait_seconds)) {
    stop(
      "Ollama did not become ready within ", wait_seconds,
      " seconds. Log: ", log_file,
      call. = FALSE
    )
  }

  message("Ollama started (pid ", state$pid, "). Log: ", log_file)
  invisible(list(
    pid = state$pid,
    base_url = base_url,
    log_file = log_file,
    already_running = FALSE
  ))
}

#' @keywords internal
.llm_ollama_spawn <- function(ollama_bin, log_file) {
  if (requireNamespace("processx", quietly = TRUE)) {
    proc <- processx::process$new(
      command = ollama_bin,
      args = "serve",
      stdout = log_file,
      stderr = log_file,
      cleanup = FALSE
    )
    return(list(
      pid = proc$get_pid(),
      process = proc,
      log_file = log_file,
      method = "processx"
    ))
  }

  if (.Platform$OS.type != "unix") {
    stop(
      "Install the suggested package 'processx' to start Ollama from R on this OS.",
      call. = FALSE
    )
  }

  cmd <- sprintf(
    "nohup %s serve >> %s 2>&1 & echo $!",
    shQuote(ollama_bin),
    shQuote(log_file)
  )
  pid_raw <- system(cmd, intern = TRUE)
  pid <- suppressWarnings(as.integer(pid_raw[[1L]]))
  if (!is.finite(pid)) {
    stop("Failed to start Ollama background process.", call. = FALSE)
  }
  list(
    pid = pid,
    process = NULL,
    log_file = log_file,
    method = "nohup"
  )
}

#' Stop an Ollama server started by [startOllamaServer()]
#'
#' Only stops processes recorded in the current R session. Does not stop a
#' system-wide Ollama service started outside lazyGas.
#'
#' @param base_url Ollama base URL.
#' @param only_started_by_lazygas If \code{TRUE}, only stop servers started via
#'   [startOllamaServer()].
#'
#' @return Logical scalar indicating whether a process was stopped.
#' @export
stopOllamaServer <- function(base_url = NULL, only_started_by_lazygas = TRUE) {
  base_url <- .llm_base_url(base_url)
  state <- .llm_ollama_get_state(base_url)

  if (only_started_by_lazygas && is.null(state)) {
    message("No Ollama process started by lazyGas for ", base_url)
    return(invisible(FALSE))
  }

  pid <- state$pid
  if (!is.null(state$process) && requireNamespace("processx", quietly = TRUE)) {
    try(state$process$kill(), silent = TRUE)
  } else if (is.finite(pid)) {
    try(tools::pskill(pid), silent = TRUE)
  }

  .llm_ollama_clear_state(base_url)
  Sys.sleep(0.5)
  still_up <- llmHealthCheck(base_url = base_url, timeout = 2)
  if (still_up && only_started_by_lazygas) {
    message(
      "Ollama API is still reachable (likely a system service). ",
      "lazyGas stopped only the process it spawned."
    )
    return(invisible(FALSE))
  }

  message("Ollama stopped for ", base_url)
  invisible(TRUE)
}

#' Pull an Ollama model
#'
#' @param model Model name (default from \code{LAZYGAS_LLM_MODEL}).
#' @param base_url Ollama base URL (must already be running).
#' @return Invisibly \code{TRUE}.
#' @export
pullOllamaModel <- function(model = NULL, base_url = NULL) {
  model <- .llm_model(model)
  base_url <- .llm_base_url(base_url)
  if (!llmHealthCheck(base_url = base_url, timeout = 5)) {
    stop(
      "Ollama server is not running. Call startOllamaServer() or startLocalLLM().",
      call. = FALSE
    )
  }
  ollama_bin <- .ollama_cli_path()
  if (!nzchar(ollama_bin)) {
    stop(
      "Ollama CLI not found. Run configureLazyGasOllama() or install Ollama on PATH.",
      call. = FALSE
    )
  }
  message("Pulling Ollama model '", model, "' (may take several minutes)...")
  exit_code <- system2(ollama_bin, args = c("pull", model))
  if (!is.null(exit_code) && exit_code != 0L) {
    stop("ollama pull failed (exit code ", exit_code, ").", call. = FALSE)
  }
  invisible(TRUE)
}

#' Start Ollama and optionally pull a model (convenience wrapper)
#'
#' @param model Model name.
#' @param base_url Ollama base URL.
#' @param pull_model If \code{TRUE}, run \code{ollama pull} when the model is
#'   missing.
#' @param start_server If \code{TRUE}, call [startOllamaServer()] when needed.
#' @param wait_seconds Seconds to wait for the server API.
#'
#' @return Invisibly, a list with \code{base_url}, \code{model},
#'   \code{running}, and \code{model_available}.
#' @export
startLocalLLM <- function(model = NULL,
                          base_url = NULL,
                          pull_model = TRUE,
                          start_server = TRUE,
                          wait_seconds = 45) {
  model <- .llm_model(model)
  base_url <- .llm_base_url(base_url)

  if (!llmHealthCheck(base_url = base_url, timeout = 2)) {
    if (!isTRUE(start_server)) {
      stop("Ollama is not running and start_server = FALSE.", call. = FALSE)
    }
    startOllamaServer(base_url = base_url, wait_seconds = wait_seconds)
  }

  model_available <- FALSE
  if (isTRUE(pull_model)) {
    if (.ollama_model_available(model = model, base_url = base_url)) {
      model_available <- TRUE
    } else {
      pullOllamaModel(model = model, base_url = base_url)
      model_available <- TRUE
    }
  } else {
    model_available <- .ollama_model_available(model = model, base_url = base_url)
  }

  invisible(list(
    base_url = base_url,
    model = model,
    running = llmHealthCheck(base_url = base_url, timeout = 5),
    model_available = model_available
  ))
}

#' Summarize local Ollama status for UI or diagnostics
#'
#' @param base_url Ollama base URL.
#' @param model Optional model name to check availability.
#' @return Character string.
#' @export
describeOllamaStatus <- function(base_url = NULL, model = NULL) {
  base_url <- .llm_base_url(base_url)
  cli <- .ollama_cli_path()
  if (!nzchar(cli)) {
    if (nzchar(.apptainer_bin())) {
      return(paste0(
        "Ollama not configured (Apptainer: ", .apptainer_bin(), "). ",
        "Run configureLazyGasOllama()."
      ))
    }
    return("Ollama CLI not found. Run configureLazyGasOllama() or install Ollama.")
  }
  prefix <- ""
  if (.lazygas_ollama_using_apptainer()) {
    sif <- Sys.getenv("LAZYGAS_OLLAMA_SIF", unset = .lazygas_ollama_sif_default())
    prefix <- paste0("Runtime: Apptainer (", basename(sif), ")\n")
  }
  if (!llmHealthCheck(base_url = base_url, timeout = 3)) {
    return(paste0(prefix, "Not running (", base_url, ")"))
  }
  models <- tryCatch(
    listOllamaModels(base_url = base_url),
    error = function(e) character()
  )
  line <- paste0(prefix, "Running (", base_url, ")")
  if (length(models)) {
    line <- paste0(line, "\nModels: ", paste(models, collapse = ", "))
  }
  if (!is.null(model) && nzchar(model)) {
    model <- .llm_model(model)
    avail <- .ollama_model_available(model = model, base_url = base_url)
    line <- paste0(
      line,
      "\nSelected model '", model, "': ",
      if (avail) "available" else "not installed (use pullOllamaModel())"
    )
  }
  line
}
