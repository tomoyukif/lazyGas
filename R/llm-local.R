################################################################################
#' Check connectivity to a local Ollama server
#'
#' @param base_url Ollama base URL.
#' @param timeout Connection timeout in seconds.
#' @return Logical scalar.
#' @export
llmHealthCheck <- function(base_url = NULL, timeout = 5) {
  base_url <- .llm_base_url(base_url)
  resp <- .llm_http_get(
    path = "/api/tags",
    base_url = base_url,
    timeout = timeout
  )
  isTRUE(!is.null(resp))
}

#' Send a chat completion request to a local Ollama server
#'
#' @param messages List of lists with \code{role} and \code{content}.
#' @param model Model name (default from \code{LAZYGAS_LLM_MODEL}).
#' @param base_url Ollama base URL.
#' @param json_mode If \code{TRUE}, request JSON-formatted output.
#' @param timeout Request timeout in seconds.
#'
#' @return Assistant message content as a character string.
#' @export
llmChat <- function(messages,
                    model = NULL,
                    base_url = NULL,
                    json_mode = FALSE,
                    timeout = 120) {
  if (!is.list(messages) || !length(messages)) {
    stop("'messages' must be a non-empty list.", call. = FALSE)
  }
  model <- .llm_model(model)
  base_url <- .llm_base_url(base_url)

  body <- list(
    model = model,
    messages = messages,
    stream = FALSE
  )
  if (isTRUE(json_mode)) {
    body$format <- "json"
  }

  resp <- .llm_http_post(
    path = "/api/chat",
    body = body,
    base_url = base_url,
    timeout = timeout
  )
  if (is.null(resp) || is.null(resp$message$content)) {
    stop("Empty response from local LLM.", call. = FALSE)
  }
  as.character(resp$message$content)
}

#' One constrained-output scrutiny pass: prior answer + rules → LLM
#'
#' Used when a restricted generation step produced invalid output. Sends the
#' original constraint rules, the previous model output, and problem notes so
#' the LLM can correct once. Does not loop.
#'
#' @param llm List with \code{model}, \code{base_url}, optional \code{timeout}.
#' @param rules Constraint / schema text (system-level rules to obey).
#' @param prior_output Previous model output (JSON or prose).
#' @param problems Character vector of validation failures.
#' @param extra_user Optional extra user context (e.g. allowed numbers).
#' @return Corrected content string, or an \code{error} condition on failure.
#' @keywords internal
.llm_scrutiny_pass <- function(llm,
                               rules,
                               prior_output,
                               problems = character(),
                               extra_user = NULL,
                               timeout = NULL) {
  if (is.null(llm) || is.null(llm$model)) {
    return(simpleError("llm_scrutiny: missing llm settings"))
  }
  sys <- paste(
    "You audit and correct a previous model answer that violated constraints.",
    "Obey the constraint rules exactly. Fix every listed problem.",
    "Keep anything that already complies. Return only the corrected answer",
    "(same format as required by the rules; JSON when rules require JSON).",
    "Do not add commentary outside that format."
  )
  user <- paste0(
    "Constraint rules:\n",
    as.character(rules)[1L],
    "\n\nProblems found:\n",
    if (length(problems)) {
      paste(paste0("- ", problems), collapse = "\n")
    } else {
      "- (unspecified validation failure)"
    },
    "\n\nPrevious model output to correct:\n",
    as.character(prior_output)[1L]
  )
  if (!is.null(extra_user) && nzchar(as.character(extra_user)[1L])) {
    user <- paste0(user, "\n\nAdditional context:\n", as.character(extra_user)[1L])
  }
  to <- timeout %||% llm$timeout %||% 120
  tryCatch(
    llmChat(
      messages = list(
        list(role = "system", content = sys),
        list(role = "user", content = user)
      ),
      model = llm$model,
      base_url = llm$base_url,
      json_mode = TRUE,
      timeout = to
    ),
    error = function(e) e
  )
}

#' Explain ranked phenotype candidates using a local LLM
#'
#' Generates a Markdown summary of the top-ranked genes using only evidence
#' collected during ranking. When \code{use_llm = FALSE}, a template-based
#' report is returned.
#'
#' For Phase 1, prefer the public entry point [llm_report()], which ranks
#' candidates, verifies claims, and saves under \code{lazygas/ai_report/}.
#'
#' @param rank_result Output of [rankPhenotypeCandidates()].
#' @param query A \code{PhenotypeQuery} object (optional; inferred from
#'   \code{rank_result$query_id} when \code{object} is given).
#' @param object \code{LazyGas} object to load saved queries.
#' @param top_n Number of top genes to include in the explanation.
#' @param use_llm Use local LLM when available.
#' @param language Response language: \code{"ja"} or \code{"en"}.
#' @param model,base_url LLM settings passed to [llmChat()].
#' @param timeout Request timeout in seconds (default 600).
#'
#' @return A character string (Markdown).
#' @export
explainPhenotypeCandidates <- function(rank_result,
                                       query = NULL,
                                       object = NULL,
                                       top_n = 10L,
                                       use_llm = TRUE,
                                       language = c("ja", "en"),
                                       model = NULL,
                                       base_url = NULL,
                                       timeout = NULL) {
  language <- match.arg(language)
  if (!is.data.frame(rank_result) || nrow(rank_result) == 0L) {
    stop("'rank_result' must be a non-empty data.frame.", call. = FALSE)
  }

  meta <- attr(rank_result, "phenotypeRank")
  if (is.null(query) && !is.null(meta$query_id) && !is.null(object)) {
    query <- .store_read_phenotype_query(object, meta$query_id)
  }
  if (is.null(query)) {
    query <- list(
      trait_text = "phenotype query",
      query_id = if (!is.null(meta)) meta$query_id else "unknown"
    )
    class(query) <- c("PhenotypeQuery", "list")
  }

  llm <- .ai_config_llm_settings(model = model, base_url = base_url, timeout = timeout)
  top_n <- as.integer(top_n)[1L]
  sub <- rank_result[seq_len(min(top_n, nrow(rank_result))), , drop = FALSE]
  evidence_bundle <- .phenotype_evidence_bundle(rank_result = sub)

  if (!isTRUE(use_llm)) {
    return(.explain_phenotype_template(
      query = query,
      evidence_bundle = evidence_bundle,
      language = language
    ))
  }

  if (!llmHealthCheck(base_url = llm$base_url, timeout = 5)) {
    stop(
      "Local LLM is unreachable at ", llm$base_url,
      ". Start the SSH tunnel / Ollama, or set use_llm = FALSE for the template path.",
      call. = FALSE
    )
  }

  lang_note <- if (language == "ja") {
    "回答は日本語で書いてください。"
  } else {
    "Write the response in English."
  }

  system_msg <- paste(
    "You are a plant genetics assistant.",
    "Use ONLY the JSON evidence provided.",
    "Do not recommend genes not in the list.",
    "Cite evidence snippets when explaining each gene.",
    "Always include: (1) keyword/annotation match notes,",
    "(2) SnpEff impact notes when present,",
    "(3) a short validity comment for each gene.",
    "For keyword/annotation notes:",
    "- Use matched_keywords, matched_synonym, matched_related, unmatched_keywords,",
    "  and matched_context from annotation.details when present.",
    "- Link synonym/related hits to from_user. Treat related hits as indirect",
    "  lexical links (not direct phenotype proof). Treat unmatched_keywords as possible",
    "  counter-evidence (user phrases whose whole family missed).",
    "- Do not list missed synonym/related strings. Do not mention hit rates or high/low.",
    "- Never print numeric scores (composite_score, score_*, keyword_score).",
    "- Do NOT say that a 'semantic match' was confirmed or cite LSA/embedding similarity.",
    "- Ignore matches of prepositions, conjunctions, and other non-biological function words (e.g. in, of, to, and); do not treat them as evidence.",
    "For validity comments: SnpEff/GWAS contradictions only; do not restate keywords.",
    lang_note
  )

  user_msg <- jsonlite::toJSON(
    list(
      phenotype_query = unclass(query),
      ranked_genes = evidence_bundle
    ),
    auto_unbox = TRUE,
    pretty = TRUE
  )

  llmChat(
    messages = list(
      list(role = "system", content = system_msg),
      list(role = "user", content = user_msg)
    ),
    model = llm$model,
    base_url = llm$base_url,
    timeout = llm$timeout
  )
}

#' Generate a Phase 2 B2/B4 AI candidate-gene HTML report
#'
#' Ranks candidates (with required SnpEff + credible-set / finemap channel),
#' then writes an HTML report: basic info, per-peak credible-set comment,
#' scrollable score table, and evidence frames. Report chrome (headings,
#' tables, CS comments) is always English. \code{language = "ja"} only
#' translates LLM-generated per-gene prose (keywords / expression / validity)
#' when \code{use_llm = TRUE}; technical terms and input strings stay unchanged.
#' Ranking scores appear only in the table. PIP / distance / -log10P /
#' SnpEff counts are written by code. Evidence verification is not included
#' yet (Phase 2 B5 deferred).
#'
#' @param object A \code{LazyGas} object with candidate genes.
#' @param pheno Phenotype name or index.
#' @param query Character phenotype description or [PhenotypeQuery]. If
#'   \code{NULL}, uses \code{pheno} name as the query text.
#' @param language \code{"en"} (default) or \code{"ja"}. \code{"ja"} requests
#'   translation of LLM evidence prose only; chrome stays English. No
#'   translation when \code{use_llm = FALSE}.
#' @param top_n Number of top genes for the report body
#'   (default 10). Ranking itself keeps every scored gene.
#' @param sources,weights Ranking channels. \code{finemap} is always added for
#'   this report (weight default 0.5).
#' @param use_llm Use local LLM for evidence prose (and optionally query parse).
#' @param query_use_llm Pass \code{use_llm} to [phenotypeQuery()] when
#'   \code{query} is character.
#' @param model,base_url,timeout LLM connection (default timeout 600s).
#' @param out_dir Optional report directory override (default [aiReportDir()]).
#' @param rank_result Optional precomputed [rankPhenotypeCandidates()] output,
#'   or a path to a saved \code{.rds}/\code{.csv} rank table. Skips re-ranking
#'   when provided (still requires GFF + credible sets for reported peaks;
#'   SnpEff is optional).
#' @param peak_id Optional peak ID(s) to report. When \code{rank_result} is
#'   \code{NULL}, candidates are ranked within each peak; default \code{NULL}
#'   uses the top five peaks by lead \code{negLog10P}.
#' @param use_llm_relevance Annotation LLM relevance during ranking.
#' @param download_tenor Download missing TENOR expression CSVs.
#' @param candidate Optional candidate \code{data.frame} override (must include
#'   \code{Gene_ID} / peak columns as produced by \code{listCandidate}). When
#'   \code{NULL}, candidates are loaded from the companion store.
#' @param gff Required \code{GRanges} gene models (same object as
#'   \code{listCandidate(..., gff = ...)}). Credible sets are also required;
#'   SnpEff is optional.
#' @param work_dir Optional Phase 2 E working directory. When set, evidence uses
#'   the gene×source E pipeline and reads \code{config.yaml}
#'   (\code{mapping_mode}, \code{top_n}). When \code{NULL}, legacy peak-batch
#'   evidence prose is used.
#' @param ann Optional functional annotation table for E0/E1/E5–E7.
#' @param snpeff Optional SnpEff table for E2 / GWAS gene assignment.
#' @param expr_mat Optional expression matrix (genes × samples) for E4.
#' @param interpro Optional InterProScan interval table for E2 domain overlap.
#' @param reinterpret If \code{TRUE} (default), build Phase-2 F hypothesis
#'   cards and select evidence genes via \code{gene_select}.
#' @param reinterpret_rules,reinterpret_views Passed to
#'   [reinterpretPeakCandidates()].
#' @param gene_select Evidence gene set when \code{reinterpret = TRUE}:
#'   \code{"union_champions"} (default), \code{"composite_top_n"}, or
#'   \code{"union_all"}.
#' @param ... Passed to [rankPhenotypeCandidates()] when ranking.
#'   \code{top_n} and \code{save} in \code{...} are ignored.
#'
#' @return A list with \code{html}, \code{markdown} (same HTML string for
#'   compatibility), \code{path}, \code{meta_path}, \code{rank_csv},
#'   \code{rank_rds}, \code{verification} (always \code{NULL} until B5),
#'   \code{rank_result}, \code{reinterpret}, and \code{query}.
#' @export
#'
#' @seealso [rankPhenotypeCandidates()], [reinterpretPeakCandidates()],
#'   [explainPhenotypeCandidates()]
llm_report <- function(object,
                       pheno,
                       query = NULL,
                       language = c("en", "ja"),
                       top_n = 10L,
                       sources = phase1DefaultSources(),
                       weights = phase1DefaultWeights(),
                       use_llm = TRUE,
                       query_use_llm = TRUE,
                       model = NULL,
                       base_url = NULL,
                       timeout = NULL,
                       out_dir = NULL,
                       rank_result = NULL,
                       peak_id = NULL,
                       use_llm_relevance = NULL,
                       download_tenor = TRUE,
                       candidate = NULL,
                       gff = NULL,
                       work_dir = NULL,
                       ann = NULL,
                       snpeff = NULL,
                       expr_mat = NULL,
                       interpro = NULL,
                       reinterpret = TRUE,
                       reinterpret_rules = reinterpretDefaultRules(),
                       reinterpret_views = reinterpretDefaultViews(),
                       gene_select = c("union_champions", "composite_top_n", "union_all"),
                       ...) {
  language <- match.arg(language)
  gene_select <- match.arg(gene_select)
  reinterpret <- isTRUE(reinterpret)
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }
  use_e_path <- !is.null(work_dir) && nzchar(as.character(work_dir)[1L])
  e_cfg <- NULL
  if (use_e_path) {
    work_dir <- path.expand(as.character(work_dir)[1L])
    if (!dir.exists(work_dir)) {
      dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    }
    cfg_path <- file.path(work_dir, "config.yaml")
    if (!file.exists(cfg_path)) {
      stop(
        "work_dir is set but config.yaml is missing. ",
        "Call writeLazyGasExploreConfig(work_dir, mapping_mode=...) first.",
        call. = FALSE
      )
    }
    e_cfg <- readLazyGasExploreConfig(work_dir)
  }
  if (is.character(rank_result) && length(rank_result) == 1L &&
      nzchar(rank_result)) {
    rank_result <- .read_rank_result(rank_result)
  } else if (!is.null(rank_result) && !is.data.frame(rank_result)) {
    stop(
      "'rank_result' must be NULL, a data.frame, or a .rds/.csv path.",
      call. = FALSE
    )
  }

  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  llm <- .ai_config_llm_settings(model = model, base_url = base_url, timeout = timeout)
  sw <- .report_sources_with_finemap(sources = sources, weights = weights)
  sources <- sw$sources
  weights <- sw$weights

  if (is.null(query)) {
    query <- pheno_name
  }
  if (is.character(query)) {
    query <- phenotypeQuery(
      text = query,
      use_llm = isTRUE(query_use_llm) && isTRUE(use_llm),
      llm_model = llm$model,
      llm_base_url = llm$base_url,
      llm_timeout = llm$timeout,
      object = object,
      save = TRUE
    )
  } else if (!inherits(query, "PhenotypeQuery")) {
    stop("'query' must be NULL, character, or PhenotypeQuery.", call. = FALSE)
  }

  ranked_in_this_call <- FALSE
  if (is.null(rank_result)) {
    peaks_plan <- .resolve_report_peaks(
      object = object,
      pheno_name = pheno_name,
      peak_id = peak_id
    )
    message(
      "llm_report: ranking ", nrow(peaks_plan), " peak(s): ",
      paste(peaks_plan$peak_ID, collapse = ", ")
    )
    .report_require_cs_gff(
      object = object,
      pheno_name = pheno_name,
      peak_ids = peaks_plan$peak_ID,
      gff = gff
    )
    candidate_all <- if (!is.null(candidate)) {
      candidate
    } else {
      lazyData(
        object = object,
        dataset = "candidate",
        pheno = pheno_name
      )
    }
    if (is.null(candidate_all) || nrow(candidate_all) == 0L) {
      stop("No candidate data found for phenotype '", pheno_name, "'.",
           call. = FALSE)
    }
    if (!use_e_path && !any(c("HIGH", "MODERATE") %in% names(candidate_all))) {
      stop(
        "Candidate table lacks SnpEff impact columns. ",
        "Run listCandidate(..., snpeff = ...) first.",
        call. = FALSE
      )
    }
    extra_rank <- list(...)
    extra_rank$top_n <- NULL
    extra_rank$save <- NULL
    rank_args <- list(
      object = object,
      pheno = pheno_name,
      query = query,
      sources = sources,
      weights = weights,
      save = FALSE,
      use_llm_relevance = use_llm_relevance,
      download_tenor = download_tenor,
      query_use_llm = FALSE,
      llm_model = llm$model,
      llm_base_url = llm$base_url,
      llm_timeout = llm$timeout
    )
    rank_args <- c(rank_args, extra_rank)
    rank_parts <- list()
    for (i in seq_len(nrow(peaks_plan))) {
      pid <- peaks_plan$peak_ID[i]
      message(
        "llm_report: rankPhenotypeCandidates peak ", pid,
        " (", i, "/", nrow(peaks_plan), ")..."
      )
      cand_peak <- candidate_all[
        candidate_all$peak_ID == pid,
        ,
        drop = FALSE
      ]
      if (nrow(cand_peak) == 0L) {
        message("llm_report: peak ", pid, " has no candidates; skip")
        next
      }
      if ("Gene_ID" %in% names(cand_peak)) {
        gid_ok <- !is.na(cand_peak$Gene_ID) & nzchar(as.character(cand_peak$Gene_ID))
        if (!any(gid_ok)) {
          message("llm_report: peak ", pid, " has no Gene_ID rows; skip")
          next
        }
        cand_peak <- cand_peak[gid_ok, , drop = FALSE]
      }
      rank_peak <- tryCatch(
        do.call(
          rankPhenotypeCandidates,
          c(
            rank_args,
            list(
              candidate = cand_peak,
              top_n = NULL
            )
          )
        ),
        error = function(e) {
          warning(
            "Ranking failed for peak ", pid, ": ", conditionMessage(e),
            " — skipping this peak.",
            call. = FALSE
          )
          NULL
        }
      )
      if (!is.null(rank_peak)) {
        message(
          "llm_report: ranked peak ", pid, " — ",
          nrow(rank_peak), " gene(s)"
        )
        rank_parts[[length(rank_parts) + 1L]] <- rank_peak
      }
    }
    if (!length(rank_parts)) {
      stop(
        "Ranking failed for all requested peaks; no report written.",
        call. = FALSE
      )
    }
    rank_result <- .bind_rank_results(rank_parts)
    ranked_in_this_call <- TRUE
    message("llm_report: ranking done")
  } else {
    peaks_plan <- .resolve_report_peaks(
      object = object,
      pheno_name = pheno_name,
      peak_id = if (!is.null(peak_id)) {
        peak_id
      } else if ("peak_ID" %in% names(rank_result)) {
        unique(rank_result$peak_ID)
      } else {
        NULL
      }
    )
    .report_require_cs_gff(
      object = object,
      pheno_name = pheno_name,
      peak_ids = peaks_plan$peak_ID,
      gff = gff
    )
  }

  basic_html <- .report_basic_info_html(
    object = object,
    pheno_name = pheno_name,
    language = language
  )
  if (use_e_path) {
    basic_html <- paste0(
      basic_html,
      "<p><em>B5 numeric grounding checks digits against each frame payload only; ",
      "non-numeric overclaims are not verified. ",
      "Rewrite failures are listed under each peak section.</em></p>\n"
    )
  }
  reinterpret_obj <- NULL
  if (reinterpret && !is.null(rank_result) && nrow(rank_result) > 0L) {
    message("llm_report: reinterpretPeakCandidates...")
    reinterpret_obj <- tryCatch(
      reinterpretPeakCandidates(
        ranked = rank_result,
        query = query,
        peak_id = peaks_plan$peak_ID,
        views = reinterpret_views,
        rules = reinterpret_rules,
        save = TRUE,
        object = object,
        pheno = pheno_name
      ),
      error = function(e) {
        warning(
          "reinterpretPeakCandidates failed: ", conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
    if (!is.null(reinterpret_obj)) {
      message("llm_report: reinterpret done")
    }
  }
  e_b5_all <- list()
  sections <- character()
  e_gff_tables <- NULL
  if (use_e_path) {
    protein_map_store_once <- tryCatch(
      lazyData(object = object, dataset = "gene_protein_map"),
      error = function(e) NULL
    )
    e_gff_tables <- .gff_e_tables(gff, protein_map = protein_map_store_once)
  }
  for (i in seq_len(nrow(peaks_plan))) {
    pid <- peaks_plan$peak_ID[i]
    message(
      "llm_report: peak section ", pid,
      " (", i, "/", nrow(peaks_plan), ")..."
    )
    rank_peak <- if (!is.null(rank_result) && "peak_ID" %in% names(rank_result)) {
      rank_result[as.character(rank_result$peak_ID) == as.character(pid), , drop = FALSE]
    } else if (!is.null(rank_result) && nrow(peaks_plan) == 1L) {
      rank_result
    } else {
      data.frame()
    }
    if (!is.null(rank_peak) && nrow(rank_peak) > 0L &&
        "composite_score" %in% names(rank_peak)) {
      rank_peak <- rank_peak[order(-rank_peak$composite_score, rank_peak$Gene_ID), ,
                             drop = FALSE]
    }
    reinterpret_peak <- NULL
    evidence_rank <- NULL
    gene_set <- character()
    if (!is.null(reinterpret_obj)) {
      reinterpret_peak <- reinterpret_obj$peaks[[as.character(pid)]]
      top_ids <- if (nrow(rank_peak)) {
        utils::head(as.character(rank_peak$Gene_ID), as.integer(top_n)[1L])
      } else {
        character()
      }
      gene_set <- geneSetForReport(
        reinterpret_obj,
        peak_id = pid,
        mode = gene_select,
        composite_ids = top_ids
      )
      if (length(gene_set) && nrow(rank_peak)) {
        evidence_rank <- rank_peak[
          as.character(rank_peak$Gene_ID) %in% gene_set,
          ,
          drop = FALSE
        ]
      }
    }
    e_run <- NULL
    if (use_e_path) {
      message("llm_report: E evidence path for peak ", pid, "...")
      cred <- tryCatch(
        lazyData(
          object = object,
          dataset = "credible_set",
          pheno = pheno_name,
          kind = paste0("peak_", pid)
        ),
        error = function(e) NULL
      )
      summ <- attr(cred, "summary")
      resolution <- if (!is.null(summ$resolution)) {
        as.character(summ$resolution)[1L]
      } else {
        "low"
      }
      simple <- tryCatch(
        lazyData(object = object, dataset = "simple_candidate", pheno = pheno_name),
        error = function(e) NULL
      )
      if (!is.null(simple) && nrow(simple) && "peak_ID" %in% names(simple)) {
        simple <- simple[as.character(simple$peak_ID) == as.character(pid), , drop = FALSE]
      }
      if (is.null(simple) || !nrow(simple)) {
        simple <- if (!is.null(rank_peak) && nrow(rank_peak)) {
          rank_peak
        } else if (!is.null(candidate)) {
          candidate[as.character(candidate$peak_ID) == as.character(pid), , drop = FALSE]
        } else {
          data.frame()
        }
      }
      if (length(gene_set) && nrow(simple) && "Gene_ID" %in% names(simple)) {
        # Prefer union_champions / gene_select set; fall back if empty intersection
        simple_f <- simple[as.character(simple$Gene_ID) %in% gene_set, , drop = FALSE]
        if (nrow(simple_f)) {
          simple <- simple_f
        }
      }
      if (nrow(simple) && !"Chr" %in% names(simple) && "Gene_chr" %in% names(simple)) {
        simple$Chr <- simple$Gene_chr
      }
      # E top_n from config; legacy report top_n only for candidate table display
      e_top <- e_cfg$top_n %||% Inf
      kw_hits <- .keyword_hits_from_rank_result(
        if (!is.null(evidence_rank) && nrow(evidence_rank)) {
          evidence_rank
        } else {
          rank_peak
        }
      )
      snpeff_peak <- .snpeff_for_peak(
        snpeff = snpeff,
        object = object,
        pheno_name = pheno_name,
        peak_id = pid
      )
      message(
        "llm_report: SnpEff peak-block markers for peak ", pid, ": ",
        nrow(snpeff_peak), " annotation row(s)"
      )
      e_run <- .report_e_run_peak(
        work_dir = work_dir,
        mapping_mode = e_cfg$mapping_mode,
        simple_candidates = simple,
        credible_set = .e_enrich_cs_coords(
          cred,
          object = object,
          pheno_name = pheno_name,
          peak_id = pid
        ),
        gff = gff,
        ann = ann,
        snpeff = snpeff_peak,
        interpro = interpro,
        expr_mat = expr_mat,
        trait_text = query$trait_text %||% query$text %||% pheno_name,
        keyword_hits_by_gene = kw_hits,
        peak_id = pid,
        use_llm = use_llm,
        llm = llm,
        resolution = resolution,
        top_n = e_top,
        cfg = e_cfg,
        gff_windows = e_gff_tables$windows,
        protein_map = e_gff_tables$protein_map
      )
      if (length(e_run$b5_failures)) {
        e_b5_all <- c(e_b5_all, e_run$b5_failures)
      }
    }
    sections <- c(
      sections,
      .report_peak_section_html(
        object = object,
        pheno_name = pheno_name,
        peak_id = pid,
        chr = peaks_plan$Chr[i],
        pos = peaks_plan$Pos[i],
        region_start = peaks_plan$region_start[i],
        region_end = peaks_plan$region_end[i],
        rank_peak = rank_peak,
        sources = sources,
        query = query,
        top_n = top_n,
        language = language,
        use_llm = use_llm,
        llm = llm,
        e_run = e_run,
        reinterpret_peak = reinterpret_peak,
        evidence_rank = evidence_rank
      )
    )
  }

  html_body <- paste(c(basic_html, sections), collapse = "\n")
  # Evidence verification omitted for now (Phase 2 B5 deferred as section;
  # per-source B5 grounding runs inside E when work_dir is set).
  html <- .report_html_document(
    body_html = html_body,
    title = paste("lazyGas AI report —", pheno_name)
  )

  report_path <- NULL
  meta_path <- NULL
  rank_csv_path <- NULL
  rank_rds_path <- NULL
  report_dir <- aiReportDir(object = object, out_dir = out_dir)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe_pheno <- .store_safe_name(pheno_name)
  base <- paste0(safe_pheno, "_", stamp)
  report_path <- file.path(report_dir, paste0(base, ".html"))
  meta_path <- file.path(report_dir, paste0(base, ".meta.json"))
  message("llm_report: writing HTML ", report_path)
  writeLines(html, report_path, useBytes = TRUE)
  if (!is.null(rank_result) && nrow(rank_result) > 0L) {
    rank_csv_path <- file.path(report_dir, paste0(base, "_rank.csv"))
    rank_rds_path <- file.path(report_dir, paste0(base, "_rank.rds"))
    .write_rank_result_files(
      rank_result = rank_result,
      csv_path = rank_csv_path,
      rds_path = rank_rds_path
    )
    if (isTRUE(ranked_in_this_call)) {
      .store_write_phenotype_rank(
        object = object,
        rank_result = rank_result,
        pheno_name = pheno_name,
        query_id = query$query_id
      )
    }
  }
  jsonlite::write_json(
    list(
      pheno = pheno_name,
      query = unclass(query),
      peak_id = peak_id,
      language = language,
      top_n = as.integer(top_n)[1L],
      sources = sources,
      weights = as.list(weights),
      use_llm = isTRUE(use_llm),
      model = llm$model,
      base_url = llm$base_url,
      created_at = stamp,
      report_file = basename(report_path),
      rank_csv = if (is.null(rank_csv_path)) NULL else basename(rank_csv_path),
      rank_rds = if (is.null(rank_rds_path)) NULL else basename(rank_rds_path),
      rank_meta = attr(rank_result, "phenotypeRank"),
      report_skeleton = if (use_e_path) "phase2_E" else "phase2_B4",
      reinterpret = reinterpret,
      gene_select = gene_select,
      work_dir = if (use_e_path) work_dir else NULL,
      mapping_mode = if (use_e_path) e_cfg$mapping_mode else NULL,
      b5_failures = if (length(e_b5_all)) names(e_b5_all) else list()
    ),
    meta_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  if (use_e_path) {
    row_sd_med <- if (!is.null(expr_mat)) {
      .e4_row_sd_median(expr_mat)
    } else {
      NA_real_
    }
    z_warn <- if (!is.null(expr_mat)) {
      .e4_warn_row_z(
        expr_mat,
        e_cfg$row_sd_z_warn_lo %||% 0.8,
        e_cfg$row_sd_z_warn_hi %||% 1.2
      )
    } else {
      NULL
    }
    .write_run_metadata(
      work_dir,
      list(
        pheno = pheno_name,
        mapping_mode = e_cfg$mapping_mode,
        top_n = e_cfg$top_n,
        tau_high = e_cfg$tau_high %||% 0.8,
        tau_moderate = e_cfg$tau_moderate %||% 0.4,
        relative_high = e_cfg$relative_high %||% 0.8,
        relative_moderate = e_cfg$relative_moderate %||% 0.4,
        matrix_scale = e_cfg$matrix_scale %||% "abundance",
        row_sd_median = row_sd_med,
        row_z_warning = z_warn,
        row_sd_z_warn_lo = e_cfg$row_sd_z_warn_lo %||% 0.8,
        row_sd_z_warn_hi = e_cfg$row_sd_z_warn_hi %||% 1.2,
        use_llm = isTRUE(use_llm),
        model = llm$model,
        report_file = basename(report_path),
        b5_failures = as.list(names(e_b5_all)),
        b5_numeric_only = TRUE,
        b5_note = paste(
          "B5 verifies numeric tokens against payloads only;",
          "non-numeric overclaims are not checked."
        ),
        created_at = stamp,
        config = e_cfg
      )
    )
  }

  list(
    html = html,
    markdown = html,
    path = report_path,
    meta_path = meta_path,
    rank_csv = rank_csv_path,
    rank_rds = rank_rds_path,
    verification = NULL,
    rank_result = rank_result,
    reinterpret = reinterpret_obj,
    query = query
  )
}

#' Bind per-peak rank tables and restore phenotypeRank metadata
#' @keywords internal
.bind_rank_results <- function(rank_parts) {
  out <- do.call(rbind, rank_parts)
  rownames(out) <- NULL
  meta <- NULL
  for (part in rank_parts) {
    meta <- attr(part, "phenotypeRank")
    if (!is.null(meta)) {
      break
    }
  }
  if (is.null(meta)) {
    meta <- list()
  }
  meta$n_genes <- nrow(out)
  attr(out, "phenotypeRank") <- meta
  out
}

#' Write a rank table as CSV (inspection) and RDS (reuse)
#' @keywords internal
.write_rank_result_files <- function(rank_result, csv_path, rds_path) {
  utils::write.csv(
    rank_result,
    file = csv_path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  saveRDS(rank_result, rds_path)
  invisible(list(csv = csv_path, rds = rds_path))
}

#' Load a saved rank table from .rds or .csv
#' @keywords internal
.read_rank_result <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(trimws(path))) {
    stop("'rank_result' path must be a single non-empty character string.",
         call. = FALSE)
  }
  path <- path.expand(trimws(path))
  if (!file.exists(path)) {
    stop("rank_result file not found: ", path, call. = FALSE)
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    out <- readRDS(path)
  } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    out <- utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8"
    )
  } else {
    stop("rank_result file must be .rds or .csv: ", path, call. = FALSE)
  }
  if (!is.data.frame(out) || nrow(out) == 0L) {
    stop("rank_result file did not contain a non-empty data.frame: ", path,
         call. = FALSE)
  }
  out
}

#' Format genomic positions for report headings
#' @keywords internal
.report_fmt_pos <- function(x) {
  format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

#' Basic-info Markdown block for llm_report (Phase 2 B1)
#' @keywords internal
.report_basic_info_md <- function(object, pheno_name, language = c("ja", "en")) {
  language <- match.arg(language)
  gds <- .store_gds_fn(object)
  companion <- .store_path(object)
  if (is.null(companion) || !nzchar(as.character(companion)[1L])) {
    companion <- "(none)"
  } else {
    companion <- normalizePath(as.character(companion)[1L], winslash = "/",
                               mustWork = FALSE)
  }
  n_sam <- as.integer(nsam(object))[1L]
  n_mar <- as.integer(nmar(object))[1L]
  pheno_df <- getPheno(object)$pheno
  pheno_vec <- if (!is.null(pheno_df) && pheno_name %in% names(pheno_df)) {
    pheno_df[[pheno_name]]
  } else {
    rep(NA, n_sam)
  }
  n_obs <- sum(!is.na(pheno_vec))
  pct <- if (is.finite(n_sam) && n_sam > 0L) {
    100 * n_obs / n_sam
  } else {
    NA_real_
  }

  if (identical(language, "ja")) {
    paste(
      c(
        "## 基礎情報",
        "",
        "- 入力ファイル:",
        paste0("  - GDS: `", gds, "`"),
        paste0("  - companion: `", companion, "`"),
        paste0("- サンプル数: ", n_sam),
        paste0("- マーカー数: ", n_mar),
        paste0("- 表現型: ", pheno_name),
        sprintf("- 表現型非欠損: %s / %s（%.1f%%）", n_obs, n_sam, pct),
        ""
      ),
      collapse = "\n"
    )
  } else {
    paste(
      c(
        "## Basic information",
        "",
        "- Input files:",
        paste0("  - GDS: `", gds, "`"),
        paste0("  - companion: `", companion, "`"),
        paste0("- Samples: ", n_sam),
        paste0("- Markers: ", n_mar),
        paste0("- Phenotype: ", pheno_name),
        sprintf("- Non-missing phenotype: %s / %s (%.1f%%)", n_obs, n_sam, pct),
        ""
      ),
      collapse = "\n"
    )
  }
}

#' Per-peak Markdown skeleton (Phase 2 B1 placeholders)
#' @keywords internal
.report_peak_section_md <- function(peak_id,
                                    chr,
                                    pos,
                                    region_start,
                                    region_end,
                                    language = c("ja", "en")) {
  language <- match.arg(language)
  header <- sprintf(
    "## Peak %s — %s:%s（%s-%s）",
    as.character(peak_id)[1L],
    as.character(chr)[1L],
    .report_fmt_pos(pos),
    .report_fmt_pos(region_start),
    .report_fmt_pos(region_end)
  )
  if (identical(language, "ja")) {
    paste(
      c(
        header,
        "",
        "### 候補一覧",
        "",
        "（未実装）",
        "",
        "### 根拠",
        "",
        "（未実装）",
        ""
      ),
      collapse = "\n"
    )
  } else {
    paste(
      c(
        header,
        "",
        "### Candidate list",
        "",
        "（未実装）",
        "",
        "### Evidence",
        "",
        "（未実装）",
        ""
      ),
      collapse = "\n"
    )
  }
}

#' Region start/end per peak_ID from a peakcall table
#' @keywords internal
.peakcall_region_bounds <- function(peakcall) {
  if (is.null(peakcall) || !nrow(peakcall) || !"peak_ID" %in% names(peakcall) ||
      !"Pos" %in% names(peakcall)) {
    return(data.frame(
      peak_ID = integer(),
      region_start = numeric(),
      region_end = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  parts <- split(peakcall, as.character(peakcall$peak_ID))
  do.call(
    rbind,
    lapply(parts, function(block) {
      pos <- as.numeric(block$Pos)
      pos <- pos[is.finite(pos)]
      data.frame(
        peak_ID = block$peak_ID[1L],
        region_start = if (length(pos)) min(pos) else NA_real_,
        region_end = if (length(pos)) max(pos) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
  )
}

#' Resolve peak table rows for per-peak LLM reports
#' @keywords internal
.resolve_report_peaks <- function(object, pheno_name, peak_id = NULL) {
  peakcall <- .get_peakcall(object = object, pheno_name = pheno_name, recalc = TRUE)
  if (is.null(peakcall) || nrow(peakcall) == 0L) {
    peakcall <- .get_peakcall(object = object, pheno_name = pheno_name, recalc = FALSE)
  }
  if (is.null(peakcall) || nrow(peakcall) == 0L) {
    stop("No peakcall data for phenotype '", pheno_name, "'.", call. = FALSE)
  }

  regions <- .peakcall_region_bounds(peakcall)

  leads <- peakcall[
    peakcall$peak_variant_ID == peakcall$variant_ID,
    ,
    drop = FALSE
  ]
  if (nrow(leads) == 0L) {
    leads <- do.call(rbind, lapply(split(peakcall, peakcall$peak_ID), function(block) {
      block[which.max(block$negLog10P), , drop = FALSE]
    }))
  }
  ord <- order(-leads$negLog10P, leads$peak_ID)
  leads <- leads[ord, , drop = FALSE]
  leads <- leads[!duplicated(leads$peak_ID), , drop = FALSE]
  out <- data.frame(
    peak_ID = leads$peak_ID,
    Chr = leads$Chr,
    Pos = leads$Pos,
    stringsAsFactors = FALSE
  )
  out <- merge(out, regions, by = "peak_ID", all.x = TRUE, sort = FALSE)
  # restore lead negLog10P order after merge
  out <- out[match(as.character(leads$peak_ID), as.character(out$peak_ID)), , drop = FALSE]
  rownames(out) <- NULL

  if (is.null(peak_id)) {
    out <- out[seq_len(min(5L, nrow(out))), , drop = FALSE]
    return(out)
  }

  peak_id <- as.character(peak_id)
  avail <- as.character(out$peak_ID)
  miss <- setdiff(peak_id, avail)
  if (length(miss)) {
    stop(
      "peak_id not found for phenotype '", pheno_name, "': ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  out <- out[match(peak_id, avail), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.verify_report_against_evidence <- function(markdown,
                                            rank_result,
                                            top_n = 10L,
                                            language = "ja") {
  # Deferred: not wired into llm_report() (Phase 2 B5 unresolved).
  top_n <- as.integer(top_n)[1L]
  sub <- rank_result[seq_len(min(top_n, nrow(rank_result))), , drop = FALSE]
  bundle <- .phenotype_evidence_bundle(sub)

  # Collect searchable evidence tokens per gene
  evidence_text <- setNames(character(nrow(sub)), sub$Gene_ID)
  for (i in seq_along(bundle)) {
    item <- bundle[[i]]
    bits <- character()
    for (src in names(item$evidence %||% list())) {
      block <- item$evidence[[src]]
      if (!is.null(block$snippets)) {
        bits <- c(bits, as.character(unlist(block$snippets)))
      }
      if (!is.null(block$details)) {
        bits <- c(bits, as.character(unlist(block$details)))
      }
    }
    evidence_text[[item$Gene_ID]] <- tolower(paste(bits, collapse = " "))
  }

  # Sentence-like claims
  text <- gsub("\r", "", markdown)
  parts <- unlist(strsplit(text, "(?<=[。．\\.\\!\\?\\n])\\s*", perl = TRUE))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts) & nchar(parts) > 12L]

  gene_ids <- as.character(sub$Gene_ID)
  claims <- list()
  n_checked <- 0L
  n_supported <- 0L
  n_marked <- 0L

  for (sent in parts) {
    hit_genes <- gene_ids[vapply(gene_ids, function(g) {
      grepl(g, sent, fixed = TRUE)
    }, logical(1))]
    if (!length(hit_genes)) {
      next
    }
    n_checked <- n_checked + 1L
    supported <- FALSE
    for (g in hit_genes) {
      ev <- evidence_text[[g]]
      # Token overlap: any alphanumeric token length >= 4 from claim in evidence
      toks <- unique(unlist(regmatches(
        tolower(sent),
        gregexpr("[a-z0-9_]{4,}", tolower(sent), perl = TRUE)
      )))
      toks <- setdiff(toks, tolower(gene_ids))
      if (!length(toks)) {
        # Gene mentioned with score-like context counts as supported if gene in rank
        supported <- TRUE
        break
      }
      if (any(vapply(toks, function(t) grepl(t, ev, fixed = TRUE), logical(1)))) {
        supported <- TRUE
        break
      }
    }
    if (supported) {
      n_supported <- n_supported + 1L
      claims[[length(claims) + 1L]] <- list(
        text = sent,
        genes = hit_genes,
        status = "supported"
      )
    } else {
      n_marked <- n_marked + 1L
      claims[[length(claims) + 1L]] <- list(
        text = sent,
        genes = hit_genes,
        status = "unsupported_marked"
      )
    }
  }

  congruence <- if (n_checked > 0L) n_supported / n_checked else NA_real_
  list(
    n_claims_checked = n_checked,
    n_supported = n_supported,
    n_unsupported_marked = n_marked,
    congruence_rate = congruence,
    claims = claims
  )
}

.append_verification_section <- function(markdown, verification, language = "ja") {
  rate <- verification$congruence_rate
  rate_txt <- if (is.finite(rate)) sprintf("%.1f%%", 100 * rate) else "NA"
  marked <- Filter(function(x) identical(x$status, "unsupported_marked"),
                   verification$claims %||% list())

  if (language == "ja") {
    hdr <- c(
      "",
      "---",
      "",
      "## 根拠検証",
      "",
      sprintf("- 照合クレーム数: %s", verification$n_claims_checked),
      sprintf("- 支持されたクレーム: %s", verification$n_supported),
      sprintf("- 根拠なし（マークのみ）: %s", verification$n_unsupported_marked),
      sprintf("- 整合率: %s", rate_txt),
      ""
    )
    if (length(marked)) {
      hdr <- c(hdr, "### 根拠なしとしてマークした記述", "")
      for (m in marked) {
        hdr <- c(hdr, paste0("- ⚠️ ", m$text))
      }
      hdr <- c(hdr, "")
    }
  } else {
    hdr <- c(
      "",
      "---",
      "",
      "## Evidence verification",
      "",
      sprintf("- Claims checked: %s", verification$n_claims_checked),
      sprintf("- Supported: %s", verification$n_supported),
      sprintf("- Unsupported (marked only): %s", verification$n_unsupported_marked),
      sprintf("- Congruence rate: %s", rate_txt),
      ""
    )
    if (length(marked)) {
      hdr <- c(hdr, "### Marked unsupported claims", "")
      for (m in marked) {
        hdr <- c(hdr, paste0("- ⚠️ ", m$text))
      }
      hdr <- c(hdr, "")
    }
  }
  paste0(markdown, "\n", paste(hdr, collapse = "\n"))
}

#' Answer a follow-up question about ranked phenotype candidates
#'
#' Unlike [explainPhenotypeCandidates()], this function responds to a specific
#' user question (e.g. gene comparison) instead of dumping a full evidence
#' report. Uses a local LLM when available; otherwise a rule-based answer.
#'
#' @param rank_result Output of [rankPhenotypeCandidates()].
#' @param question User question (character).
#' @param query Optional \code{PhenotypeQuery} object.
#' @param object Optional \code{LazyGas} object to load saved queries.
#' @param top_n Number of ranked genes to include as context.
#' @param use_llm Use local LLM when available.
#' @param language \code{"ja"} or \code{"en"}.
#' @param model,base_url LLM settings.
#' @param chat_history Optional prior chat messages (list of \code{role}/\code{content}).
#'
#' @return A character string answering the question.
#' @export
answerPhenotypeQuestion <- function(rank_result,
                                    question,
                                    query = NULL,
                                    object = NULL,
                                    top_n = 10L,
                                    use_llm = TRUE,
                                    language = c("ja", "en"),
                                    model = NULL,
                                    base_url = NULL,
                                    chat_history = NULL) {
  language <- match.arg(language)
  if (!is.character(question) || length(question) != 1L || !nzchar(trimws(question))) {
    stop("'question' must be a non-empty character string.", call. = FALSE)
  }
  question <- trimws(question)

  if (!is.data.frame(rank_result) || nrow(rank_result) == 0L) {
    stop("'rank_result' must be a non-empty data.frame.", call. = FALSE)
  }

  meta <- attr(rank_result, "phenotypeRank")
  if (is.null(query) && !is.null(meta$query_id) && !is.null(object)) {
    query <- .store_read_phenotype_query(object, meta$query_id)
  }
  if (is.null(query)) {
    query <- list(
      trait_text = "phenotype query",
      query_id = if (!is.null(meta)) meta$query_id else "unknown"
    )
    class(query) <- c("PhenotypeQuery", "list")
  }

  top_n <- as.integer(top_n)[1L]
  sub <- rank_result[seq_len(min(top_n, nrow(rank_result))), , drop = FALSE]
  evidence_bundle <- .phenotype_evidence_bundle(rank_result = sub)

  if (isTRUE(use_llm) && llmHealthCheck(base_url = base_url, timeout = 5)) {
    lang_note <- if (language == "ja") {
      paste(
        "ユーザーの質問に直接答えてください。",
        "全文レポートや箇条書きの証拠ダンプは不要です。",
        "比較質問なら結論を最初に述べ、根拠を短く説明してください。",
        "日本語で回答してください。"
      )
    } else {
      paste(
        "Answer the user's question directly.",
        "Do not dump a full evidence report.",
        "For comparison questions, state the conclusion first, then brief reasons.",
        "Write in English."
      )
    }

    messages <- list(
      list(
        role = "system",
        content = paste(
          "You are a plant genetics assistant.",
          "Use ONLY the JSON evidence provided.",
          "Do not recommend genes not in the list.",
          lang_note
        )
      )
    )
    if (!is.null(chat_history) && length(chat_history)) {
      for (msg in chat_history) {
        if (!is.null(msg$role) && !is.null(msg$content) && nzchar(msg$content)) {
          messages[[length(messages) + 1L]] <- list(
            role = msg$role,
            content = as.character(msg$content)
          )
        }
      }
    }
    messages[[length(messages) + 1L]] <- list(
      role = "user",
      content = jsonlite::toJSON(
        list(
          question = question,
          phenotype_query = unclass(query),
          ranked_genes = evidence_bundle
        ),
        auto_unbox = TRUE,
        pretty = TRUE
      )
    )

    return(tryCatch(
      llmChat(
        messages = messages,
        model = model,
        base_url = base_url,
        timeout = 180
      ),
      error = function(e) {
        warning("LLM answer failed; using rule-based reply.", call. = FALSE)
        .answer_phenotype_question_rules(
          question = question,
          query = query,
          evidence_bundle = evidence_bundle,
          language = language
        )
      }
    ))
  }

  .answer_phenotype_question_rules(
    question = question,
    query = query,
    evidence_bundle = evidence_bundle,
    language = language
  )
}

.phenotype_evidence_bundle <- function(rank_result) {
  lapply(seq_len(nrow(rank_result)), function(i) {
    row <- rank_result[i, , drop = FALSE]
    ev <- tryCatch(
      jsonlite::fromJSON(row$evidence_json[1L], simplifyVector = FALSE),
      error = function(e) list()
    )
    list(
      rank = i,
      Gene_ID = row$Gene_ID[1L],
      Name = if ("Name" %in% names(row)) row$Name[1L] else NA_character_,
      composite_score = row$composite_score[1L],
      negLog10P = if ("negLog10P" %in% names(row)) row$negLog10P[1L] else NA_real_,
      dist2peak = if ("dist2peak" %in% names(row)) row$dist2peak[1L] else NA_real_,
      scores = list(
        annotation = if ("score_annotation" %in% names(row)) row$score_annotation[1L] else NA_real_,
        snpeff = if ("score_snpeff" %in% names(row)) row$score_snpeff[1L] else NA_real_,
        gwas = if ("score_gwas" %in% names(row)) row$score_gwas[1L] else NA_real_,
        expression = if ("score_expression" %in% names(row)) row$score_expression[1L] else NA_real_,
        literature = if ("score_literature" %in% names(row)) row$score_literature[1L] else NA_real_,
        ortholog = if ("score_ortholog" %in% names(row)) row$score_ortholog[1L] else NA_real_
      ),
      evidence = ev
    )
  })
}

.is_compare_question <- function(question) {
  q <- tolower(question)
  grepl(
    "どちら|どっち|より|比較|信頼|おすすめ|上位|強い|弱い|which|more reliable|better candidate|compare|versus|\\bvs\\b",
    q,
    perl = TRUE
  )
}

.is_ranking_question <- function(question) {
  q <- tolower(question)
  grepl(
    "なぜ|理由|why|how come|rank|ランク|スコア|根拠",
    q,
    perl = TRUE
  )
}

.evidence_snippet_text <- function(item, source) {
  block <- item$evidence[[source]]
  if (is.null(block) || !length(block$snippets)) {
    return(character())
  }
  as.character(unlist(block$snippets))
}

.answer_phenotype_question_rules <- function(question,
                                            query,
                                            evidence_bundle,
                                            language) {
  if (!length(evidence_bundle)) {
    return(if (language == "ja") {
      "候補遺伝子データがありません。先に Rank candidates を実行してください。"
    } else {
      "No ranked genes available. Run Rank candidates first."
    })
  }

  if (.is_compare_question(question)) {
    return(.answer_compare_genes(
      evidence_bundle = evidence_bundle,
      query = query,
      language = language
    ))
  }

  if (.is_ranking_question(question)) {
    return(.answer_ranking_rationale(
      query = query,
      evidence_bundle = evidence_bundle,
      language = language
    ))
  }

  top <- evidence_bundle[[1L]]
  if (language == "ja") {
    paste0(
      "ご質問「", question, "」について:\n",
      "現在のランキング1位は **", top$Gene_ID, "**（総合スコア ",
      sprintf("%.3f", top$composite_score), "）です。",
      "比較の質問には「どちらがより信頼できる？」のように聞いてください。"
    )
  } else {
    paste0(
      "Regarding your question \"", question, "\":\n",
      "The top-ranked gene is **", top$Gene_ID, "** (composite score ",
      sprintf("%.3f", top$composite_score), ").",
      " For comparisons, try asking which gene is more reliable."
    )
  }
}

.answer_compare_genes <- function(evidence_bundle, language, query = NULL) {
  if (length(evidence_bundle) < 2L) {
    g <- evidence_bundle[[1L]]
    if (language == "ja") {
      return(paste0(
        "比較対象は **", g$Gene_ID, "** のみです（総合スコア ",
        sprintf("%.3f", g$composite_score), "）。他候補がランキングに含まれていません。"
      ))
    }
    return(paste0(
      "Only **", g$Gene_ID, "** is available (composite ",
      sprintf("%.3f", g$composite_score), ")."
    ))
  }

  a <- evidence_bundle[[1L]]
  b <- evidence_bundle[[2L]]
  winner <- if (a$composite_score >= b$composite_score) a else b
  loser <- if (identical(winner, a)) b else a

  reasons <- character()
  if (language == "ja") {
    if (winner$scores$annotation > loser$scores$annotation + 0.05) {
      ann_w <- paste(.evidence_snippet_text(winner, "annotation"), collapse = " ")
      trait <- if (!is.null(query$trait_text)) query$trait_text else "表現型"
      reasons <- c(
        reasons,
        paste0(
          "- **機能アノテーション**: ", winner$Gene_ID, " の方が「", trait,
          "」に一致",
          if (nzchar(ann_w)) paste0("（", ann_w, "）") else ""
        )
      )
    }
    if (!is.na(winner$dist2peak) && !is.na(loser$dist2peak)) {
      reasons <- c(
        reasons,
        paste0(
          "- **GWAS 位置**: ", winner$Gene_ID, " dist2peak=",
          format(winner$dist2peak, big.mark = ","), " bp vs ",
          loser$Gene_ID, " ", format(loser$dist2peak, big.mark = ","), " bp"
        )
      )
    } else if (!is.na(winner$negLog10P)) {
      reasons <- c(
        reasons,
        paste0(
          "- **GWAS**: negLog10P ", sprintf("%.2f", winner$negLog10P),
          " vs ", sprintf("%.2f", loser$negLog10P)
        )
      )
    }
    if (winner$scores$expression > loser$scores$expression + 0.01) {
      expr_w <- paste(.evidence_snippet_text(winner, "expression"), collapse = "; ")
      expr_l <- paste(.evidence_snippet_text(loser, "expression"), collapse = "; ")
      reasons <- c(
        reasons,
        paste0("- **発現**: ", winner$Gene_ID, " — ", expr_w, " / ", loser$Gene_ID, " — ", expr_l)
      )
    }
    if (winner$scores$ortholog > loser$scores$ortholog + 0.01) {
      ortho_w <- paste(.evidence_snippet_text(winner, "ortholog"), collapse = ", ")
      ortho_l <- paste(.evidence_snippet_text(loser, "ortholog"), collapse = ", ")
      reasons <- c(
        reasons,
        paste0("- **オルソログ**: ", winner$Gene_ID, " (", ortho_w, ") vs ",
               loser$Gene_ID, " (", ortho_l, ")")
      )
    }
    if (winner$scores$literature > 0 && loser$scores$literature > 0) {
      reasons <- c(
        reasons,
        "- **文献**: 両遺伝子とも汎用的な PubMed ヒットが付いており、候補判断への寄与は限定的です"
      )
    }

    hdr <- paste0(
      "**", winner$Gene_ID, "** の方が候補としてより信頼できます（総合スコア ",
      sprintf("%.3f", winner$composite_score), " vs ",
      sprintf("%.3f", loser$composite_score), "）。\n\n",
      if (length(reasons)) "主な理由:\n" else "",
      paste(reasons, collapse = "\n")
    )
    return(hdr)
  }

  # English
  if (winner$scores$annotation > loser$scores$annotation + 0.05) {
    reasons <- c(reasons, paste0(
      "- **Annotation**: ", winner$Gene_ID, " matches the phenotype query better."
    ))
  }
  if (!is.na(winner$negLog10P)) {
    reasons <- c(reasons, paste0(
      "- **GWAS**: negLog10P ", sprintf("%.2f", winner$negLog10P),
      " vs ", sprintf("%.2f", loser$negLog10P)
    ))
  }
  if (winner$scores$ortholog > loser$scores$ortholog + 0.01) {
    reasons <- c(reasons, paste0(
      "- **Ortholog**: ", winner$Gene_ID, " has a stronger ortholog score."
    ))
  }
  paste0(
    "**", winner$Gene_ID, "** is the more reliable candidate (composite ",
    sprintf("%.3f", winner$composite_score), " vs ",
    sprintf("%.3f", loser$composite_score), ").\n\n",
    if (length(reasons)) paste(paste(reasons, collapse = "\n")) else ""
  )
}

.answer_ranking_rationale <- function(query, evidence_bundle, language) {
  top <- evidence_bundle[[1L]]
  lines <- character()
  if (language == "ja") {
    lines <- c(
      lines,
      paste0("**", top$Gene_ID, "** が1位なのは、総合スコア（",
             sprintf("%.3f", top$composite_score), "）が最も高いためです。")
    )
    if (top$scores$annotation > 0) {
      ann <- paste(.evidence_snippet_text(top, "annotation"), collapse = " ")
      lines <- c(lines, paste0("- アノテーション: ", ann))
    }
    if (!is.na(top$negLog10P)) {
      lines <- c(lines, paste0("- GWAS 有意性: negLog10P = ", sprintf("%.2f", top$negLog10P)))
    }
    return(paste(lines, collapse = "\n"))
  }
  paste0(
    "**", top$Gene_ID, "** ranks first with composite score ",
    sprintf("%.3f", top$composite_score), "."
  )
}

.explain_phenotype_template <- function(query, evidence_bundle, language) {
  hdr <- if (language == "ja") {
    paste0("## 表現型探索結果: ", query$trait_text, "\n\n")
  } else {
    paste0("## Phenotype exploration: ", query$trait_text, "\n\n")
  }

  lines <- character()
  for (item in evidence_bundle) {
    title <- if (language == "ja") {
      sprintf("### %d. %s (総合スコア %.3f)", item$rank, item$Gene_ID, item$composite_score)
    } else {
      sprintf("### %d. %s (composite %.3f)", item$rank, item$Gene_ID, item$composite_score)
    }
    lines <- c(lines, title)

    ev <- item$evidence
    for (src in names(ev)) {
      block <- ev[[src]]
      if (is.null(block) || !length(block$snippets)) {
        next
      }
      snip <- block$snippets
      if (length(snip) > 3L) {
        snip <- snip[seq_len(3L)]
      }
      lines <- c(
        lines,
        sprintf("- **%s** (score %.2f): %s", src, block$score %||% 0, paste(snip, collapse = "; "))
      )
    }
    lines <- c(lines, "")
  }

  paste0(hdr, paste(lines, collapse = "\n"))
}

.llm_model <- function(model) {
  settings <- .ai_config_llm_settings(model = model, base_url = NULL, timeout = NULL)
  settings$model
}

.llm_base_url <- function(base_url) {
  settings <- .ai_config_llm_settings(model = NULL, base_url = base_url, timeout = NULL)
  settings$base_url
}

.llm_http_get <- function(path, base_url, timeout = 10) {
  url <- paste0(base_url, path)
  if (requireNamespace("httr2", quietly = TRUE)) {
    req <- httr2::request(url)
    req <- httr2::req_timeout(req, timeout)
    resp <- tryCatch(httr2::req_perform(req), error = function(e) NULL)
    if (is.null(resp)) {
      return(NULL)
    }
    return(httr2::resp_body_json(resp))
  }

  .llm_http_get_base(url, timeout)
}

.llm_http_post <- function(path, body, base_url, timeout = 120) {
  url <- paste0(base_url, path)
  payload <- jsonlite::toJSON(body, auto_unbox = TRUE)

  if (requireNamespace("httr2", quietly = TRUE)) {
    req <- httr2::request(url)
    req <- httr2::req_timeout(req, timeout)
    req <- httr2::req_body_raw(req, payload, "application/json")
    resp <- httr2::req_perform(req)
    return(httr2::resp_body_json(resp))
  }

  .llm_http_post_base(url, payload, timeout)
}

.llm_http_get_base <- function(url, timeout) {
  con <- utils::url(url, open = "rt")
  on.exit(close(con), add = TRUE)
  if (timeout > 0) {
    setTimeLimit(elapsed = timeout, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE), add = TRUE)
  }
  txt <- readLines(con, warn = FALSE)
  if (!length(txt)) {
    return(NULL)
  }
  jsonlite::fromJSON(paste(txt, collapse = "\n"), simplifyVector = FALSE)
}

.llm_http_post_base <- function(url, payload, timeout) {
  if (requireNamespace("curl", quietly = TRUE)) {
    h <- curl::new_handle()
    curl::handle_setheaders(h, "Content-Type" = "application/json")
    curl::handle_setopt(h, postfields = payload)
    if (timeout > 0) {
      curl::handle_setopt(h, timeout = timeout)
    }
    txt <- curl::curl_fetch_memory(url, handle = h)$content
    return(jsonlite::fromJSON(rawToChar(txt), simplifyVector = FALSE))
  }

  stop(
    "HTTP requests require 'httr2' or 'curl'. Install with install.packages(\"httr2\").",
    call. = FALSE
  )
}
