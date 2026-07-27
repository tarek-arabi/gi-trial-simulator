#' Parameter packs
#'
#' A parameter pack is a YAML file describing one clinical scenario: its arms,
#' its endpoints with control and treatment event rates, the published sources
#' each value came from, and the design defaults a trialist would start from.
#' Packs are data, not code, so they can be versioned, cited and audited
#' independently of the engine that consumes them.
#'
#' Packs are found on a search path: the directories in
#' `getOption("gitrialsim.pack_paths")`, then those in the
#' `GITRIALSIM_PACK_PATH` environment variable (separated by the platform's path
#' separator, `.Platform$path.sep`),
#' then the packs shipped inside the installed package. Earlier directories win,
#' so a pack held outside this repository can extend or override a shipped one
#' without modifying the package.
#'
#' @name parameter-packs
NULL

#' Directories searched for parameter packs
#'
#' @return Character vector of existing directories, in search order.
#' @seealso [parameter-packs]
#' @examples
#' pack_search_path()
#' @export
pack_search_path <- function() {
  from_option <- getOption("gitrialsim.pack_paths", character())
  from_env <- Sys.getenv("GITRIALSIM_PACK_PATH", unset = "")
  from_env <- if (nzchar(from_env)) {
    strsplit(from_env, .Platform$path.sep, fixed = TRUE)[[1]]
  } else {
    character()
  }
  shipped <- system.file("parameters", package = "gitrialsim")
  paths <- c(as.character(from_option), from_env, shipped)
  paths <- paths[nzchar(paths)]
  unique(paths[dir.exists(paths)])
}

#' List available parameter packs
#'
#' @return A data frame with one row per pack: `id`, `title`, `version`,
#'   `provenance`, `restricted` and the `path` it was found at. Where two
#'   directories supply the same `id`, the one earlier on the search path wins.
#' @seealso [pack_search_path()], [load_pack()]
#' @examples
#' list_packs()
#' @export
list_packs <- function() {
  files <- unlist(lapply(
    pack_search_path(),
    list.files,
    pattern = "\\.ya?ml$", full.names = TRUE
  ))
  if (length(files) == 0L) {
    return(data.frame(
      id = character(), title = character(), version = character(),
      provenance = character(), restricted = logical(), path = character(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(files, function(f) {
    p <- try(yaml::read_yaml(f), silent = TRUE)
    if (inherits(p, "try-error") || is.null(p$id)) {
      return(NULL)
    }
    data.frame(
      id = p$id,
      title = p$title %||% NA_character_,
      version = as.character(p$version %||% NA_character_),
      provenance = p$provenance %||% NA_character_,
      restricted = isTRUE(p$restricted),
      path = f,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  out <- out[!duplicated(out$id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Load a parameter pack
#'
#' @param id Pack identifier (as reported by [list_packs()]), or a path to a
#'   YAML file.
#' @return An object of class `gi_pack`.
#' @seealso [list_packs()], [validate_pack()], [scenario()]
#' @examples
#' pack <- load_pack("ercp_acute_cholangitis")
#' pack$title
#' @export
load_pack <- function(id) {
  stopifnot(is.character(id), length(id) == 1L, nzchar(id))
  path <- if (file.exists(id) && !dir.exists(id)) {
    id
  } else {
    available <- list_packs()
    hit <- available$path[available$id == id]
    if (length(hit) == 0L) {
      stop(
        "No parameter pack with id '", id, "'. Available: ",
        if (nrow(available)) paste(available$id, collapse = ", ") else "none",
        ". Searched: ", paste(pack_search_path(), collapse = ", "),
        call. = FALSE
      )
    }
    hit[[1L]]
  }
  pack <- yaml::read_yaml(path)
  pack$path <- path
  class(pack) <- c("gi_pack", "list")
  validate_pack(pack)
}

#' Validate a parameter pack
#'
#' Checks the structural and numerical invariants a pack must satisfy before
#' any design is built on it: required fields present, event rates strictly
#' inside (0, 1), every endpoint citing a source that the pack declares, and no
#' endpoint whose direction contradicts its own rates.
#'
#' @param pack A `gi_pack`, as returned by [load_pack()].
#' @return The pack, invisibly, if valid. Otherwise an error describing every
#'   problem found.
#' @examples
#' validate_pack(load_pack("hrs_terlipressin"))
#' @export
validate_pack <- function(pack) {
  stopifnot(inherits(pack, "gi_pack"))
  problems <- character()

  for (field in c("id", "title", "version", "provenance", "endpoints", "sources")) {
    if (is.null(pack[[field]])) {
      problems <- c(problems, paste0("missing required field '", field, "'"))
    }
  }
  if (length(problems)) {
    stop("Invalid parameter pack:\n  - ", paste(problems, collapse = "\n  - "), call. = FALSE)
  }

  if (!pack$provenance %in% c("published_literature", "rwd_aggregate", "proprietary")) {
    problems <- c(problems, paste0(
      "provenance '", pack$provenance,
      "' is not one of published_literature, rwd_aggregate, proprietary"
    ))
  }

  for (key in names(pack$endpoints)) {
    ep <- pack$endpoints[[key]]
    where <- paste0("endpoint '", key, "'")

    if (!identical(ep$type, "binary")) {
      problems <- c(problems, paste0(where, ": only binary endpoints are supported"))
      next
    }
    for (rate_field in c("control_rate", "treatment_rate")) {
      r <- ep[[rate_field]]
      if (is.null(r) || !is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0 || r >= 1) {
        problems <- c(problems, paste0(
          where, ": ", rate_field, " must be a single number strictly between 0 and 1"
        ))
      }
    }
    if (is.null(ep$direction) ||
      !ep$direction %in% c("lower_is_better", "higher_is_better")) {
      problems <- c(problems, paste0(
        where, ": direction must be lower_is_better or higher_is_better"
      ))
    }
    if (is.null(ep$source) || !ep$source %in% names(pack$sources)) {
      problems <- c(problems, paste0(
        where, ": source '", ep$source %||% "<none>", "' is not declared in sources"
      ))
    }
  }

  for (key in names(pack$priors %||% list())) {
    pr <- pack$priors[[key]]
    where <- paste0("prior '", key, "'")
    if (is.null(pr$source) || !pr$source %in% names(pack$sources)) {
      problems <- c(problems, paste0(where, ": source is not declared in sources"))
    }
    if (!is.null(pr$ci_lower) && !is.null(pr$ci_upper) &&
      pr$ci_lower >= pr$ci_upper) {
      problems <- c(problems, paste0(where, ": ci_lower must be below ci_upper"))
    }
  }

  if (length(problems)) {
    stop(
      "Invalid parameter pack '", pack$id, "':\n  - ",
      paste(problems, collapse = "\n  - "),
      call. = FALSE
    )
  }
  invisible(pack)
}

#' Build a scenario from a pack endpoint
#'
#' A scenario is the minimal object every design function in this package
#' consumes: two event rates, the direction of benefit, and enough provenance
#' to reconstruct where the rates came from.
#'
#' @param pack A `gi_pack` or a pack id.
#' @param endpoint Endpoint key. Defaults to the endpoint whose `role` is
#'   `primary`.
#' @param control_rate,treatment_rate Optional overrides, for sensitivity
#'   analysis across values the published sources do not supply.
#' @return An object of class `gi_scenario`.
#' @examples
#' scenario("ercp_acute_cholangitis")
#' scenario("ercp_acute_cholangitis", treatment_rate = 0.05)
#' @export
scenario <- function(pack, endpoint = NULL, control_rate = NULL, treatment_rate = NULL) {
  if (is.character(pack)) pack <- load_pack(pack)
  stopifnot(inherits(pack, "gi_pack"))

  if (is.null(endpoint)) {
    roles <- vapply(pack$endpoints, function(e) e$role %||% "", character(1))
    primary <- names(pack$endpoints)[roles == "primary"]
    if (length(primary) != 1L) {
      stop(
        "Pack '", pack$id, "' does not have exactly one primary endpoint; ",
        "name one of: ", paste(names(pack$endpoints), collapse = ", "),
        call. = FALSE
      )
    }
    endpoint <- primary
  }
  ep <- pack$endpoints[[endpoint]]
  if (is.null(ep)) {
    stop(
      "Pack '", pack$id, "' has no endpoint '", endpoint, "'. Available: ",
      paste(names(pack$endpoints), collapse = ", "),
      call. = FALSE
    )
  }

  p_control <- control_rate %||% ep$control_rate
  p_treatment <- treatment_rate %||% ep$treatment_rate
  for (r in c(p_control, p_treatment)) {
    if (!is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0 || r >= 1) {
      stop("Event rates must be single numbers strictly between 0 and 1.", call. = FALSE)
    }
  }

  structure(
    list(
      pack_id = pack$id,
      pack_version = as.character(pack$version),
      endpoint = endpoint,
      label = ep$label %||% endpoint,
      control_rate = p_control,
      treatment_rate = p_treatment,
      direction = ep$direction,
      control_arm = pack$arms$control$label %||% "Control",
      treatment_arm = pack$arms$treatment$label %||% "Treatment",
      source = pack$sources[[ep$source]]$citation %||% NA_character_,
      overridden = !is.null(control_rate) || !is.null(treatment_rate),
      defaults = pack$design_defaults %||% list()
    ),
    class = c("gi_scenario", "list")
  )
}

#' @export
print.gi_pack <- function(x, ...) {
  cat("<gi_pack> ", x$id, " v", as.character(x$version), "\n", sep = "")
  cat(x$title, "\n", sep = "")
  cat("provenance: ", x$provenance,
    if (isTRUE(x$restricted)) "  [RESTRICTED - not redistributable]" else "",
    "\n",
    sep = ""
  )
  cat("endpoints: ", paste(names(x$endpoints), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
print.gi_scenario <- function(x, ...) {
  cat("<gi_scenario> ", x$pack_id, " / ", x$endpoint, "\n", sep = "")
  cat(x$label, "\n", sep = "")
  cat(sprintf(
    "  %-28s %.4f\n  %-28s %.4f\n",
    x$control_arm, x$control_rate, x$treatment_arm, x$treatment_rate
  ))
  cat("  direction: ", x$direction,
    if (isTRUE(x$overridden)) "   (rates overridden)" else "", "\n",
    sep = ""
  )
  invisible(x)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
