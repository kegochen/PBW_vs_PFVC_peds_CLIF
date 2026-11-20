# =============================================================================
# Script 00: Local single-site workflow orchestration
# =============================================================================
# Sourcing this file defines functions only. It does not execute the pipeline.
# Invoke only after execution is explicitly authorized.
#
# Analysis scripts run sequentially in separate R processes using the project
# library and base R library, with user and system profiles disabled.
# Reporting is a separate explicit invocation against saved diagnostics.
#
# Existing script outputs may be overwritten by the analysis scripts.
# This runner does not install packages, delete outputs, pool sites, or use Git.

pbwpfvc_runner_context <- function(
    repo_root = Sys.getenv("PBWPFVC_REPO_ROOT", unset = ""),
    config_path = Sys.getenv("PBWPFVC_CONFIG_PATH", unset = "")
) {
  if (
    !is.character(repo_root) || length(repo_root) != 1L ||
    is.na(repo_root)
  ) {
    stop("repo_root must be one character value.")
  }
  
  if (!nzchar(repo_root)) {
    repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    
    repeat {
      markers_present <- file.exists(file.path(repo_root, "utils", "config.R")) &&
        dir.exists(file.path(repo_root, "code")) &&
        dir.exists(file.path(repo_root, "config"))
      
      if (markers_present) break
      
      parent <- dirname(repo_root)
      
      if (identical(parent, repo_root)) {
        stop("Repository not found. Supply repo_root or PBWPFVC_REPO_ROOT.")
      }
      
      repo_root <- parent
    }
  } else {
    repo_root <- normalizePath(
      path.expand(repo_root), winslash = "/", mustWork = TRUE
    )
  }
  
  if (
    !is.character(config_path) || length(config_path) != 1L ||
    is.na(config_path)
  ) {
    stop("config_path must be one character value.")
  }
  
  if (!nzchar(config_path)) {
    config_path <- file.path("config", "config.json")
  }
  
  loader <- file.path(repo_root, "utils", "config.R")
  
  if (!file.exists(loader) || dir.exists(loader)) {
    stop("Configuration loader not found: ", loader)
  }
  
  loader_environment <- new.env(parent = baseenv())
  sys.source(loader, envir = loader_environment)
  
  config <- loader_environment$load_config(
    repo_root = repo_root,
    config_path = config_path
  )
  
  if (!requireNamespace("renv", quietly = TRUE)) {
    stop("renv is required. No installation was attempted.")
  }
  
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop("callr is required. No installation was attempted.")
  }
  
  project_library <- renv::paths$library(project = repo_root)
  
  if (!dir.exists(project_library)) {
    stop("Project library not found: ", project_library)
  }
  
  list(
    config = config,
    libpath = unique(c(
      normalizePath(project_library, winslash = "/", mustWork = TRUE),
      normalizePath(.Library, winslash = "/", mustWork = TRUE)
    )),
    environment = c(
      PBWPFVC_REPO_ROOT = config$repo_root,
      PBWPFVC_CONFIG_PATH = config$config_path,
      R_LIBS_SITE = ""
    )
  )
}

run_pbwpfvc_analysis <- function(
    repo_root = Sys.getenv("PBWPFVC_REPO_ROOT", unset = ""),
    config_path = Sys.getenv("PBWPFVC_CONFIG_PATH", unset = "")
) {
  context <- pbwpfvc_runner_context(repo_root, config_path)
  
  script_names <- c(
    "01_cohort_identification_v3.R",
    "02_quality_checks_v3.R",
    "03_variable_derivation_v3.R",
    "04_analysis_v3.R",
    "05_normalization_analysis_v3.R"
  )
  
  script_paths <- file.path(context$config$repo_root, "code", script_names)
  
  missing <- script_paths[
    !file.exists(script_paths) | dir.exists(script_paths)
  ]
  
  if (length(missing) > 0L) {
    stop(
      "Required analysis scripts unavailable: ",
      paste(missing, collapse = "; ")
    )
  }
  
  message("Analysis starting: ", context$config$site_name)

  warning_root <- context$config$model_diagnostic_dir
  if (!dir.exists(warning_root)) {
    if (!dir.create(warning_root, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create analysis warning destination: ", warning_root)
    }
  }
  warning_run_dir <- tempfile(
    pattern = "analysis_warning_run_", tmpdir = warning_root
  )
  if (!dir.create(warning_run_dir, recursive = FALSE)) {
    stop("Could not create retained analysis warning-run directory.")
  }
  warning_run_dir <- normalizePath(
    warning_run_dir, winslash = "/", mustWork = TRUE
  )
  message("Retained analysis warnings: ", warning_run_dir)

  completed_scripts <- character()
  
  for (index in seq_along(script_paths)) {
    message("Starting: ", script_names[[index]])
    
    result <- tryCatch(
      callr::r(
        func = function(script_path, warning_path) {
          # Capture state stays local when the script clears its global workspace.
          warning_number <- 0L
          header <- data.frame(
            warning_number = integer(), timestamp_utc = character(),
            script = character(), condition_class = character(),
            message = character(), call = character(),
            stringsAsFactors = FALSE
          )
          connection <- file(warning_path, open = "wt", encoding = "UTF-8")
          tryCatch(
            write.table(
              header, connection, sep = ",", row.names = FALSE,
              col.names = TRUE, quote = TRUE
            ),
            finally = close(connection)
          )

          retain_warning <- function(condition) {
            warning_number <<- warning_number + 1L
            warning_call <- conditionCall(condition)
            row <- data.frame(
              warning_number = warning_number,
              timestamp_utc = format(
                Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"
              ),
              script = basename(script_path),
              condition_class = paste(class(condition), collapse = "/"),
              message = conditionMessage(condition),
              call = if (is.null(warning_call)) "" else
                paste(deparse(warning_call), collapse = " "),
              stringsAsFactors = FALSE
            )
            tryCatch(
              {
                connection <- file(
                  warning_path, open = "at", encoding = "UTF-8"
                )
                tryCatch(
                  write.table(
                    row, connection, sep = ",", row.names = FALSE,
                    col.names = FALSE, quote = TRUE
                  ),
                  finally = close(connection)
                )
              },
              error = function(error) {
                stop(
                  "Warning retention failed: ", warning_path,
                  ". Original warning: ", conditionMessage(condition),
                  ". Retention error: ", conditionMessage(error),
                  call. = FALSE
                )
              }
            )
            # Do not invoke muffleWarning: existing warning behavior continues.
            invisible(NULL)
          }

          withCallingHandlers(
            base::sys.source(script_path, envir = globalenv()),
            warning = retain_warning
          )
          TRUE
        },
        args = list(
          script_path = script_paths[[index]],
          warning_path = file.path(
            warning_run_dir,
            paste0(tools::file_path_sans_ext(script_names[[index]]), "_warnings.csv")
          )
        ),
        libpath = context$libpath,
        wd = context$config$repo_root,
        system_profile = FALSE,
        user_profile = FALSE,
        env = context$environment,
        show = TRUE
      ),
      error = function(error) {
        stop(
          "Analysis stopped at ", script_names[[index]], ". ",
          "Previously returned scripts: ",
          if (length(completed_scripts)) {
            paste(completed_scripts, collapse = ", ")
          } else {
            "none"
          },
          ". Existing or partial outputs were retained. ",
          "No later scripts or reports were started. Error: ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )
    
    if (!identical(result, TRUE)) {
      stop("Unexpected subprocess result from: ", script_names[[index]])
    }
    
    completed_scripts <- c(completed_scripts, script_names[[index]])
    invisible(gc())
    message("Completed: ", script_names[[index]])
  }
  
  message(
    "Analysis complete: ", length(completed_scripts),
    " scripts returned successfully."
  )
  
  invisible(list(
    site = context$config$site_name,
    completed_scripts = completed_scripts,
    reporting_started = FALSE
  ))
}

run_pbwpfvc_reports <- function(
    repo_root = Sys.getenv("PBWPFVC_REPO_ROOT", unset = ""),
    config_path = Sys.getenv("PBWPFVC_CONFIG_PATH", unset = "")
) {
  context <- pbwpfvc_runner_context(repo_root, config_path)
  
  utility <- file.path(
    context$config$utils_dir, "render_diagnostic_reports.R"
  )
  
  if (!file.exists(utility) || dir.exists(utility)) {
    stop("Reporting utility not found: ", utility)
  }
  
  # Render saved diagnostics and update reporting status without rerunning analysis.
  callr::r(
    func = function(utility_path, config) {
      reporting_environment <- new.env(parent = globalenv())
      base::sys.source(utility_path, envir = reporting_environment)
      reporting_environment$render_diagnostic_reports(config)
    },
    args = list(utility_path = utility, config = context$config),
    libpath = context$libpath,
    wd = context$config$repo_root,
    system_profile = FALSE,
    user_profile = FALSE,
    env = context$environment,
    show = TRUE
  )
}
