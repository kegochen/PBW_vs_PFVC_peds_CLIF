# Validated site configuration; loading does not create output directories.
load_config <- function(repo_root, config_path = file.path("config", "config.json")) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required. Install the project dependencies before loading configuration.")
  }
  
  validate_text <- function(value, label) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
      stop(label, " must be one nonempty character value.")
    }
    value
  }
  
  is_absolute_path <- function(path) {
    grepl("^([A-Za-z]:[/\\\\]|[/\\\\]|~)", path)
  }
  
  repo_root <- validate_text(repo_root, "repo_root")
  repo_root <- normalizePath(path.expand(repo_root), winslash = "/", mustWork = TRUE)
  
  if (!dir.exists(repo_root)) {
    stop("repo_root must be a directory.")
  }
  
  resolve_path <- function(path, label, must_exist = FALSE) {
    path <- path.expand(validate_text(path, label))
    
    if (!is_absolute_path(path)) {
      path <- file.path(repo_root, path)
    }
    
    if (must_exist && !file.exists(path)) {
      stop(label, " does not exist: ", path)
    }
    
    normalizePath(path, winslash = "/", mustWork = must_exist)
  }
  
  config_path <- resolve_path(config_path, "config_path", must_exist = TRUE)
  
  if (dir.exists(config_path)) {
    stop("config_path must identify a JSON file.")
  }
  
  config <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
  
  required_keys <- c("site_name", "tables_path", "file_type", "output_root")
  optional_keys <- "supplemental_height"
  documentation_keys <- c("_comment", "_supplemental_height_example")
  allowed_keys <- c(required_keys, optional_keys, documentation_keys)
  
  if (!is.list(config) || is.null(names(config)) || anyDuplicated(names(config))) {
    stop("Configuration must be a JSON object with unique named keys.")
  }
  
  missing_keys <- setdiff(required_keys, names(config))
  unknown_keys <- setdiff(names(config), allowed_keys)
  
  if (length(missing_keys) > 0L) {
    stop("Missing configuration keys: ", paste(missing_keys, collapse = ", "))
  }
  
  if (length(unknown_keys) > 0L) {
    stop("Unknown configuration keys: ", paste(unknown_keys, collapse = ", "))
  }
  
  # Documentation keys are never interpreted as configuration inputs.
  config <- config[!(names(config) %in% documentation_keys)]
  
  for (key in required_keys) {
    config[[key]] <- validate_text(config[[key]], key)
  }
  
  if (
    !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", config$site_name) ||
    config$site_name %in% c("YOUR_SITE", "Your_Site_Name")
  ) {
    stop("site_name must be a non-placeholder label using letters, digits, underscores, or hyphens.")
  }
  
  if (!config$file_type %in% c("parquet", "csv", "fst")) {
    stop("file_type must be exactly one of parquet, csv, or fst.")
  }
  
  config$tables_path <- resolve_path(config$tables_path, "tables_path", must_exist = TRUE)
  
  if (!dir.exists(config$tables_path)) {
    stop("tables_path must be a directory.")
  }
  
  config$output_root <- resolve_path(config$output_root, "output_root")
  
  if (file.exists(config$output_root) && !dir.exists(config$output_root)) {
    stop("output_root exists but is not a directory.")
  }
  
  comparison_path <- function(path) {
    path <- sub("/+$", "", gsub("\\\\", "/", path))
    
    if (.Platform$OS.type == "windows") {
      tolower(path)
    } else {
      path
    }
  }
  
  tables_comparison <- comparison_path(config$tables_path)
  output_comparison <- comparison_path(config$output_root)
  
  if (
    tables_comparison == output_comparison ||
    startsWith(output_comparison, paste0(tables_comparison, "/")) ||
    startsWith(tables_comparison, paste0(output_comparison, "/"))
  ) {
    stop("tables_path and output_root must be separate, non-nested directories.")
  }
  
  supplemental_height <- config[["supplemental_height"]]
  
  if (!is.null(supplemental_height)) {
    if (
      !is.list(supplemental_height) ||
      is.null(names(supplemental_height)) ||
      anyDuplicated(names(supplemental_height))
    ) {
      stop("supplemental_height must be null or a JSON object with unique named keys.")
    }
    
    supplemental_required_keys <- c("format", "path")
    supplemental_missing_keys <- setdiff(supplemental_required_keys, names(supplemental_height))
    supplemental_unknown_keys <- setdiff(names(supplemental_height), supplemental_required_keys)
    
    if (length(supplemental_missing_keys) > 0L) {
      stop(
        "Missing supplemental_height keys: ",
        paste(supplemental_missing_keys, collapse = ", ")
      )
    }
    
    if (length(supplemental_unknown_keys) > 0L) {
      stop(
        "Unknown supplemental_height keys: ",
        paste(supplemental_unknown_keys, collapse = ", ")
      )
    }
    
    supplemental_height$format <- validate_text(
      supplemental_height$format,
      "supplemental_height.format"
    )
    
    if (!supplemental_height$format %in% c("michigan_crosswalk", "height_csv")) {
      stop("supplemental_height.format must be exactly one of michigan_crosswalk or height_csv.")
    }
    
    supplemental_height$path <- resolve_path(
      supplemental_height$path,
      "supplemental_height.path",
      must_exist = FALSE
    )
    
    if (dir.exists(supplemental_height$path)) {
      stop("supplemental_height.path must identify a file, not a directory.")
    }
  }
  
  # Single-bracket assignment preserves an explicit NULL entry in the returned list.
  config["supplemental_height"] <- list(supplemental_height)
  
  config$repo_root <- repo_root
  config$config_path <- config_path
  config$site_output_dir <- file.path(config$output_root, paste0(config$site_name, "_output"))
  config$output_dir <- file.path(config$site_output_dir, "intermediate")
  config$diagnostic_dir <- file.path(config$site_output_dir, "diagnostics")
  config$cohort_diagnostic_dir <- file.path(config$diagnostic_dir, "cohort")
  config$table1_diagnostic_dir <- file.path(config$diagnostic_dir, "table1")
  config$model_diagnostic_dir <- file.path(config$diagnostic_dir, "models")
  config$final_dir <- file.path(config$site_output_dir, "final")
  config$table_dir <- file.path(config$final_dir, "tables")
  config$figure_dir <- file.path(config$final_dir, "figures")
  config$report_dir <- file.path(config$final_dir, "reports")
  config$optional_diagnostic_export_dir <- file.path(config$final_dir, "optional_diagnostics")
  config$utils_dir <- file.path(repo_root, "utils")
  config$thresholds_dir <- file.path(repo_root, "outlier-thresholds")
  config$growth_reference_dir <- file.path(repo_root, "reference-data", "growth")
  
  required_directories <- c(config$utils_dir, config$thresholds_dir, config$growth_reference_dir)
  
  if (any(!dir.exists(required_directories))) {
    stop(
      "Missing repository directories: ",
      paste(required_directories[!dir.exists(required_directories)], collapse = ", ")
    )
  }
  
  message("Configuration validated for site: ", config$site_name)
  config
}