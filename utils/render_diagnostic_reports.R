# =============================================================================
# Separate diagnostic PDF reporting
# =============================================================================
# Defining this function does not render reports.
# Invoke only when report execution is explicitly authorized.
#
# Reports execute QMD display code against saved diagnostics; no analysis scripts
# are sourced and no models are refitted.
#
# Each invocation retains a unique run directory and backups. No files are
# explicitly deleted. Final PDFs are published only after all renders validate.
# The saved diagnostics are updated after successful rendering and validation;
# report content therefore describes the pre-render diagnostic snapshot.

render_diagnostic_reports <- function(config) {
  reporting_operation <- "Validate reporting configuration and prerequisites"
  reporting_run_directory <- NULL
  reporting_completed <- FALSE
  
  on.exit({
    if (!reporting_completed) {
      message(
        "Reporting did not reach complete success. Last operation attempted: ",
        reporting_operation
      )
      if (!is.null(reporting_run_directory)) {
        message("Retained reporting run directory: ", reporting_run_directory)
      }
      message(
        "Some files or diagnostic status fields may already have changed. ",
        "The last attempted operation is not proof of success or failure of its write. ",
        "Standalone normalization diagnostics and model QC may disagree. ",
        "No automatic retry, rollback, or deletion was performed. ",
        "Inspect the original error, affected files, and retained backups before retrying."
      )
    }
  }, add = TRUE)
  
  required_config <- c(
    "repo_root", "site_name", "utils_dir", "report_dir",
    "cohort_diagnostic_dir", "table1_diagnostic_dir", "model_diagnostic_dir"
  )
  
  if (!is.list(config)) {
    stop("config must be the validated list returned by load_config().")
  }
  
  for (name in required_config) {
    value <- config[[name]]
    if (
      !is.character(value) || length(value) != 1L ||
      is.na(value) || !nzchar(value)
    ) {
      stop("Missing or invalid configuration field: ", name)
    }
  }
  
  if (!requireNamespace("quarto", quietly = TRUE)) {
    stop("The project quarto R package is required. No installation was attempted.")
  }
  
  quarto_executable <- quarto::quarto_path()
  
  if (is.null(quarto_executable) || !file.exists(quarto_executable)) {
    stop("Quarto CLI is unavailable. No installation was attempted.")
  }
  
  if (!any(nzchar(Sys.which(c("pdflatex", "xelatex", "lualatex"))))) {
    stop("No TeX executable was found. No installation was attempted.")
  }
  
  # This routing policy requires standalone documents, not inherited projects.
  directory <- normalizePath(config$utils_dir, winslash = "/", mustWork = TRUE)
  
  repeat {
    candidates <- file.path(directory, c("_quarto.yml", "_quarto.yaml"))
    found <- candidates[file.exists(candidates)]
    
    if (length(found) > 0L) {
      stop(
        "Quarto project configuration now exists; review routing before rendering: ",
        paste(found, collapse = "; ")
      )
    }
    
    parent <- dirname(directory)
    if (identical(parent, directory)) break
    directory <- parent
  }
  
  require_file <- function(path) {
    if (!file.exists(path) || dir.exists(path)) {
      stop("Required file is unavailable: ", path)
    }
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  
  validate_checks <- function(checks, label) {
    if (
      !is.list(checks) || length(checks) == 0L ||
      is.null(names(checks)) || anyNA(names(checks)) ||
      any(!nzchar(names(checks))) || anyDuplicated(names(checks))
    ) {
      stop(label, " must contain nonempty, uniquely named checks.")
    }
    
    valid <- vapply(
      checks,
      function(value) {
        is.logical(value) && length(value) == 1L && !is.na(value)
      },
      logical(1)
    )
    
    if (!all(valid)) {
      stop(label, " contains malformed checks.")
    }
    
    passed <- vapply(checks, isTRUE, logical(1))
    
    if (!all(passed)) {
      stop(
        label, " contains failed checks: ",
        paste(names(checks)[!passed], collapse = ", ")
      )
    }
  }
  
  save_readback <- function(object, path) {
    reporting_operation <<- paste("Write RDS:", path)
    saveRDS(object, path, compress = "xz")
    
    reporting_operation <<- paste("Read back and compare complete RDS:", path)
    if (!identical(readRDS(path), object)) {
      stop("Complete RDS read-back failed: ", path)
    }
  }
  
  copy_verified <- function(from, to, overwrite = FALSE) {
    reporting_operation <<- paste("Copy file:", from, "->", to)
    
    if (!isTRUE(file.copy(from, to, overwrite = overwrite))) {
      stop("File copy failed: ", from, " -> ", to)
    }
    
    reporting_operation <<- paste("Verify copied file checksum:", to)
    
    if (!identical(
      unname(tools::md5sum(from)),
      unname(tools::md5sum(to))
    )) {
      stop("Copied file failed checksum comparison: ", to)
    }
  }
  
  validate_pdf <- function(path) {
    require_file(path)
    size <- file.info(path, extra_cols = FALSE)$size[[1]]
    
    if (!is.finite(size) || size <= 5L) {
      stop("PDF is empty or too small: ", path)
    }
    
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    signature <- readBin(connection, what = "raw", n = 5L)
    
    if (!identical(signature, charToRaw("%PDF-"))) {
      stop("Output does not have a PDF signature: ", path)
    }
    
    invisible(TRUE)
  }
  
  inspect_manifest <- function(manifest) {
    if (
      !is.data.frame(manifest) ||
      !all(c("output", "path") %in% names(manifest)) ||
      !is.character(manifest$path) || anyNA(manifest$path) ||
      any(!nzchar(manifest$path)) || anyDuplicated(manifest$path)
    ) {
      stop("Output manifest has an invalid structure or duplicate paths.")
    }
    
    manifest$exists <- file.exists(manifest$path) & !dir.exists(manifest$path)
    manifest$size_bytes <- vapply(
      manifest$path,
      function(path) {
        if (file.exists(path) && !dir.exists(path)) {
          file.info(path, extra_cols = FALSE)$size[[1]]
        } else {
          NA_real_
        }
      },
      numeric(1)
    )
    manifest$nonempty <- manifest$exists &
      is.finite(manifest$size_bytes) & manifest$size_bytes > 0
    manifest$extension <- tools::file_ext(manifest$path)
    
    if (!all(manifest$nonempty)) {
      stop("One or more manifest artifacts are missing or empty.")
    }
    
    manifest
  }
  
  validate_types <- function(manifest, expected_n) {
    if (
      nrow(manifest) != sum(expected_n) ||
      any(!manifest$extension %in% names(expected_n))
    ) {
      stop("Manifest file count or extensions do not match the contract.")
    }
    
    observed <- vapply(
      names(expected_n),
      function(extension) sum(manifest$extension == extension),
      integer(1)
    )
    
    if (!identical(unname(observed), unname(expected_n))) {
      stop("Manifest output-type counts do not match the contract.")
    }
    
    data.frame(
      extension = names(expected_n),
      expected_n = unname(expected_n),
      observed_n = unname(observed),
      count_matches = TRUE,
      stringsAsFactors = FALSE
    )
  }
  
  site <- config$site_name
  paths <- list(
    cohort = require_file(file.path(
      config$cohort_diagnostic_dir, "cohort_diagnostics.rds"
    )),
    table1 = require_file(file.path(
      config$table1_diagnostic_dir, paste0("table1_qc_", site, ".rds")
    )),
    model = require_file(file.path(
      config$model_diagnostic_dir, paste0("model_qc_", site, ".rds")
    )),
    normalization = require_file(file.path(
      config$model_diagnostic_dir, "normalization_analysis_diagnostics.rds"
    ))
  )
  
  table1 <- readRDS(paths$table1)
  model <- readRDS(paths$model)
  normalization <- readRDS(paths$normalization)
  
  validate_checks(table1$checks, "Table 1 checks")
  validate_checks(model$checks, "Model checks")
  validate_checks(normalization$checks, "Normalization checks")
  
  if (
    !identical(model$metadata$site, site) ||
    !identical(normalization$metadata$site, site) ||
    !identical(normalization$metadata$computational_status, "PASS")
  ) {
    stop("Diagnostic site or normalization computational status mismatch.")
  }
  
  analysis_manifest_before <- normalization$output_manifest
  analysis_type_counts_before <- normalization$output_type_counts
  analysis_manifest_checked <- inspect_manifest(analysis_manifest_before)
  validate_types(analysis_manifest_checked, c(csv = 14L, pdf = 5L))
  
  source_paths <- c(
    cohort = file.path(config$utils_dir, "cohort_diagnostic_report.qmd"),
    analysis = file.path(config$utils_dir, "analysis_diagnostic_report.qmd"),
    normalization = file.path(
      config$utils_dir, "normalization_diagnostic_report.qmd"
    )
  )
  source_paths <- vapply(source_paths, require_file, character(1))
  
  report_names <- c(
    cohort = paste0("cohort_diagnostic_report_", site, ".pdf"),
    analysis = paste0("analysis_diagnostic_report_", site, ".pdf"),
    normalization = paste0("normalization_diagnostic_report_", site, ".pdf")
  )
  
  final_paths <- file.path(config$report_dir, unname(report_names))
  names(final_paths) <- names(report_names)
  
  if (
    !identical(
      normalizePath(
        normalization$reporting$intended_output,
        winslash = "/", mustWork = FALSE
      ),
      normalizePath(
        final_paths[["normalization"]],
        winslash = "/", mustWork = FALSE
      )
    ) ||
    !identical(normalization$reporting$report_inclusive_expected_file_count, 20L)
  ) {
    stop("Saved normalization reporting destination or count differs from contract.")
  }
  
  parameters <- list(
    cohort = list(diagnostic_file = paths$cohort, site_name = site),
    analysis = list(
      table1_qc_file = paths$table1,
      model_qc_file = paths$model,
      site_name = site
    ),
    normalization = list(
      diagnostic_file = paths$normalization, site_name = site
    )
  )
  
  dir.create(config$report_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(config$report_dir)) stop("Could not create report directory.")
  
  run_dir <- tempfile(pattern = "report_run_", tmpdir = config$report_dir)
  if (!dir.create(run_dir, recursive = FALSE)) {
    stop("Could not create a unique reporting run directory.")
  }
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  reporting_run_directory <- run_dir
  
  backup_dir <- file.path(run_dir, "backups")
  if (!dir.create(backup_dir)) stop("Could not create reporting backup directory.")
  
  started_at <- Sys.time()
  input_hashes <- tools::md5sum(unlist(paths, use.names = FALSE))
  source_hashes <- tools::md5sum(source_paths)
  
  # Preserve the original diagnostic artifacts before any reporting update.
  for (name in names(paths)) {
    copy_verified(
      paths[[name]],
      file.path(backup_dir, paste0(name, "_diagnostics_before_reporting.rds"))
    )
  }
  
  staged_paths <- setNames(character(length(report_names)), names(report_names))
  
  for (name in names(report_names)) {
    destination <- file.path(run_dir, name)
    if (!dir.create(destination)) stop("Could not create report staging directory.")
    
    reporting_operation <- paste("Render PDF:", source_paths[[name]])
    quarto::quarto_render(
      input = source_paths[[name]],
      output_format = "pdf",
      output_file = report_names[[name]],
      execute = TRUE,
      execute_params = parameters[[name]],
      execute_dir = config$repo_root,
      cache = FALSE,
      use_freezer = FALSE,
      quarto_args = c(
        "--output-dir",
        normalizePath(destination, winslash = "/", mustWork = TRUE),
        "--no-clean"
      ),
      quiet = FALSE,
      as_job = FALSE
    )
    
    staged_paths[[name]] <- file.path(destination, report_names[[name]])
    validate_pdf(staged_paths[[name]])
  }
  
  # Reject changes to saved inputs or QMD sources during rendering.
  if (
    !identical(tools::md5sum(unlist(paths, use.names = FALSE)), input_hashes) ||
    !identical(tools::md5sum(source_paths), source_hashes)
  ) {
    stop("Inputs or report sources changed during rendering; reports were not published.")
  }
  
  for (name in names(final_paths)) {
    if (dir.exists(final_paths[[name]])) {
      stop("Final report destination is a directory: ", final_paths[[name]])
    }
    if (file.exists(final_paths[[name]])) {
      copy_verified(
        final_paths[[name]],
        file.path(backup_dir, paste0("previous_", report_names[[name]]))
      )
    }
  }
  
  for (name in names(final_paths)) {
    copy_verified(staged_paths[[name]], final_paths[[name]], overwrite = TRUE)
    validate_pdf(final_paths[[name]])
  }
  
  full_manifest <- rbind(
    as.data.frame(analysis_manifest_checked),
    as.data.frame(analysis_manifest_checked[1L, , drop = FALSE])
  )
  report_row <- nrow(full_manifest)
  full_manifest$output[report_row] <- "Normalization diagnostic report"
  full_manifest$path[report_row] <- final_paths[["normalization"]]
  full_manifest <- inspect_manifest(full_manifest)
  full_type_counts <- validate_types(full_manifest, c(csv = 14L, pdf = 6L))
  
  # Re-read the latest model artifact before replacing only its Script05 section.
  latest_model <- readRDS(paths$model)
  section <- "normalization_analysis_script05"
  existing_section <- latest_model[[section]]
  
  if (
    !identical(latest_model$metadata$site, site) ||
    (!is.null(existing_section) && (
      !is.list(existing_section) ||
      !identical(existing_section$metadata$site, site) ||
      !identical(
        existing_section$metadata$script,
        normalization$metadata$script
      )
    ))
  ) {
    stop("Latest model QC or its existing Script05 section has an identity mismatch.")
  }
  
  preserved_sections <- setdiff(names(latest_model), section)
  preserved_object <- latest_model[preserved_sections]
  
  completed_at <- Sys.time()
  normalization$report_inclusive_output_manifest <- full_manifest
  normalization$report_inclusive_output_type_counts <- full_type_counts
  normalization$reporting$status <- "COMPLETE"
  normalization$reporting$started_at <- started_at
  normalization$reporting$completed_at <- completed_at
  normalization$reporting$run_directory <- run_dir
  normalization$reporting$rendered_by <- "utils/render_diagnostic_reports.R"
  normalization$reporting$report_rendered_by_script05 <- FALSE
  normalization$reporting$validation <- list(
    pdf_signature_valid = TRUE,
    published_copy_matches_rendered_pdf = TRUE,
    report_inclusive_manifest_valid = TRUE,
    input_hashes_before_render = input_hashes,
    qmd_hashes_before_render = source_hashes,
    note = paste(
      "Checks establish file generation and manifest integrity.",
      "They do not establish scientific validity or disclosure approval.",
      "The PDF displays the pre-render diagnostic snapshot."
    )
  )
  normalization$reporting$note <- paste(
    "Separate reporting completed and the report-inclusive manifest validated.",
    "The original Script05 analysis manifest and analysis statuses are preserved.",
    "Retained run files and backups are not additional analysis-manifest artifacts."
  )
  normalization$metadata$reporting_status <- "COMPLETE"
  normalization$metadata$report_input_stage <-
    "Saved analysis diagnostics rendered; report-inclusive manifest validated"
  normalization$metadata$updated_at <- completed_at
  normalization$final_summary$reporting_status <- "COMPLETE"
  normalization$final_summary$scope <-
    "Script 05 analysis artifacts and diagnostics; separate reporting completed"
  
  if (
    !identical(normalization$output_manifest, analysis_manifest_before) ||
    !identical(normalization$output_type_counts, analysis_type_counts_before)
  ) {
    stop("The original analysis manifest was altered.")
  }
  
  save_readback(normalization, paths$normalization)
  latest_model[[section]] <- readRDS(paths$normalization)
  
  if (!identical(latest_model[preserved_sections], preserved_object)) {
    stop("Unrelated model QC sections were altered before synchronization.")
  }
  
  save_readback(latest_model, paths$model)
  model_readback <- readRDS(paths$model)
  
  if (
    !identical(model_readback[[section]], readRDS(paths$normalization)) ||
    !identical(model_readback[preserved_sections], preserved_object)
  ) {
    stop("Normalization synchronization or section preservation failed.")
  }
  
  result <- list(
    status = "COMPLETE",
    site = site,
    started_at = started_at,
    completed_at = completed_at,
    reports = final_paths,
    retained_run_directory = run_dir,
    normalization_diagnostics = paths$normalization,
    synchronized_model_diagnostics = paths$model,
    report_inclusive_output_manifest = full_manifest
  )
  
  save_readback(result, file.path(run_dir, "reporting_result.rds"))
  reporting_completed <- TRUE
  
  message("Diagnostic reporting complete for site: ", site)
  message("Validated PDFs: ", paste(final_paths, collapse = "; "))
  message("Retained run files and backups: ", run_dir)
  message("Normalization reporting status saved and synchronized into model QC.")
  
  invisible(result)
}