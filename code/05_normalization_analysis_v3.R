# =============================================================================
# Script 05: PBW vs PFVC normalization of physiologic injury metrics
# Discordance and exploratory prognostic utility
# Pediatric CLIF PBW vs PFVC study
# =============================================================================
# Run after Scripts 03 and 04. Script 05 consumes the cross-sectional analysis
# dataset and validated upstream diagnostics. Cross-site pooling is deferred.
#
# PART 1: Descriptive normalization discordance
# For positive values, the following identities hold:
#   log(Ers * PFVC) - log(Ers * PBW) = -log(PBW / PFVC)
#   log(MP / PFVC) - log(MP / PBW) = log(PBW / PFVC)
#
# These identities describe the numerical relationship between normalizers.
# PBW/PFVC is a unit-bearing conversion factor; it does not establish measured
# lung volume, correction of prediction error, or clinical superiority.
# Changing a size normalizer does not by itself resolve age-related recoil,
# pressure-measurement limitations, or other sources of measurement error.
#
# Descriptive analyses examine demographic patterns and changes between
# independently constructed, tie-preserving, cohort-relative metric terciles.
# These categories are not validated clinical injury thresholds.
#
# PART 2: Exploratory associations with 60-day mortality
# Candidate models compare mechanic-only, PBW-locked, PFVC-locked, and separate
# mechanic-plus-size specifications across seven alternative adjustment terms.
# Common complete-case samples are used within mechanic families.
# Further analyses evaluate demographic modification, age/PFVC interaction
# ladders, balanced single-modifier models, and encompassing size models.
#
# AIC and apparent in-sample discrimination are exploratory comparisons.
# Computational stability does not establish statistical reliability, causal
# identification, clinical validity, or a preferred normalization policy.
# Existing sparse-event interpretation and curve-display criteria are retained.
#
# Driving-pressure provenance and validity are checked without redefining the
# upstream outcome. Metric-specific positive/finite requirements determine the
# corresponding analysis frames; they do not change the upstream study cohort.
#
# Inputs:
#   Configured intermediate analysis_cross_sectional.parquet
#   Configured cohort diagnostic RDS and Script 04 model diagnostic RDS
#
# Outputs:
#   Existing CSV artifacts under the configured final/tables directory
#   Existing analysis PDFs under the configured final/figures directory
#   Detailed normalization diagnostics under local diagnostics/models
#
# Diagnostic-report rendering and report-inclusive manifest validation
# are separate from analysis execution.
# =============================================================================

rm(list = ls())

library(tidyverse)
library(arrow)
library(here)
library(broom)
library(splines)
library(patchwork)

# =============================================================================
# Load validated repository configuration
# =============================================================================
# Script 05 consumes the Script 03 analysis dataset and completed Script 04 QC.
# Cross-site pooling is a separate, deferred workflow.

repo_root_override <- Sys.getenv("PBWPFVC_REPO_ROOT", unset = "")

if (nzchar(repo_root_override)) {
  repo_root <- normalizePath(path.expand(repo_root_override), winslash = "/", mustWork = TRUE)
} else {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  
  repeat {
    repository_markers_present <- file.exists(file.path(repo_root, "utils", "config.R")) &&
      dir.exists(file.path(repo_root, "code")) &&
      dir.exists(file.path(repo_root, "config"))
    
    if (repository_markers_present) {
      break
    }
    
    parent_directory <- dirname(repo_root)
    
    if (identical(parent_directory, repo_root)) {
      stop(
        "Repository root could not be located. Run from the repository directory or set PBWPFVC_REPO_ROOT."
      )
    }
    
    repo_root <- parent_directory
  }
}

config_loader_path <- file.path(repo_root, "utils", "config.R")

if (!file.exists(config_loader_path)) {
  stop("Repository configuration loader not found: ", config_loader_path)
}

source(config_loader_path)

config_path_override <- Sys.getenv("PBWPFVC_CONFIG_PATH", unset = "")
config_path <- if (nzchar(config_path_override)) config_path_override else file.path("config", "config.json")

config <- load_config(repo_root = repo_root, config_path = config_path)

site_name <- config$site_name
output_dir <- config$output_dir
final_dir <- config$final_dir
table_dir <- config$table_dir
figure_dir <- config$figure_dir
report_dir <- config$report_dir
model_diagnostic_dir <- config$model_diagnostic_dir

required_script05_input_directories <- c(
  output_dir,
  config$cohort_diagnostic_dir,
  model_diagnostic_dir
)

missing_script05_input_directories <- required_script05_input_directories[
  !dir.exists(required_script05_input_directories)
]

if (length(missing_script05_input_directories) > 0L) {
  stop(
    "Script 05 prerequisite directory(s) not found: ",
    paste(missing_script05_input_directories, collapse = "; "),
    ". Complete Scripts 01-04 before executing Script 05."
  )
}

script05_deliverable_directories <- c(table_dir, figure_dir, report_dir)

for (directory in script05_deliverable_directories) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

missing_script05_deliverable_directories <- script05_deliverable_directories[
  !dir.exists(script05_deliverable_directories)
]

if (length(missing_script05_deliverable_directories) > 0L) {
  stop(
    "Script 05 could not create deliverable directory(s): ",
    paste(missing_script05_deliverable_directories, collapse = "; ")
  )
}

rm(directory)

cross_sectional_path <- file.path(output_dir, "analysis_cross_sectional.parquet")

normalization_analysis_diagnostics_path <- file.path(
  model_diagnostic_dir,
  "normalization_analysis_diagnostics.rds"
)

if (!file.exists(cross_sectional_path) || dir.exists(cross_sectional_path)) {
  stop("Script 05 input must identify an existing analysis file: ", cross_sectional_path)
}

cross_sectional <- read_parquet(cross_sectional_path)

stopifnot(nrow(cross_sectional) > 0L)

cohort_diagnostics_path <- file.path(
  config$cohort_diagnostic_dir,
  "cohort_diagnostics.rds"
)

script04_model_qc_path <- file.path(
  model_diagnostic_dir,
  paste0("model_qc_", site_name, ".rds")
)

message("Script 05 starting | site: ", site_name)

if (!file.exists(cohort_diagnostics_path)) {
  stop("The combined Scripts 01-03 diagnostic RDS was not found: ", cohort_diagnostics_path)
}

if (!file.exists(script04_model_qc_path)) {
  stop("The Script 04 model QC RDS was not found: ", script04_model_qc_path)
}

cohort_diagnostics <- readRDS(cohort_diagnostics_path)
script04_model_qc <- readRDS(script04_model_qc_path)

if (!is.list(cohort_diagnostics) || !"quality_checks_script03" %in% names(cohort_diagnostics)) {
  stop("The combined diagnostic RDS does not contain a valid quality_checks_script03 section.")
}

if (!is.list(script04_model_qc) || !is.list(script04_model_qc$checks)) {
  stop("The Script 04 model QC artifact does not contain a valid checks list.")
}

required_script04_sections <- c("metadata", "upstream_input_contract", "unified_results")
missing_script04_sections <- setdiff(required_script04_sections, names(script04_model_qc))

if (length(missing_script04_sections) > 0L) {
  stop(
    "The Script 04 model QC artifact is missing required section(s): ",
    paste(missing_script04_sections, collapse = ", ")
  )
}

if (
  !is.data.frame(script04_model_qc$upstream_input_contract$summary) ||
  !is.list(script04_model_qc$unified_results$final_artifact) ||
  !is.data.frame(script04_model_qc$unified_results$final_artifact$duplicate_keys)
) {
  stop("The Script 04 model QC artifact does not contain the expected nested structures.")
}

script03_qc <- cohort_diagnostics$quality_checks_script03

script04_cross_sectional_summary <- script04_model_qc$upstream_input_contract$summary %>%
  filter(source == "Cross-sectional analysis")

script04_cross_sectional_n_matches <- nrow(script04_cross_sectional_summary) == 1L &&
  script04_cross_sectional_summary$n_rows[[1]] == nrow(cross_sectional)

script04_unified_artifact <- script04_model_qc$unified_results$final_artifact

script04_unified_artifact_checks <- list(
  unified_row_count_positive = is.numeric(script04_unified_artifact$row_count) &&
    length(script04_unified_artifact$row_count) == 1L &&
    script04_unified_artifact$row_count > 0L,
  unified_analysis_count_positive = is.numeric(script04_unified_artifact$analysis_count) &&
    length(script04_unified_artifact$analysis_count) == 1L &&
    script04_unified_artifact$analysis_count > 0L,
  unified_duplicate_keys_empty = nrow(script04_unified_artifact$duplicate_keys) == 0L
)

script05_upstream_checks <- list(
  script03_overall_status_passed = identical(script03_qc$metadata$overall_status, "PASS"),
  script04_runtime_checks_passed = isTRUE(all(unlist(script04_model_qc$checks, use.names = FALSE))),
  script04_site_matches = identical(script04_model_qc$metadata$site, site_name),
  script04_source_file_matches = identical(
    normalizePath(script04_model_qc$metadata$source_file, mustWork = TRUE),
    normalizePath(cross_sectional_path, mustWork = TRUE)
  ),
  script04_cross_sectional_n_matches = script04_cross_sectional_n_matches,
  script04_unified_results_have_no_duplicate_keys =
    script04_unified_artifact_checks$unified_duplicate_keys_empty,
  script04_unified_row_count_positive =
    script04_unified_artifact_checks$unified_row_count_positive,
  script04_unified_analysis_count_positive =
    script04_unified_artifact_checks$unified_analysis_count_positive
)

stopifnot(isTRUE(all(unlist(script05_upstream_checks, use.names = FALSE))))

normalization_analysis_diagnostics <- list(
  metadata = list(
    script = "Script 05: PBW versus PFVC normalization of physiologic injury metrics",
    site = site_name, generated_at = Sys.time(), source_file = cross_sectional_path,
    outcome = "60-day mortality", outcome_variable = "mortality_event_60"
  ),
  warnings = list(),
  checks = list()
)

duplicate_normalization_check_names <- unique(
  names(normalization_analysis_diagnostics$checks)[duplicated(names(normalization_analysis_diagnostics$checks))]
)

if (length(duplicate_normalization_check_names) > 0L) {
  stop(
    "normalization_analysis_diagnostics contains duplicate check name(s): ",
    paste(duplicate_normalization_check_names, collapse = ", ")
  )
}

normalization_analysis_diagnostics$upstream_input_contract <- list(
  cohort_diagnostics_path = cohort_diagnostics_path,
  script04_model_qc_path = script04_model_qc_path,
  script04_unified_result_summary = script04_model_qc$unified_results$final_artifact,
  checks = script05_upstream_checks
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  script05_upstream_checks
)

save_normalization_diagnostics <- function(diagnostics, path = normalization_analysis_diagnostics_path) {
  stopifnot(
    is.list(diagnostics),
    is.list(diagnostics$checks),
    length(diagnostics$checks) > 0L,
    all(purrr::map_lgl(diagnostics$checks, ~ is.logical(.x) && length(.x) == 1L && !is.na(.x)))
  )
  
  duplicate_check_names <- unique(names(diagnostics$checks)[duplicated(names(diagnostics$checks))])
  
  if (length(duplicate_check_names) > 0L) {
    stop(
      "normalization_analysis_diagnostics contains duplicate check name(s): ",
      paste(duplicate_check_names, collapse = ", ")
    )
  }
  
  saveRDS(diagnostics, path, compress = "xz")
  
  stopifnot(file.exists(path))
  
  diagnostics_readback <- readRDS(path)
  
  stopifnot(
    identical(diagnostics_readback, diagnostics),
    identical(diagnostics_readback$metadata$site, diagnostics$metadata$site),
    identical(names(diagnostics_readback$checks), names(diagnostics$checks)),
    isTRUE(all(unlist(diagnostics_readback$checks, use.names = FALSE)))
  )
  
  invisible(path)
}

# In-sample C-statistic calculated from the observations retained by the fitted model. It is optimistic in absolute terms
# but can support relative comparison when candidate models are fitted to the same outcome observations.
auc_from_model <- function(model) {
  outcome <- model.response(model.frame(model))
  predicted_probability <- fitted(model)
  
  stopifnot(
    length(outcome) == nobs(model),
    length(predicted_probability) == nobs(model),
    all(outcome %in% c(0L, 1L)),
    all(is.finite(predicted_probability)),
    all(predicted_probability >= 0 & predicted_probability <= 1)
  )
  
  n_events <- sum(outcome == 1L)
  n_nonevents <- sum(outcome == 0L)
  
  if (n_events == 0L || n_nonevents == 0L) {
    return(NA_real_)
  }
  
  probability_ranks <- rank(predicted_probability, ties.method = "average")
  
  (sum(probability_ranks[outcome == 1L]) - n_events * (n_events + 1) / 2) / (n_events * n_nonevents)
}

check_subset <- function(data, condition, predicate, empty_ok = TRUE) {
  stopifnot(
    is.data.frame(data),
    is.logical(condition),
    length(condition) == nrow(data),
    !any(is.na(condition)),
    is.function(predicate),
    is.logical(empty_ok),
    length(empty_ok) == 1L,
    !is.na(empty_ok)
  )
  
  subset_data <- data[condition, , drop = FALSE]
  
  if (nrow(subset_data) == 0L) {
    return(empty_ok)
  }
  
  result <- predicate(subset_data)
  
  if (!is.logical(result) || length(result) != 1L || is.na(result)) {
    stop("A subset-validation predicate did not return one nonmissing logical value.")
  }
  
  result
}

# Maximum per-column variance inflation factor for a model design matrix. A one-column design has VIF = 1. A constant
# column or perfectly explained column returns Inf.
vif_max <- function(model) {
  design_matrix <- model.matrix(model)
  design_matrix <- design_matrix[, colnames(design_matrix) != "(Intercept)", drop = FALSE]
  
  if (ncol(design_matrix) == 0L) {
    return(NA_real_)
  }
  
  if (ncol(design_matrix) == 1L) {
    return(1)
  }
  
  vif_values <- vapply(seq_len(ncol(design_matrix)), function(column_index) {
    response <- design_matrix[, column_index]
    predictors <- design_matrix[, -column_index, drop = FALSE]
    auxiliary_design <- cbind(`(Intercept)` = 1, predictors)
    auxiliary_fit <- lm.fit(x = auxiliary_design, y = response)
    
    if (auxiliary_fit$rank < ncol(auxiliary_design) || any(!is.finite(auxiliary_fit$residuals))) {
      return(Inf)
    }
    
    residual_sum_squares <- sum(auxiliary_fit$residuals^2)
    total_sum_squares <- sum((response - mean(response))^2)
    
    if (!is.finite(total_sum_squares) || total_sum_squares <= 0) {
      return(Inf)
    }
    
    r_squared <- 1 - residual_sum_squares / total_sum_squares
    
    if (!is.finite(r_squared) || r_squared >= 1) Inf else 1 / (1 - r_squared)
  }, numeric(1))
  
  max(vif_values)
}

# -----------------------------------------------------------------------------
# QC + modelling frames
# -----------------------------------------------------------------------------

# These are alternative adjustment specifications. Weight/PFVC variables are size ratios, not VT measurements.
# Each outcome model retains exactly one of these terms, alongside its existing mechanic and adjustment covariates.
normalization_adjustment_specs <- tribble(
  ~normalization_spec, ~normalization_variable, ~normalization_label, ~normalization_kind,
  "vt_mclaren_pbw", "vtpbw", "VT/McLaren PBW", "Tidal-volume normalization",
  "vt_nchs_weight", "vtnchs", "VT/NCHS weight", "Tidal-volume normalization",
  "vt_actual_weight", "vtabw", "VT/actual body weight", "Tidal-volume normalization",
  "mclaren_pbw_pfvc", "pbwpfvc", "McLaren PBW/PFVC", "Predicted-size ratio",
  "nchs_weight_pfvc", "nchspfvc", "NCHS weight/PFVC", "Predicted-size ratio",
  "actual_weight_pfvc", "abwpfvc", "Actual body weight/PFVC", "Measured-weight/predicted-volume ratio",
  "vt_pfvc", "vtpfvc", "VT/PFVC", "Tidal-volume normalization"
)

normalization_adjustment_vars <- normalization_adjustment_specs$normalization_variable

normalization_specification_checks <- list(
  normalization_registry_has_seven_specifications = nrow(normalization_adjustment_specs) == 7L,
  normalization_specification_ids_unique = !anyDuplicated(normalization_adjustment_specs$normalization_spec),
  normalization_adjustment_variables_unique = !anyDuplicated(normalization_adjustment_vars),
  normalization_registry_fields_complete = all(complete.cases(normalization_adjustment_specs)),
  normalization_registry_fields_nonempty = all(purrr::map_lgl(
    normalization_adjustment_specs,
    ~ all(nzchar(.x))
  )),
  normalization_registry_matches_requested_variables = setequal(
    normalization_adjustment_vars,
    c("vtpbw", "vtnchs", "vtabw", "pbwpfvc", "nchspfvc", "abwpfvc", "vtpfvc")
  )
)

stopifnot(isTRUE(all(unlist(normalization_specification_checks, use.names = FALSE))))

normalization_required_vars <- unique(c(
  "hospitalization_id", "age_at_admission", "sex_category", "race_category", "pbw", "pbw_nchs", "abw",
  "pfvc", "pbwpfvc",
  normalization_adjustment_vars,
  "bmi", "sofa_total", "sf_ratio", "mortality_event_60", "dp", "dp_raw", "dp_source", "dp_type",
  "dp_invalid_reason", "ers", "ers_pbw", "ers_pfvc", "mechanical_power", "mp_pbw", "mp_pfvc"
))

missing_normalization_vars <- setdiff(normalization_required_vars, names(cross_sectional))

if (length(missing_normalization_vars) > 0L) {
  stop("Script 05 is missing required variable(s): ", paste(missing_normalization_vars, collapse = ", "))
}

expected_sex_levels <- c("Male", "Female")
expected_race_levels <- c("WHITE", "BLACK", "OTHER")

observed_sex_values <- sort(unique(na.omit(as.character(cross_sectional$sex_category))))
observed_race_values <- sort(unique(na.omit(as.character(cross_sectional$race_category))))
observed_mortality_values <- sort(unique(na.omit(cross_sectional$mortality_event_60)))

unexpected_sex_values <- setdiff(observed_sex_values, expected_sex_levels)
unexpected_race_values <- setdiff(observed_race_values, expected_race_levels)

numeric_input_vars <- unique(c(
  "age_at_admission", "pbw", "pbw_nchs", "abw", "pfvc", "pbwpfvc", normalization_adjustment_vars,
  "bmi", "sofa_total", "sf_ratio", "mortality_event_60", "dp", "dp_raw",
  "ers", "ers_pbw", "ers_pfvc", "mechanical_power", "mp_pbw", "mp_pfvc"
))

input_variable_qc <- purrr::map_dfr(normalization_required_vars, function(variable) {
  values <- cross_sectional[[variable]]
  finite_values <- if (is.numeric(values)) values[is.finite(values)] else numeric()
  
  tibble(
    variable = variable, class = paste(class(values), collapse = "/"), n_total = length(values),
    n_nonmissing = sum(!is.na(values)), n_missing = sum(is.na(values)),
    percent_missing = 100 * mean(is.na(values)), n_unique_nonmissing = n_distinct(values, na.rm = TRUE),
    n_infinite = if (is.numeric(values)) sum(is.infinite(values), na.rm = TRUE) else NA_integer_,
    minimum = if (length(finite_values) > 0L) min(finite_values) else NA_real_,
    maximum = if (length(finite_values) > 0L) max(finite_values) else NA_real_
  )
})

input_checks <- list(
  input_file_exists = file.exists(cross_sectional_path),
  cross_sectional_nonempty = nrow(cross_sectional) > 0L,
  hospitalization_ids_complete = !any(is.na(cross_sectional$hospitalization_id)),
  one_row_per_hospitalization = !anyDuplicated(cross_sectional$hospitalization_id),
  required_variables_present = length(missing_normalization_vars) == 0L,
  numeric_variables_numeric = all(purrr::map_lgl(cross_sectional[numeric_input_vars], is.numeric)),
  numeric_variables_have_no_infinite_values = !any(input_variable_qc$n_infinite > 0L, na.rm = TRUE),
  sex_values_valid = length(unexpected_sex_values) == 0L,
  race_values_valid = length(unexpected_race_values) == 0L,
  mortality_complete = !any(is.na(cross_sectional$mortality_event_60)),
  mortality_coding_valid = all(observed_mortality_values %in% c(0L, 1L)),
  mortality_has_events = any(cross_sectional$mortality_event_60 == 1L),
  mortality_has_nonevents = any(cross_sectional$mortality_event_60 == 0L)
)

if (!isTRUE(all(unlist(input_checks, use.names = FALSE)))) {
  print(input_variable_qc, n = Inf, width = Inf)
  print(input_checks)
}

stopifnot(isTRUE(all(unlist(input_checks, use.names = FALSE))))

cross_sectional <- cross_sectional %>%
  mutate(
    .source_row = row_number(),
    sex_category = factor(sex_category, levels = expected_sex_levels),
    race_category = factor(race_category, levels = expected_race_levels),
    age = as.numeric(age_at_admission),
    sf10 = sf_ratio / 10
  )

factor_configuration_checks <- list(
  sex_factor_levels_valid = identical(levels(cross_sectional$sex_category), expected_sex_levels),
  race_factor_levels_valid = identical(levels(cross_sectional$race_category), expected_race_levels)
)

category_coverage_qc <- tibble(
  variable = c("sex_category", "race_category"),
  expected_levels = c(
    paste(expected_sex_levels, collapse = ", "),
    paste(expected_race_levels, collapse = ", ")
  ),
  observed_values = c(
    paste(observed_sex_values, collapse = ", "),
    paste(observed_race_values, collapse = ", ")
  ),
  all_expected_categories_observed = c(
    setequal(observed_sex_values, expected_sex_levels),
    setequal(observed_race_values, expected_race_levels)
  )
)

if (!isTRUE(all(unlist(factor_configuration_checks, use.names = FALSE)))) {
  print(category_coverage_qc, n = Inf, width = Inf)
  print(factor_configuration_checks)
}

stopifnot(isTRUE(all(unlist(factor_configuration_checks, use.names = FALSE))))

sf10_compared <- is.finite(cross_sectional$sf10) & is.finite(cross_sectional$sf_ratio)

sf10_max_absolute_difference <- if (any(sf10_compared)) {
  max(abs(cross_sectional$sf10[sf10_compared] - cross_sectional$sf_ratio[sf10_compared] / 10))
} else {
  NA_real_
}

sf10_checks <- list(
  sf10_has_comparable_values = any(sf10_compared),
  sf10_derivation_valid = is.finite(sf10_max_absolute_difference) &&
    sf10_max_absolute_difference <= 1e-12
)

age_compared <- is.finite(cross_sectional$age) & is.finite(cross_sectional$age_at_admission)

age_max_absolute_difference <- if (any(age_compared)) {
  max(abs(cross_sectional$age[age_compared] - cross_sectional$age_at_admission[age_compared]))
} else {
  NA_real_
}

age_derivation_checks <- list(
  age_has_comparable_values = any(age_compared),
  age_derivation_valid = is.finite(age_max_absolute_difference) &&
    age_max_absolute_difference <= 1e-12
)

stopifnot(
  isTRUE(all(unlist(sf10_checks, use.names = FALSE))),
  isTRUE(all(unlist(age_derivation_checks, use.names = FALSE)))
)

normalization_analysis_diagnostics$input <- list(
  required_variables = normalization_required_vars,
  missing_variables = missing_normalization_vars,
  variable_summary = input_variable_qc,
  observed_sex_values = observed_sex_values,
  observed_race_values = observed_race_values,
  observed_mortality_values = observed_mortality_values,
  category_coverage = category_coverage_qc,
  sf10_derivation = tibble(
    n_compared = sum(sf10_compared),
    maximum_absolute_difference = sf10_max_absolute_difference
  ),
  age_derivation = tibble(
    n_compared = sum(age_compared),
    maximum_absolute_difference = age_max_absolute_difference
  )
)

normalization_analysis_diagnostics$normalization_specifications <- list(
  registry = normalization_adjustment_specs,
  interpretation = paste0(
    "Each requested variable defines a separate alternative adjustment specification in the outcome models. ",
    "Weight/PFVC variables are size ratios, not tidal-volume measurements. ",
    "The PBW-versus-PFVC mechanics normalizers and descriptive reclassification definitions are unchanged."
  ),
  checks = normalization_specification_checks
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  normalization_specification_checks,
  input_checks,
  factor_configuration_checks,
  sf10_checks,
  age_derivation_checks
)

check_derived_variable <- function(observed, expected, derivation, relative_tolerance = 1e-8) {
  expected_finite <- is.finite(expected)
  observed_finite <- is.finite(observed)
  compared <- expected_finite & observed_finite
  absolute_difference <- abs(observed[compared] - expected[compared])
  relative_difference <- absolute_difference / pmax(abs(expected[compared]), .Machine$double.eps)
  
  tibble(
    derivation = derivation, n_expected = sum(expected_finite), n_observed = sum(observed_finite),
    n_compared = sum(compared), n_missing_when_expected = sum(expected_finite & !observed_finite),
    n_observed_when_not_expected = sum(!expected_finite & observed_finite),
    maximum_absolute_difference = if (any(compared)) max(absolute_difference) else NA_real_,
    maximum_relative_difference = if (any(compared)) max(relative_difference) else NA_real_,
    relative_tolerance = relative_tolerance,
    passed = any(compared) && !any(expected_finite & !observed_finite) &&
      !any(!expected_finite & observed_finite) && all(relative_difference <= relative_tolerance)
  )
}

derivation_qc <- bind_rows(
  check_derived_variable(
    observed = cross_sectional$pbwpfvc, expected = cross_sectional$pbw / cross_sectional$pfvc,
    derivation = "pbwpfvc = pbw / pfvc"
  ),
  check_derived_variable(
    observed = cross_sectional$ers_pbw, expected = cross_sectional$ers * cross_sectional$pbw,
    derivation = "ers_pbw = ers * pbw"
  ),
  check_derived_variable(
    observed = cross_sectional$ers_pfvc, expected = cross_sectional$ers * cross_sectional$pfvc,
    derivation = "ers_pfvc = ers * pfvc"
  ),
  check_derived_variable(
    observed = cross_sectional$mp_pbw, expected = cross_sectional$mechanical_power / cross_sectional$pbw,
    derivation = "mp_pbw = mechanical_power / pbw"
  ),
  check_derived_variable(
    observed = cross_sectional$mp_pfvc, expected = cross_sectional$mechanical_power / cross_sectional$pfvc,
    derivation = "mp_pfvc = mechanical_power / pfvc"
  )
)

derivation_checks <- list(
  derived_variables_consistent = all(derivation_qc$passed)
)

if (!isTRUE(all(unlist(derivation_checks, use.names = FALSE)))) {
  print(derivation_qc, n = Inf, width = Inf)
  print(derivation_checks)
}

stopifnot(isTRUE(all(unlist(derivation_checks, use.names = FALSE))))

normalization_analysis_diagnostics$derivations <- derivation_qc
normalization_analysis_diagnostics$checks <- c(normalization_analysis_diagnostics$checks, derivation_checks)

dp_qc <- cross_sectional %>%
  summarise(
    n_total = n(), n_dp_raw_nonmissing = sum(!is.na(dp_raw)),
    n_dp_raw_nonpositive = sum(dp_raw <= 0, na.rm = TRUE), n_dp_valid = sum(!is.na(dp)),
    n_dp_missing = sum(is.na(dp)), n_plateau_derived = sum(dp_source == "plateau_obs", na.rm = TRUE),
    n_pip_derived = sum(dp_source %in% c("pip_obs", "pip_set"), na.rm = TRUE),
    n_nonpositive_raw_retained_as_valid = sum(!is.na(dp_raw) & dp_raw <= 0 & !is.na(dp)),
    n_valid_dp_without_source = sum(!is.na(dp) & is.na(dp_source)),
    n_valid_dp_without_type = sum(!is.na(dp) & is.na(dp_type))
  )

dp_values_match_raw <- with(
  cross_sectional,
  all(is.na(dp) | (is.finite(dp_raw) & abs(dp - dp_raw) <= 1e-12))
)

dp_checks <- list(
  valid_dp_values_positive = all(cross_sectional$dp > 0, na.rm = TRUE),
  valid_dp_matches_raw_dp = dp_values_match_raw,
  nonpositive_raw_dp_not_retained = dp_qc$n_nonpositive_raw_retained_as_valid == 0L,
  valid_dp_has_source = dp_qc$n_valid_dp_without_source == 0L,
  valid_dp_has_type = dp_qc$n_valid_dp_without_type == 0L
)

if (!isTRUE(all(unlist(dp_checks, use.names = FALSE)))) {
  print(dp_qc, width = Inf)
  print(dp_checks)
}

stopifnot(isTRUE(all(unlist(dp_checks, use.names = FALSE))))

# The general discordance cohort requires only valid PBW and PFVC. Metric-specific requirements are imposed in the
# corresponding metric frames, while model covariates are handled later through explicit complete-case samples.
base <- cross_sectional %>%
  filter(is.finite(pbw), pbw > 0, is.finite(pfvc), pfvc > 0, is.finite(pbwpfvc), pbwpfvc > 0)

dp_data <- base %>%
  filter(is.finite(dp), dp > 0)

ers_data <- base %>%
  filter(is.finite(ers), ers > 0, is.finite(ers_pbw), ers_pbw > 0, is.finite(ers_pfvc), ers_pfvc > 0)

mp_data <- base %>%
  filter(
    is.finite(mechanical_power), mechanical_power > 0,
    is.finite(mp_pbw), mp_pbw > 0, is.finite(mp_pfvc), mp_pfvc > 0
  )

summarise_analysis_frame <- function(data, frame_name) {
  tibble(
    frame = frame_name, n = nrow(data), n_hospitalizations = n_distinct(data$hospitalization_id),
    n_deaths_60 = sum(data$mortality_event_60 == 1L, na.rm = TRUE),
    n_alive_60 = sum(data$mortality_event_60 == 0L, na.rm = TRUE),
    n_missing_mortality_60 = sum(is.na(data$mortality_event_60)),
    mortality_60_percent = 100 * mean(data$mortality_event_60 == 1L, na.rm = TRUE)
  )
}

analysis_frame_qc <- bind_rows(
  summarise_analysis_frame(base, "PBW/PFVC base"),
  summarise_analysis_frame(dp_data, "Estimated driving pressure"),
  summarise_analysis_frame(ers_data, "Elastance"),
  summarise_analysis_frame(mp_data, "Mechanical power")
)

frame_exclusion_qc <- tibble(
  frame = c("Estimated driving pressure", "Elastance", "Mechanical power"),
  n_base = nrow(base), n_retained = c(nrow(dp_data), nrow(ers_data), nrow(mp_data))
) %>%
  mutate(n_excluded = n_base - n_retained, percent_retained = 100 * n_retained / n_base)

frame_checks <- list(
  base_frame_nonempty = nrow(base) > 0L,
  dp_frame_nonempty = nrow(dp_data) > 0L,
  ers_frame_nonempty = nrow(ers_data) > 0L,
  mp_frame_nonempty = nrow(mp_data) > 0L,
  base_hospitalization_ids_unique = !anyDuplicated(base$hospitalization_id),
  dp_hospitalization_ids_unique = !anyDuplicated(dp_data$hospitalization_id),
  ers_hospitalization_ids_unique = !anyDuplicated(ers_data$hospitalization_id),
  mp_hospitalization_ids_unique = !anyDuplicated(mp_data$hospitalization_id)
)

if (!isTRUE(all(unlist(frame_checks, use.names = FALSE)))) {
  print(analysis_frame_qc, n = Inf, width = Inf)
  print(frame_exclusion_qc, n = Inf, width = Inf)
  print(frame_checks)
}

stopifnot(isTRUE(all(unlist(frame_checks, use.names = FALSE))))

normalization_analysis_diagnostics$driving_pressure <- dp_qc
normalization_analysis_diagnostics$analysis_frames <- list(
  summary = analysis_frame_qc,
  exclusions = frame_exclusion_qc
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  dp_checks,
  frame_checks
)

dp_valid_n <- dp_qc$n_dp_valid
dp_pip_percent <- if (dp_valid_n > 0L) 100 * dp_qc$n_pip_derived / dp_valid_n else NA_real_

normalization_analysis_diagnostics$warnings$initial <- tibble(
  warning = c(
    "Driving-pressure estimates are predominantly PIP-derived",
    "Mechanical-power analysis has sparse 60-day mortality events",
    "Mechanical-power availability is limited"
  ),
  value = c(dp_pip_percent, analysis_frame_qc$n_deaths_60[analysis_frame_qc$frame == "Mechanical power"],
            frame_exclusion_qc$percent_retained[frame_exclusion_qc$frame == "Mechanical power"]),
  unit = c("Percent of valid DP values", "Deaths in the metric frame", "Percent of the PBW/PFVC base frame"),
  interpretation = c(
    paste0(
      dp_qc$n_pip_derived, " of ", dp_valid_n,
      " valid estimated driving-pressure values are PIP-derived fallbacks; ",
      dp_qc$n_plateau_derived, " is plateau-derived."
    ),
    paste0(
      "The broad mechanical-power frame contains ",
      analysis_frame_qc$n_deaths_60[analysis_frame_qc$frame == "Mechanical power"],
      " deaths. Model-specific complete-case samples may contain fewer events."
    ),
    paste0(
      "Mechanical power is available for ", nrow(mp_data), " of ", nrow(base),
      " patients, so mechanical-power findings may reflect a selected subset."
    )
  )
)

message(
  "Script 05 | interpretation cautions: PIP-derived driving pressure, ",
  "sparse mortality events, and selected mechanical-power availability; ",
  "details retained in normalization diagnostics."
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

message("Script 05 | deriving descriptive normalization discordance")


# Descriptive method comparisons only; these objects do not alter the Part 2 outcome-model frames.
descriptive_weight_methods <- tribble(
  ~weight_method, ~weight_label, ~weight_variable, ~ratio_variable, ~method_color,
  "mclaren", "McLaren PBW", "pbw", "pbwpfvc", "#0072B2",
  "nchs", "NCHS predicted weight", "pbw_nchs", "nchspfvc", "#D55E00",
  "abw", "Actual body weight", "abw", "abwpfvc", "#009E73"
)

descriptive_method_data <- purrr::pmap(
  descriptive_weight_methods,
  function(weight_method, weight_label, weight_variable, ratio_variable, method_color) {
    base %>%
      mutate(
        .method_weight = .data[[weight_variable]],
        .method_ratio = .data[[ratio_variable]]
      ) %>%
      filter(
        is.finite(.method_weight), .method_weight > 0,
        is.finite(.method_ratio), .method_ratio > 0
      ) %>%
      mutate(
        ers_weight_descriptive = ers * .method_weight,
        mp_weight_descriptive = mechanical_power / .method_weight
      )
  }
)

names(descriptive_method_data) <- descriptive_weight_methods$weight_method

descriptive_method_qc <- purrr::imap_dfr(
  descriptive_method_data,
  function(data, method_id) {
    registry_row <- descriptive_weight_methods %>%
      filter(.data$weight_method == .env$method_id)
    
    expected_ratio <- data$.method_weight / data$pfvc
    relative_difference <- abs(data$.method_ratio - expected_ratio) /
      pmax(abs(expected_ratio), .Machine$double.eps)
    
    tibble(
      weight_method = method_id,
      weight_label = registry_row$weight_label[[1]],
      n_base = nrow(base),
      n_valid_weight_and_ratio = nrow(data),
      n_excluded_weight_or_ratio = nrow(base) - nrow(data),
      n_positive_elastance = sum(
        is.finite(data$ers) & data$ers > 0 &
          is.finite(data$ers_weight_descriptive) & data$ers_weight_descriptive > 0 &
          is.finite(data$ers_pfvc) & data$ers_pfvc > 0
      ),
      n_positive_mechanical_power = sum(
        is.finite(data$mechanical_power) & data$mechanical_power > 0 &
          is.finite(data$mp_weight_descriptive) & data$mp_weight_descriptive > 0 &
          is.finite(data$mp_pfvc) & data$mp_pfvc > 0
      ),
      maximum_ratio_relative_difference = if (nrow(data) > 0L) max(relative_difference) else NA_real_,
      ratio_derivation_valid = nrow(data) > 0L && all(relative_difference <= 1e-8),
      hospitalization_ids_unique = !anyDuplicated(data$hospitalization_id)
    )
  }
)

descriptive_method_exclusions <- purrr::imap_dfr(
  descriptive_method_data,
  function(data, method_id) {
    base %>%
      filter(!.source_row %in% data$.source_row) %>%
      transmute(
        weight_method = method_id,
        source_row = .source_row,
        hospitalization_id,
        reason = "Method weight or weight/PFVC ratio is missing, nonfinite, or nonpositive"
      )
  }
)

descriptive_method_checks <- list(
  descriptive_methods_registry_has_three_methods = nrow(descriptive_weight_methods) == 3L,
  descriptive_methods_registry_ids_unique = !anyDuplicated(descriptive_weight_methods$weight_method),
  descriptive_methods_data_names_match_registry = identical(
    names(descriptive_method_data),
    descriptive_weight_methods$weight_method
  ),
  descriptive_methods_have_valid_weight_ratio_rows = all(descriptive_method_qc$n_valid_weight_and_ratio > 0L),
  descriptive_methods_ratio_derivations_valid = all(descriptive_method_qc$ratio_derivation_valid),
  descriptive_methods_hospitalization_ids_unique = all(descriptive_method_qc$hospitalization_ids_unique),
  descriptive_methods_exclusion_counts_reconcile = all(purrr::map2_lgl(
    descriptive_method_qc$weight_method,
    descriptive_method_qc$n_excluded_weight_or_ratio,
    ~ sum(descriptive_method_exclusions$weight_method == .x) == .y
  ))
)

if (!isTRUE(all(unlist(descriptive_method_checks, use.names = FALSE)))) {
  print(descriptive_method_qc, n = Inf, width = Inf)
  print(descriptive_method_checks)
}

stopifnot(isTRUE(all(unlist(descriptive_method_checks, use.names = FALSE))))

normalization_analysis_diagnostics$descriptive_weight_comparison <- list(
  registry = descriptive_weight_methods,
  sample_summary = descriptive_method_qc,
  excluded_rows = descriptive_method_exclusions,
  derivations = c(
    elastance = "ers_weight_descriptive = ers * method weight",
    mechanical_power = "mp_weight_descriptive = mechanical_power / method weight"
  ),
  sample_policy = paste0(
    "Method-specific available samples within the existing PBW/PFVC descriptive base. ",
    "Missing method weights do not exclude patients from other methods, the upstream cohort, or Part 2 models. ",
    "Metric-specific positive/finite requirements are applied during reclassification."
  ),
  source_rows = purrr::map(
    descriptive_method_data,
    ~ .x %>% select(.source_row, hospitalization_id)
  ),
  checks = descriptive_method_checks
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  descriptive_method_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# PART 1 -- Physiologic discordance
# =============================================================================
# For multiplicative normalization, PBW-versus-PFVC discordance is determined exactly by PBW/PFVC:
#
#   log(Ers x PFVC) - log(Ers x PBW) = -log(PBW/PFVC)
#   log(MP / PFVC) - log(MP / PBW) = log(PBW/PFVC)
#
# The overall analysis uses all patients with valid positive PBW and PFVC. Metric-specific identity checks use their
# corresponding metric frames and therefore do not impose unrelated driving-pressure or model-covariate requirements.

discordance_summary_path <- file.path(table_dir, paste0("norm_discordance_summary_", site_name, ".csv"))
discordance_demographics_path <- file.path(table_dir, paste0("norm_discordance_demographics_", site_name, ".csv"))

# PBW/PFVC is a unit-bearing conversion factor in kg/L. Dividing it by 10 gives the numerical ratio between VT/PFVC
# expressed as percent and VT/PBW expressed as mL/kg. The scaled ratio is used only for interpretation; it is not a
# biological claim that 10 kg/L is a normal or ideal PBW/PFVC value.
VT_NORMALIZATION_CONVERSION <- 10

disc_overall <- base %>%
  summarise(
    site = site_name, metric = "PBW/PFVC", unit = "kg/L", n = n(),
    median = median(pbwpfvc), q25 = quantile(pbwpfvc, 0.25, names = FALSE),
    q75 = quantile(pbwpfvc, 0.75, names = FALSE), min = min(pbwpfvc), max = max(pbwpfvc),
    median_vt_scale_ratio = median(pbwpfvc / VT_NORMALIZATION_CONVERSION),
    q25_vt_scale_ratio = quantile(pbwpfvc / VT_NORMALIZATION_CONVERSION, 0.25, names = FALSE),
    q75_vt_scale_ratio = quantile(pbwpfvc / VT_NORMALIZATION_CONVERSION, 0.75, names = FALSE),
    median_vtpfvc_at_6_ml_kg = 6 * median(pbwpfvc / VT_NORMALIZATION_CONVERSION),
    q25_vtpfvc_at_6_ml_kg = 6 * quantile(pbwpfvc / VT_NORMALIZATION_CONVERSION, 0.25, names = FALSE),
    q75_vtpfvc_at_6_ml_kg = 6 * quantile(pbwpfvc / VT_NORMALIZATION_CONVERSION, 0.75, names = FALSE)
  )

# Verify the exact mathematical relationships motivating this analysis.
discordance_identity_qc <- bind_rows(
  ers_data %>%
    summarise(
      metric = "Normalized elastance",
      identity = "log(10 * ers_pfvc) - log(ers_pbw) = -log(pbwpfvc / 10)",
      n_compared = n(),
      maximum_absolute_difference = max(abs(
        log(VT_NORMALIZATION_CONVERSION * ers_pfvc) - log(ers_pbw) +
          log(pbwpfvc / VT_NORMALIZATION_CONVERSION)
      ))
    ),
  mp_data %>%
    summarise(
      metric = "Mechanical power",
      identity = "log(mp_pfvc / 10) - log(mp_pbw) = log(pbwpfvc / 10)",
      n_compared = n(),
      maximum_absolute_difference = max(abs(
        log(mp_pfvc / VT_NORMALIZATION_CONVERSION) - log(mp_pbw) -
          log(pbwpfvc / VT_NORMALIZATION_CONVERSION)
      ))
    )
) %>%
  mutate(passed = n_compared > 0L & maximum_absolute_difference <= 1e-10)

# Use an explicit complete-case sample even though these fields are expected to be complete in the Script 03 cohort.
disc_model_vars <- c("pbwpfvc", "age", "sex_category", "race_category")

disc_model_data <- base %>%
  select(.source_row, hospitalization_id, all_of(disc_model_vars)) %>%
  drop_na(all_of(disc_model_vars))

stopifnot(
  nrow(disc_model_data) > 0L,
  n_distinct(disc_model_data$pbwpfvc) >= 2L,
  n_distinct(disc_model_data$age) >= 2L
)

disc_model <- lm(
  log(pbwpfvc) ~ age + sex_category + race_category,
  data = disc_model_data, x = TRUE, y = TRUE, model = TRUE
)

disc_model_matrix <- model.matrix(disc_model)
disc_model_coefficient_table <- summary(disc_model)$coefficients
disc_model_aliased <- is.na(coef(disc_model))

disc_model_qc <- tibble(
  formula = paste(deparse(formula(disc_model)), collapse = " "),
  n_source = nrow(base), n_obs = nobs(disc_model), n_excluded = n_source - n_obs,
  percent_analyzed = 100 * n_obs / n_source, residual_df = df.residual(disc_model),
  model_matrix_columns = ncol(disc_model_matrix), model_matrix_rank = qr(disc_model_matrix)$rank,
  full_rank = qr(disc_model_matrix)$rank == ncol(disc_model_matrix),
  n_aliased_coefficients = sum(disc_model_aliased),
  coefficients_finite = all(is.finite(coef(disc_model)[!disc_model_aliased])),
  standard_errors_finite = all(is.finite(disc_model_coefficient_table[, "Std. Error"])),
  fitted_values_finite = all(is.finite(fitted(disc_model))),
  residuals_finite = all(is.finite(residuals(disc_model)))
)

# Preserve the existing coefficient columns while adding the model sample size and analysis metadata.
disc_coef <- broom::tidy(disc_model, conf.int = TRUE) %>%
  mutate(site = site_name, n_obs = nobs(disc_model), outcome = "log(PBW/PFVC)", .before = 1)

discordance_checks <- list(
  discordance_base_nonempty = nrow(base) > 0L,
  discordance_values_finite = all(is.finite(base$pbwpfvc)),
  discordance_values_positive = all(base$pbwpfvc > 0),
  discordance_summary_n_matches_base = disc_overall$n == nrow(base),
  exact_metric_discordance_identities_valid = all(discordance_identity_qc$passed),
  discordance_demographic_sample_nonempty = nrow(disc_model_data) > 0L,
  discordance_demographic_model_uses_expected_n = nobs(disc_model) == nrow(disc_model_data),
  discordance_demographic_model_matrix_full_rank = disc_model_qc$full_rank,
  discordance_demographic_model_has_no_aliased_coefficients = disc_model_qc$n_aliased_coefficients == 0L,
  discordance_demographic_coefficients_finite = disc_model_qc$coefficients_finite,
  discordance_demographic_standard_errors_finite = disc_model_qc$standard_errors_finite,
  discordance_demographic_fitted_values_finite = disc_model_qc$fitted_values_finite,
  discordance_demographic_residuals_finite = disc_model_qc$residuals_finite,
  discordance_demographic_residual_df_positive = disc_model_qc$residual_df > 0L
)

if (!isTRUE(all(unlist(discordance_checks, use.names = FALSE)))) {
  print(disc_overall, width = Inf)
  print(discordance_identity_qc, n = Inf, width = Inf)
  print(disc_model_qc, width = Inf)
  print(disc_coef, n = Inf, width = Inf)
  print(discordance_checks)
}

stopifnot(isTRUE(all(unlist(discordance_checks, use.names = FALSE))))

write_csv(disc_overall, discordance_summary_path)
write_csv(disc_coef, discordance_demographics_path)

discordance_output_qc <- tibble(
  output = c("Overall discordance summary", "Discordance demographic model"),
  path = c(discordance_summary_path, discordance_demographics_path)
) %>%
  mutate(
    exists = file.exists(path),
    size_bytes = purrr::map_dbl(path, function(output_path) {
      if (file.exists(output_path)) file.info(output_path, extra_cols = FALSE)$size[[1]] else NA_real_
    }),
    nonempty = exists & is.finite(size_bytes) & size_bytes > 0
  )

discordance_output_checks <- list(
  discordance_outputs_exist = all(discordance_output_qc$exists),
  discordance_outputs_nonempty = all(discordance_output_qc$nonempty)
)

if (!isTRUE(all(unlist(discordance_output_checks, use.names = FALSE)))) {
  print(discordance_output_qc, n = Inf, width = Inf)
  print(discordance_output_checks)
}

stopifnot(isTRUE(all(unlist(discordance_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$discordance <- list(
  overall_summary = disc_overall,
  exact_metric_identities = discordance_identity_qc,
  demographic_model = list(
    required_variables = disc_model_vars,
    source_rows = disc_model_data %>% select(.source_row, hospitalization_id),
    diagnostics = disc_model_qc,
    coefficients = disc_coef
  ),
  outputs = discordance_output_qc
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  discordance_checks,
  discordance_output_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)


# =============================================================================
# Change across cohort-relative metric terciles
# =============================================================================
# Terciles are created separately for the PBW- and PFVC-normalized forms using tie-preserving value boundaries. These are
# cohort-relative categories rather than validated clinical injury thresholds. Tied values are never split across
# terciles; group sizes may therefore differ slightly.

AGE_RECLASS_BREAKS <- c(3, 9, 13, 18)
AGE_RECLASS_LABELS <- c("3 to <9", "9 to <13", "13 to <18")

reclassification_output_path <- file.path(
  table_dir, paste0("norm_discordance_reclassification_", site_name, ".csv")
)

assign_tie_preserving_terciles <- function(values) {
  stopifnot(
    is.numeric(values),
    length(values) > 0L,
    all(is.finite(values))
  )
  
  unique_values <- sort(unique(values))
  
  if (length(unique_values) < 3L) {
    stop("At least three unique values are required for tie-preserving terciles.")
  }
  
  value_counts <- tabulate(match(values, unique_values), nbins = length(unique_values))
  cumulative_counts <- cumsum(value_counts)
  candidate_indices <- seq_len(length(unique_values) - 1L)
  target_counts <- length(values) * c(1 / 3, 2 / 3)
  
  cut_indices <- purrr::map_int(target_counts, function(target_count) {
    candidate_indices[[which.min(abs(cumulative_counts[candidate_indices] - target_count))]]
  })
  
  if (cut_indices[[1]] >= cut_indices[[2]]) {
    stop("Tie-preserving tercile boundaries could not produce three ordered groups.")
  }
  
  cut_points <- c(
    mean(unique_values[c(cut_indices[[1]], cut_indices[[1]] + 1L)]),
    mean(unique_values[c(cut_indices[[2]], cut_indices[[2]] + 1L)])
  )
  
  if (!all(is.finite(cut_points)) || cut_points[[1]] >= cut_points[[2]]) {
    stop("Tie-preserving tercile cut points are invalid.")
  }
  
  terciles <- cut(
    values,
    breaks = c(-Inf, cut_points, Inf),
    labels = 1:3,
    include.lowest = TRUE,
    right = TRUE
  ) %>%
    as.integer()
  
  if (any(is.na(terciles)) || !identical(sort(unique(terciles)), 1:3)) {
    stop("Tie-preserving tercile assignment did not produce three complete groups.")
  }
  
  ties_split <- tibble(value = values, tercile = terciles) %>%
    distinct() %>%
    count(value, name = "n_terciles") %>%
    filter(n_terciles > 1L)
  
  if (nrow(ties_split) > 0L) {
    stop("Tie-preserving tercile assignment unexpectedly split tied values.")
  }
  
  terciles
}

reclassify_metric <- function(data, pbw_var, pfvc_var, label) {
  required_vars <- c(
    ".source_row", "hospitalization_id", "age", "sex_category", "race_category", pbw_var, pfvc_var
  )
  
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0L) {
    stop("Reclassification for ", label, " is missing variable(s): ", paste(missing_vars, collapse = ", "))
  }
  
  analysis_data <- data %>%
    select(all_of(required_vars)) %>%
    filter(
      is.finite(.data[[pbw_var]]), .data[[pbw_var]] > 0,
      is.finite(.data[[pfvc_var]]), .data[[pfvc_var]] > 0
    ) %>%
    mutate(
      pbw_tertile = assign_tie_preserving_terciles(.data[[pbw_var]]),
      pfvc_tertile = assign_tie_preserving_terciles(.data[[pfvc_var]]),
      reclassified = pbw_tertile != pfvc_tertile,
      tertile_change = pfvc_tertile - pbw_tertile,
      age_group = cut(
        age, breaks = AGE_RECLASS_BREAKS, labels = AGE_RECLASS_LABELS,
        right = FALSE, include.lowest = TRUE
      )
    )
  
  if (nrow(analysis_data) == 0L) {
    stop("No valid observations are available for reclassification of ", label, ".")
  }
  
  age_outside_range <- analysis_data %>%
    filter(age < min(AGE_RECLASS_BREAKS) | age >= max(AGE_RECLASS_BREAKS))
  
  if (nrow(age_outside_range) > 0L) {
    stop("At least one age falls outside the prespecified pediatric reclassification range for ", label, ".")
  }
  
  overall_summary <- tibble(
    metric = label, group_type = "Overall", group = "All", n = nrow(analysis_data),
    pct_reclassified = 100 * mean(analysis_data$reclassified)
  )
  
  summarise_group <- function(group_variable, group_type_label) {
    analysis_data %>%
      group_by(group = as.character(.data[[group_variable]])) %>%
      summarise(n = n(), pct_reclassified = 100 * mean(reclassified), .groups = "drop") %>%
      mutate(metric = label, group_type = group_type_label)
  }
  
  grouped_summary <- bind_rows(
    summarise_group("age_group", "Age"),
    summarise_group("sex_category", "Sex"),
    summarise_group("race_category", "Race")
  )
  
  movement_summary <- tibble(
    metric = label,
    n = nrow(analysis_data),
    n_unchanged = sum(analysis_data$tertile_change == 0L),
    n_moved_lower = sum(analysis_data$tertile_change < 0L),
    n_moved_higher = sum(analysis_data$tertile_change > 0L),
    n_moved_one_tertile = sum(abs(analysis_data$tertile_change) == 1L),
    n_moved_two_tertiles = sum(abs(analysis_data$tertile_change) == 2L),
    pct_reclassified = 100 * mean(analysis_data$reclassified)
  )
  
  tertile_counts <- bind_rows(
    analysis_data %>%
      count(tertile = pbw_tertile, name = "n") %>%
      mutate(metric = label, normalization = "PBW"),
    analysis_data %>%
      count(tertile = pfvc_tertile, name = "n") %>%
      mutate(metric = label, normalization = "PFVC")
  ) %>%
    select(metric, normalization, tertile, n)
  
  tied_values_split <- bind_rows(
    analysis_data %>%
      transmute(value = .data[[pbw_var]], tertile = pbw_tertile) %>%
      distinct() %>%
      count(value, name = "n_tertiles") %>%
      filter(n_tertiles > 1L) %>%
      mutate(metric = label, normalization = "PBW"),
    analysis_data %>%
      transmute(value = .data[[pfvc_var]], tertile = pfvc_tertile) %>%
      distinct() %>%
      count(value, name = "n_tertiles") %>%
      filter(n_tertiles > 1L) %>%
      mutate(metric = label, normalization = "PFVC")
  ) %>%
    select(metric, normalization, value, n_tertiles)
  
  transition_matrix <- analysis_data %>%
    count(pbw_tertile, pfvc_tertile, name = "n") %>%
    complete(pbw_tertile = 1:3, pfvc_tertile = 1:3, fill = list(n = 0L)) %>%
    mutate(metric = label, .before = 1)
  
  checks <- list(
    analysis_sample_nonempty = nrow(analysis_data) > 0L,
    hospitalization_ids_unique = !anyDuplicated(analysis_data$hospitalization_id),
    metric_values_positive_and_finite = all(
      is.finite(analysis_data[[pbw_var]]) & analysis_data[[pbw_var]] > 0 &
        is.finite(analysis_data[[pfvc_var]]) & analysis_data[[pfvc_var]] > 0
    ),
    pbw_tertiles_complete = !any(is.na(analysis_data$pbw_tertile)),
    pfvc_tertiles_complete = !any(is.na(analysis_data$pfvc_tertile)),
    pbw_has_three_tertiles = identical(sort(unique(analysis_data$pbw_tertile)), 1:3),
    pfvc_has_three_tertiles = identical(sort(unique(analysis_data$pfvc_tertile)), 1:3),
    no_tied_values_split = nrow(tied_values_split) == 0L,
    age_groups_complete = !any(is.na(analysis_data$age_group)),
    pbw_tertile_counts_reconcile = sum(tertile_counts$n[tertile_counts$normalization == "PBW"]) ==
      nrow(analysis_data),
    pfvc_tertile_counts_reconcile = sum(tertile_counts$n[tertile_counts$normalization == "PFVC"]) ==
      nrow(analysis_data),
    transition_matrix_reconciles = sum(transition_matrix$n) == nrow(analysis_data),
    movement_counts_reconcile = movement_summary$n_unchanged + movement_summary$n_moved_lower +
      movement_summary$n_moved_higher == nrow(analysis_data),
    movement_distance_counts_reconcile = movement_summary$n_unchanged + movement_summary$n_moved_one_tertile +
      movement_summary$n_moved_two_tertiles == nrow(analysis_data),
    reclassification_percent_valid = movement_summary$pct_reclassified >= 0 &
      movement_summary$pct_reclassified <= 100
  )
  
  list(
    analysis_data = analysis_data,
    summary = bind_rows(overall_summary, grouped_summary),
    movement_summary = movement_summary,
    tertile_counts = tertile_counts,
    tied_values_split = tied_values_split,
    transition_matrix = transition_matrix,
    checks = checks
  )
}

reclassification_results <- list(
  elastance = reclassify_metric(
    data = ers_data, pbw_var = "ers_pbw", pfvc_var = "ers_pfvc",
    label = "Normalized elastance (Goligher)"
  ),
  mechanical_power = reclassify_metric(
    data = mp_data, pbw_var = "mp_pbw", pfvc_var = "mp_pfvc",
    label = "Mechanical power (Gattinoni)"
  )
)

# Additional descriptive comparisons; preserve the original McLaren results used by the existing Part 1b models.
descriptive_reclassification_results <- purrr::imap(
  descriptive_method_data,
  function(data, method_id) {
    elastance_data <- data %>%
      filter(
        is.finite(ers), ers > 0,
        is.finite(ers_weight_descriptive), ers_weight_descriptive > 0,
        is.finite(ers_pfvc), ers_pfvc > 0
      )
    
    power_data <- data %>%
      filter(
        is.finite(mechanical_power), mechanical_power > 0,
        is.finite(mp_weight_descriptive), mp_weight_descriptive > 0,
        is.finite(mp_pfvc), mp_pfvc > 0
      )
    
    list(
      elastance = reclassify_metric(
        data = elastance_data,
        pbw_var = "ers_weight_descriptive",
        pfvc_var = "ers_pfvc",
        label = "Normalized elastance"
      ),
      mechanical_power = reclassify_metric(
        data = power_data,
        pbw_var = "mp_weight_descriptive",
        pfvc_var = "mp_pfvc",
        label = "Normalized mechanical power"
      )
    )
  }
)

descriptive_reclassification_summary <- purrr::imap_dfr(
  descriptive_reclassification_results,
  function(method_results, method_id) {
    registry_row <- descriptive_weight_methods %>%
      filter(.data$weight_method == .env$method_id)
    
    purrr::map_dfr(method_results, "summary") %>%
      mutate(
        site = site_name,
        weight_method = method_id,
        weight_label = registry_row$weight_label[[1]],
        .before = 1
      ) %>%
      select(site, weight_method, weight_label, metric, group_type, group, n, pct_reclassified)
  }
)

descriptive_reclassification_check_summary <- purrr::imap_dfr(
  descriptive_reclassification_results,
  function(method_results, method_id) {
    purrr::imap_dfr(method_results, function(result, metric_id) {
      tibble(
        weight_method = method_id,
        metric_id = metric_id,
        check = names(result$checks),
        passed = as.logical(unlist(result$checks, use.names = FALSE))
      )
    })
  }
)

descriptive_mclaren_reconciliation <- purrr::map_dfr(
  c("elastance", "mechanical_power"),
  function(metric_id) {
    original_data <- reclassification_results[[metric_id]]$analysis_data
    expanded_data <- descriptive_reclassification_results$mclaren[[metric_id]]$analysis_data
    
    tibble(
      metric_id = metric_id,
      source_rows_match = identical(expanded_data$.source_row, original_data$.source_row),
      weight_terciles_match = identical(expanded_data$pbw_tertile, original_data$pbw_tertile),
      pfvc_terciles_match = identical(expanded_data$pfvc_tertile, original_data$pfvc_tertile),
      reclassification_matches = identical(expanded_data$reclassified, original_data$reclassified)
    )
  }
)

descriptive_reclassification_checks <- list(
  descriptive_reclassification_has_three_methods = identical(
    names(descriptive_reclassification_results),
    descriptive_weight_methods$weight_method
  ),
  descriptive_reclassification_metric_checks_pass = all(descriptive_reclassification_check_summary$passed),
  descriptive_reclassification_has_six_overall_rows = nrow(
    filter(descriptive_reclassification_summary, group_type == "Overall")
  ) == 6L,
  descriptive_reclassification_summary_keys_unique = !anyDuplicated(
    descriptive_reclassification_summary[c("weight_method", "metric", "group_type", "group")]
  ),
  descriptive_reclassification_mclaren_source_rows_preserved = all(
    descriptive_mclaren_reconciliation$source_rows_match
  ),
  descriptive_reclassification_mclaren_weight_terciles_preserved = all(
    descriptive_mclaren_reconciliation$weight_terciles_match
  ),
  descriptive_reclassification_mclaren_pfvc_terciles_preserved = all(
    descriptive_mclaren_reconciliation$pfvc_terciles_match
  ),
  descriptive_reclassification_mclaren_changes_preserved = all(
    descriptive_mclaren_reconciliation$reclassification_matches
  )
)

if (!isTRUE(all(unlist(descriptive_reclassification_checks, use.names = FALSE)))) {
  print(descriptive_reclassification_summary, n = Inf, width = Inf)
  print(descriptive_mclaren_reconciliation, n = Inf, width = Inf)
  print(descriptive_reclassification_checks)
}

stopifnot(isTRUE(all(unlist(descriptive_reclassification_checks, use.names = FALSE))))

normalization_analysis_diagnostics$descriptive_weight_comparison$reclassification <- list(
  interpretation = paste0(
    "Each method compares independently constructed, tie-preserving normalized-metric terciles against PFVC ",
    "normalization on the same method-specific metric sample. These are cohort-relative categories, not validated ",
    "clinical injury thresholds. The original McLaren assignments are checked for exact preservation."
  ),
  summary = descriptive_reclassification_summary,
  metric_checks = descriptive_reclassification_check_summary,
  mclaren_reconciliation = descriptive_mclaren_reconciliation,
  movement_summary = purrr::map(descriptive_reclassification_results, function(method_results) {
    purrr::map(method_results, "movement_summary")
  }),
  tertile_counts = purrr::map(descriptive_reclassification_results, function(method_results) {
    purrr::map(method_results, "tertile_counts")
  }),
  transition_matrices = purrr::map(descriptive_reclassification_results, function(method_results) {
    purrr::map(method_results, "transition_matrix")
  }),
  source_rows = purrr::map(descriptive_reclassification_results, function(method_results) {
    purrr::map(method_results, ~ .x$analysis_data %>% select(.source_row, hospitalization_id))
  }),
  checks = descriptive_reclassification_checks
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  descriptive_reclassification_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

recl_tbl <- purrr::map_dfr(reclassification_results, "summary") %>%
  select(metric, group_type, group, n, pct_reclassified) %>%
  mutate(site = site_name, .before = 1)

reclassification_movement_qc <- purrr::map_dfr(reclassification_results, "movement_summary")
reclassification_tertile_counts <- purrr::map_dfr(reclassification_results, "tertile_counts")
reclassification_tied_values <- purrr::map_dfr(reclassification_results, "tied_values_split")
reclassification_transition_matrices <- purrr::map_dfr(reclassification_results, "transition_matrix")

reclassification_checks <- list(
  elastance_reclassification_checks_passed = isTRUE(all(unlist(
    reclassification_results$elastance$checks, use.names = FALSE
  ))),
  mechanical_power_reclassification_checks_passed = isTRUE(all(unlist(
    reclassification_results$mechanical_power$checks, use.names = FALSE
  ))),
  both_metrics_summarized = n_distinct(recl_tbl$metric) == 2L,
  exactly_one_overall_row_per_metric = recl_tbl %>%
    filter(group_type == "Overall") %>%
    count(metric) %>%
    summarise(valid = n() == 2L && all(n == 1L)) %>%
    pull(valid),
  reclassification_percentages_valid = all(
    recl_tbl$pct_reclassified >= 0 & recl_tbl$pct_reclassified <= 100
  )
)

if (!isTRUE(all(unlist(reclassification_checks, use.names = FALSE)))) {
  print(recl_tbl, n = Inf, width = Inf)
  print(reclassification_movement_qc, n = Inf, width = Inf)
  print(reclassification_tertile_counts, n = Inf, width = Inf)
  print(reclassification_tied_values, n = Inf, width = Inf)
  print(reclassification_transition_matrices, n = Inf, width = Inf)
  print(reclassification_checks)
}

stopifnot(isTRUE(all(unlist(reclassification_checks, use.names = FALSE))))

write_csv(recl_tbl, reclassification_output_path)

reclassification_output_qc <- tibble(
  output = "Cohort-relative tercile-change summary",
  path = reclassification_output_path,
  exists = file.exists(reclassification_output_path),
  size_bytes = if (file.exists(reclassification_output_path)) {
    file.info(reclassification_output_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

reclassification_output_checks <- list(
  reclassification_output_exists = reclassification_output_qc$exists,
  reclassification_output_nonempty = reclassification_output_qc$nonempty
)

stopifnot(isTRUE(all(unlist(reclassification_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$reclassification <- list(
  interpretation = paste0(
    "Terciles are cohort-relative rank categories rather than validated clinical injury thresholds. ",
    "PBW and PFVC terciles were created separately using tie-preserving value boundaries."
  ),
  summary = recl_tbl,
  movement_summary = reclassification_movement_qc,
  tertile_counts = reclassification_tertile_counts,
  tied_values_split_across_tertiles = reclassification_tied_values,
  transition_matrices = reclassification_transition_matrices,
  source_rows = purrr::map(
    reclassification_results,
    ~ .x$analysis_data %>% select(.source_row, hospitalization_id)
  ),
  metric_checks = purrr::map(reclassification_results, "checks"),
  output = reclassification_output_qc
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  reclassification_checks,
  reclassification_output_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# Part 1 figure: three weight methods and method-specific reclassification
# =============================================================================

part1_figure_path <- file.path(figure_dir, paste0("norm_discordance_part1_", site_name, ".pdf"))
RECLASSIFICATION_SMALL_GROUP_N <- 20L

descriptive_method_colors <- setNames(
  descriptive_weight_methods$method_color,
  descriptive_weight_methods$weight_label
)

part1_ratio_plot_data <- purrr::imap_dfr(descriptive_method_data, function(data, method_id) {
  registry_row <- descriptive_weight_methods %>%
    filter(.data$weight_method == .env$method_id)
  
  data %>%
    transmute(
      source_row = .source_row,
      hospitalization_id,
      age,
      weight_method = method_id,
      weight_label = registry_row$weight_label[[1]],
      weight_pfvc_ratio = .method_ratio
    )
}) %>%
  filter(is.finite(age)) %>%
  mutate(weight_label = factor(weight_label, levels = descriptive_weight_methods$weight_label))

reclassification_small_groups <- descriptive_reclassification_summary %>%
  filter(group_type != "Overall", n < RECLASSIFICATION_SMALL_GROUP_N) %>%
  arrange(weight_method, metric, group_type, n)

normalization_analysis_diagnostics$warnings$reclassification <- reclassification_small_groups %>%
  transmute(
    warning = "Tercile-change subgroup has fewer than 20 patients",
    weight_method, weight_label, metric, group_type, group, n, pct_reclassified,
    interpretation = "The subgroup percentage is descriptive and may be unstable."
  )

if (nrow(normalization_analysis_diagnostics$warnings$reclassification) > 0L) {
  message(
    "Script 05 | reclassification caution: ",
    nrow(normalization_analysis_diagnostics$warnings$reclassification),
    " method/metric subgroup summaries have n < ",
    RECLASSIFICATION_SMALL_GROUP_N,
    "; unstable percentages and denominators retained in diagnostics."
  )
}

PART1_RATIO_DISPLAY_LIMITS <- c(0, 30)

part1_ratio_display_qc <- part1_ratio_plot_data %>%
  group_by(weight_method, weight_label) %>%
  summarise(
    n_observations = n(),
    n_below_display_range = sum(weight_pfvc_ratio < PART1_RATIO_DISPLAY_LIMITS[[1]]),
    n_above_display_range = sum(weight_pfvc_ratio > PART1_RATIO_DISPLAY_LIMITS[[2]]),
    n_in_display_range = sum(between(
      weight_pfvc_ratio,
      PART1_RATIO_DISPLAY_LIMITS[[1]],
      PART1_RATIO_DISPLAY_LIMITS[[2]]
    )),
    .groups = "drop"
  )

part1_ratio_display_checks <- list(
  part1_ratio_display_counts_reconcile = all(
    part1_ratio_display_qc$n_below_display_range +
      part1_ratio_display_qc$n_above_display_range +
      part1_ratio_display_qc$n_in_display_range ==
      part1_ratio_display_qc$n_observations
  ),
  part1_ratio_display_qc_covers_all_methods = setequal(
    part1_ratio_display_qc$weight_method,
    descriptive_weight_methods$weight_method
  )
)

if (!isTRUE(all(unlist(part1_ratio_display_checks, use.names = FALSE)))) {
  print(part1_ratio_display_qc, n = Inf, width = Inf)
  print(part1_ratio_display_checks)
}

stopifnot(isTRUE(all(unlist(part1_ratio_display_checks, use.names = FALSE))))

p_disc_age <- ggplot(
  part1_ratio_plot_data,
  aes(x = age, y = weight_pfvc_ratio, color = weight_label, fill = weight_label)
) +
  geom_point(alpha = 0.12, size = 0.8) +
  geom_smooth(
    method = "loess", formula = y ~ x, se = TRUE, level = 0.95,
    linewidth = 0.9, alpha = 0.12
  ) +
  scale_color_manual(values = descriptive_method_colors, name = NULL, drop = FALSE) +
  scale_fill_manual(values = descriptive_method_colors, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = seq(3, 17, by = 2)) +
  coord_cartesian(ylim = PART1_RATIO_DISPLAY_LIMITS, expand = FALSE) +
  labs(
    title = "Weight/PFVC conversion factors across age",
    subtitle = str_wrap(
      paste0(site_name, " — descriptive age curves for McLaren PBW, NCHS predicted weight, and actual body weight."),
      width = 95
    ),
    x = "Age at admission (years)",
    y = "Weight/PFVC conversion factor (kg/L)",
    caption = str_wrap(
      paste0(
        "Points represent available patient-method observations. Shading shows approximate 95% confidence bands ",
        "around the LOESS curves, not patient-level reference ranges. The visible range is 0–30 kg/L; ",
        "observations outside this range remain included in curve estimation. Samples are method-specific. ",
        "Actual body weight is measured weight, not predicted body weight."
      ),
      width = 105
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 8)),
    plot.caption = element_text(hjust = 0, margin = margin(t = 8)),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
  )

reclassification_group_levels <- c(
  "3 to <9", "9 to <13", "13 to <18", "Male", "Female", "WHITE", "BLACK", "OTHER"
)

reclassification_plot_data <- descriptive_reclassification_summary %>%
  filter(group_type != "Overall") %>%
  mutate(
    group_type = factor(group_type, levels = c("Age", "Sex", "Race")),
    group = factor(group, levels = rev(reclassification_group_levels))
  )

# Use a common percentage axis across all three method blocks.
part1_reclassification_x_max <- min(
  120,
  max(10, max(reclassification_plot_data$pct_reclassified) * 1.22)
)

part1_method_panels <- purrr::map(
  descriptive_weight_methods$weight_method,
  function(method_id) {
    registry_row <- descriptive_weight_methods %>%
      filter(.data$weight_method == .env$method_id)
    
    plot_data <- reclassification_plot_data %>%
      filter(.data$weight_method == .env$method_id)
    
    if (nrow(plot_data) == 0L) {
      stop("No demographic reclassification plot data are available for method: ", method_id)
    }
    
    ggplot(plot_data, aes(x = pct_reclassified, y = group, fill = metric)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_text(
        aes(label = paste0(sprintf("%.1f", pct_reclassified), "% (n=", n, ")")),
        position = position_dodge(width = 0.7), hjust = -0.08, size = 2.7
      ) +
      facet_grid(group_type ~ ., scales = "free_y", space = "free_y") +
      scale_fill_manual(
        values = c("Normalized elastance" = "#009E73", "Normalized mechanical power" = "#E69F00"),
        name = NULL,
        drop = FALSE
      ) +
      scale_x_continuous(
        limits = c(0, part1_reclassification_x_max),
        expand = expansion(mult = c(0, 0.02))
      ) +
      labs(
        title = paste0(registry_row$weight_label[[1]], ": changes in cohort-relative metric terciles"),
        subtitle = str_wrap(
          "Comparison with PFVC normalization. Labels show the percentage reclassified and subgroup denominator.",
          width = 95
        ),
        x = "Patients changing cohort-relative tercile (%)",
        y = NULL,
        caption = str_wrap(
          paste0(
            "Terciles are independently constructed using tie-preserving boundaries, not clinical injury thresholds. ",
            "Method-specific metric samples may differ. Subgroups with n < 20 are retained and flagged in diagnostics."
          ),
          width = 105
        )
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", color = registry_row$method_color[[1]]),
        plot.subtitle = element_text(margin = margin(b = 8)),
        plot.caption = element_text(hjust = 0, margin = margin(t = 8)),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(t = 8, r = 20, b = 8, l = 12)
      )
  }
)

names(part1_method_panels) <- descriptive_weight_methods$weight_method

# Retain reference plot objects while adding the explicitly paginated presentation.
p_recl <- part1_method_panels$mclaren
part1_figure <- p_disc_age / p_recl + patchwork::plot_layout(heights = c(1, 1))

part1_figure_pages <- list(
  page1 = part1_figure,
  page2 = part1_method_panels$nchs / part1_method_panels$abw +
    patchwork::plot_layout(heights = c(1, 1))
)

write_part1_multipage_pdf <- function(pages, path) {
  stopifnot(length(pages) == 2L)
  
  grDevices::pdf(file = path, width = 11, height = 14, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  for (page in pages) {
    print(page)
  }
  
  invisible(path)
}

write_part1_multipage_pdf(part1_figure_pages, part1_figure_path)

part1_figure_qc <- tibble(
  output = "Part 1 three-method discordance and reclassification figure",
  path = part1_figure_path,
  exists = file.exists(part1_figure_path),
  size_bytes = if (file.exists(part1_figure_path)) {
    file.info(part1_figure_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

part1_figure_checks <- list(
  part1_discordance_plot_created = !is.null(p_disc_age),
  part1_reclassification_plot_created = !is.null(p_recl),
  part1_combined_figure_created = !is.null(part1_figure),
  part1_three_method_panels_created = identical(
    names(part1_method_panels),
    descriptive_weight_methods$weight_method
  ),
  part1_two_page_objects_created = length(part1_figure_pages) == 2L,
  part1_ratio_plot_keys_unique = !anyDuplicated(part1_ratio_plot_data[c("weight_method", "source_row")]),
  part1_ratio_plot_values_valid = all(
    is.finite(part1_ratio_plot_data$age) &
      is.finite(part1_ratio_plot_data$weight_pfvc_ratio) &
      part1_ratio_plot_data$weight_pfvc_ratio > 0
  ),
  part1_reclassification_plot_groups_complete = !any(is.na(reclassification_plot_data$group)),
  part1_figure_exists = part1_figure_qc$exists,
  part1_figure_nonempty = part1_figure_qc$nonempty
)

if (!isTRUE(all(unlist(part1_figure_checks, use.names = FALSE)))) {
  print(part1_figure_qc, n = Inf, width = Inf)
  print(part1_figure_checks)
}

stopifnot(isTRUE(all(unlist(part1_figure_checks, use.names = FALSE))))

normalization_analysis_diagnostics$discordance$figure <- part1_figure_qc
normalization_analysis_diagnostics$reclassification$small_group_threshold <- RECLASSIFICATION_SMALL_GROUP_N

normalization_analysis_diagnostics$descriptive_weight_comparison$figure <- list(
  ratio_plot_data = part1_ratio_plot_data,
  ratio_display_limits = PART1_RATIO_DISPLAY_LIMITS,
  ratio_display_diagnostics = part1_ratio_display_qc,
  ratio_display_checks = part1_ratio_display_checks,
  reclassification_plot_data = reclassification_plot_data,
  layout = tibble(
    page = c(1L, 1L, 2L, 2L),
    block = c("Top", "Bottom", "Top", "Bottom"),
    content = c("Three-method age curves", "McLaren reclassification", "NCHS reclassification", "ABW reclassification")
  ),
  method_colors = descriptive_method_colors,
  width_inches = 11,
  height_inches = 14,
  small_group_threshold = RECLASSIFICATION_SMALL_GROUP_N,
  output = part1_figure_qc,
  checks = part1_figure_checks,
  visual_review_status = "PDF page count and readability require inspection; runtime checks validate page objects"
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  part1_ratio_display_checks,
  part1_figure_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

message("Script 05 | descriptive figure saved; deriving age patterns")

# =============================================================================
# PART 1b -- Continuous-age pattern and demographic heterogeneity of PBW/PFVC
# =============================================================================
# Because log(PBW/PFVC) is modeled, exponentiated age coefficients represent multiplicative changes in the conversion
# factor. These analyses describe relationships among prediction formulas and demographic inputs; they are not causal
# estimates of age, sex, or race effects on measured lung volume.

discordance_age_output_path <- file.path(table_dir, paste0("norm_discordance_age_", site_name, ".csv"))
discordance_heterogeneity_output_path <- file.path(
  table_dir, paste0("norm_discordance_age_heterogeneity_", site_name, ".csv")
)

disc_age_model_data <- base %>%
  select(.source_row, hospitalization_id, pbwpfvc, age, sex_category, race_category) %>%
  drop_na(pbwpfvc, age, sex_category, race_category)

stopifnot(
  nrow(disc_age_model_data) > 0L,
  all(is.finite(disc_age_model_data$pbwpfvc)),
  all(disc_age_model_data$pbwpfvc > 0),
  n_distinct(disc_age_model_data$age) >= 4L
)

# Compare a linear age function with a three-degree-of-freedom natural spline on the same observations.
m_disc_lin <- lm(log(pbwpfvc) ~ age, data = disc_age_model_data, x = TRUE, y = TRUE, model = TRUE)
m_disc_spl <- lm(log(pbwpfvc) ~ ns(age, df = 3), data = disc_age_model_data, x = TRUE, y = TRUE, model = TRUE)

disc_spl_comparison <- anova(m_disc_lin, m_disc_spl)
disc_nonlinearity_p <- disc_spl_comparison$`Pr(>F)`[2]

disc_age_tbl <- broom::tidy(m_disc_lin, conf.int = TRUE) %>%
  filter(term == "age") %>%
  transmute(
    site = site_name,
    test = "PBW/PFVC conversion factor versus continuous age",
    n_obs = nobs(m_disc_lin),
    age_unit = "year",
    per_year_pct = 100 * (exp(estimate) - 1),
    lo_pct = 100 * (exp(conf.low) - 1),
    hi_pct = 100 * (exp(conf.high) - 1),
    per_5_year_pct = 100 * (exp(5 * estimate) - 1),
    per_5_year_lo_pct = 100 * (exp(5 * conf.low) - 1),
    per_5_year_hi_pct = 100 * (exp(5 * conf.high) - 1),
    p_value = p.value,
    nonlinearity_p = disc_nonlinearity_p
  )

disc_age_model_qc <- bind_rows(
  tibble(
    model = "Linear age", formula = paste(deparse(formula(m_disc_lin)), collapse = " "),
    n_obs = nobs(m_disc_lin), residual_df = df.residual(m_disc_lin),
    model_matrix_columns = ncol(model.matrix(m_disc_lin)),
    model_matrix_rank = qr(model.matrix(m_disc_lin))$rank,
    full_rank = qr(model.matrix(m_disc_lin))$rank == ncol(model.matrix(m_disc_lin)),
    coefficients_finite = all(is.finite(coef(m_disc_lin))),
    standard_errors_finite = all(is.finite(summary(m_disc_lin)$coefficients[, "Std. Error"])),
    residuals_finite = all(is.finite(residuals(m_disc_lin)))
  ),
  tibble(
    model = "Spline age", formula = paste(deparse(formula(m_disc_spl)), collapse = " "),
    n_obs = nobs(m_disc_spl), residual_df = df.residual(m_disc_spl),
    model_matrix_columns = ncol(model.matrix(m_disc_spl)),
    model_matrix_rank = qr(model.matrix(m_disc_spl))$rank,
    full_rank = qr(model.matrix(m_disc_spl))$rank == ncol(model.matrix(m_disc_spl)),
    coefficients_finite = all(is.finite(coef(m_disc_spl))),
    standard_errors_finite = all(is.finite(summary(m_disc_spl)$coefficients[, "Std. Error"])),
    residuals_finite = all(is.finite(residuals(m_disc_spl)))
  )
)

# Test whether the continuous-age association differs jointly by sex or race.
m_het0 <- lm(
  log(pbwpfvc) ~ age + sex_category + race_category,
  data = disc_age_model_data, x = TRUE, y = TRUE, model = TRUE
)

m_het <- lm(
  log(pbwpfvc) ~ age * sex_category + age * race_category,
  data = disc_age_model_data, x = TRUE, y = TRUE, model = TRUE
)

het_comparison <- anova(m_het0, m_het)
linear_joint_heterogeneity_p <- het_comparison$`Pr(>F)`[2]

# Also test whether the spline-age pattern differs jointly by sex or race.
# This analysis does not assume that nonlinearity was established in the current sample.
# Spline-basis interaction coefficients are not individually interpretable.
m_het_spline0 <- lm(
  log(pbwpfvc) ~ ns(age, df = 3) + sex_category + race_category,
  data = disc_age_model_data, x = TRUE, y = TRUE, model = TRUE
)

m_het_spline <- lm(
  log(pbwpfvc) ~ ns(age, df = 3) * sex_category + ns(age, df = 3) * race_category,
  data = disc_age_model_data, x = TRUE, y = TRUE, model = TRUE
)

het_spline_comparison <- anova(m_het_spline0, m_het_spline)
spline_joint_heterogeneity_p <- het_spline_comparison$`Pr(>F)`[2]

# Interaction coefficients represent differences in the age slope from the reference group, not the total age slope
# within that group. Male and white are the reference categories.
het_tbl <- broom::tidy(m_het, conf.int = TRUE) %>%
  filter(str_detect(term, "(^age:)|(:age$)")) %>%
  transmute(
    site = site_name, n_obs = nobs(m_het), reference_sex = "Male", reference_race = "WHITE", term,
    interpretation = "Difference in linear annual log(PBW/PFVC) slope from the reference category",
    slope_difference_pct_per_year = 100 * (exp(estimate) - 1),
    lo_pct = 100 * (exp(conf.low) - 1), hi_pct = 100 * (exp(conf.high) - 1),
    p_value = p.value, linear_joint_heterogeneity_p = linear_joint_heterogeneity_p,
    spline_joint_heterogeneity_p = spline_joint_heterogeneity_p
  )

heterogeneity_models <- list(
  linear_main_effects = m_het0,
  linear_interactions = m_het,
  spline_main_effects = m_het_spline0,
  spline_interactions = m_het_spline
)

heterogeneity_model_qc <- purrr::imap_dfr(heterogeneity_models, function(fitted_model, model_name) {
  model_matrix <- model.matrix(fitted_model)
  coefficient_table <- summary(fitted_model)$coefficients
  aliased_coefficients <- is.na(coef(fitted_model))
  
  tibble(
    model = model_name, formula = paste(deparse(formula(fitted_model)), collapse = " "),
    n_obs = nobs(fitted_model), residual_df = df.residual(fitted_model),
    model_matrix_columns = ncol(model_matrix), model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    n_aliased_coefficients = sum(aliased_coefficients),
    coefficients_finite = all(is.finite(coef(fitted_model)[!aliased_coefficients])),
    standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
    residuals_finite = all(is.finite(residuals(fitted_model)))
  )
})

discordance_age_checks <- list(
  discordance_age_sample_nonempty = nrow(disc_age_model_data) > 0L,
  discordance_age_sample_ids_unique = !anyDuplicated(disc_age_model_data$hospitalization_id),
  linear_and_spline_models_use_same_n = nobs(m_disc_lin) == nobs(m_disc_spl),
  linear_and_spline_models_use_expected_n = nobs(m_disc_lin) == nrow(disc_age_model_data),
  linear_and_spline_model_matrices_full_rank = all(disc_age_model_qc$full_rank),
  linear_and_spline_coefficients_finite = all(disc_age_model_qc$coefficients_finite),
  linear_and_spline_standard_errors_finite = all(disc_age_model_qc$standard_errors_finite),
  linear_and_spline_residuals_finite = all(disc_age_model_qc$residuals_finite),
  nonlinear_comparison_has_expected_rows = nrow(disc_spl_comparison) == 2L,
  nonlinearity_p_value_valid = is.finite(disc_nonlinearity_p) &&
    disc_nonlinearity_p >= 0 && disc_nonlinearity_p <= 1,
  heterogeneity_models_use_same_n = n_distinct(heterogeneity_model_qc$n_obs) == 1L,
  heterogeneity_models_use_expected_n = all(heterogeneity_model_qc$n_obs == nrow(disc_age_model_data)),
  heterogeneity_model_matrices_full_rank = all(heterogeneity_model_qc$full_rank),
  heterogeneity_models_have_no_aliased_coefficients = all(heterogeneity_model_qc$n_aliased_coefficients == 0L),
  heterogeneity_coefficients_finite = all(heterogeneity_model_qc$coefficients_finite),
  heterogeneity_standard_errors_finite = all(heterogeneity_model_qc$standard_errors_finite),
  heterogeneity_residuals_finite = all(heterogeneity_model_qc$residuals_finite),
  linear_heterogeneity_comparison_has_expected_rows = nrow(het_comparison) == 2L,
  spline_heterogeneity_comparison_has_expected_rows = nrow(het_spline_comparison) == 2L,
  linear_joint_heterogeneity_p_value_valid = is.finite(linear_joint_heterogeneity_p) &&
    linear_joint_heterogeneity_p >= 0 && linear_joint_heterogeneity_p <= 1,
  spline_joint_heterogeneity_p_value_valid = is.finite(spline_joint_heterogeneity_p) &&
    spline_joint_heterogeneity_p >= 0 && spline_joint_heterogeneity_p <= 1,
  expected_linear_interaction_terms_exported = nrow(het_tbl) == 3L
)

if (!isTRUE(all(unlist(discordance_age_checks, use.names = FALSE)))) {
  print(disc_age_tbl, n = Inf, width = Inf)
  print(disc_age_model_qc, n = Inf, width = Inf)
  print(disc_spl_comparison)
  print(het_tbl, n = Inf, width = Inf)
  print(heterogeneity_model_qc, n = Inf, width = Inf)
  print(het_comparison)
  print(het_spline_comparison)
  print(discordance_age_checks)
}

stopifnot(isTRUE(all(unlist(discordance_age_checks, use.names = FALSE))))

write_csv(disc_age_tbl, discordance_age_output_path)
write_csv(het_tbl, discordance_heterogeneity_output_path)

discordance_age_output_qc <- tibble(
  output = c("Continuous-age discordance summary", "Age-slope heterogeneity summary"),
  path = c(discordance_age_output_path, discordance_heterogeneity_output_path)
) %>%
  mutate(
    exists = file.exists(path),
    size_bytes = purrr::map_dbl(path, function(output_path) {
      if (file.exists(output_path)) file.info(output_path, extra_cols = FALSE)$size[[1]] else NA_real_
    }),
    nonempty = exists & is.finite(size_bytes) & size_bytes > 0
  )

discordance_age_output_checks <- list(
  discordance_age_outputs_exist = all(discordance_age_output_qc$exists),
  discordance_age_outputs_nonempty = all(discordance_age_output_qc$nonempty)
)

if (!isTRUE(all(unlist(discordance_age_output_checks, use.names = FALSE)))) {
  print(discordance_age_output_qc, n = Inf, width = Inf)
  print(discordance_age_output_checks)
}

stopifnot(isTRUE(all(unlist(discordance_age_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$discordance$age_pattern <- list(
  sample_n = nrow(disc_age_model_data),
  source_rows = disc_age_model_data %>% select(.source_row, hospitalization_id),
  linear_summary = disc_age_tbl,
  model_diagnostics = disc_age_model_qc,
  linear_vs_spline_comparison = as_tibble(disc_spl_comparison, rownames = "model_number"),
  heterogeneity = list(
    coefficient_summary = het_tbl,
    model_diagnostics = heterogeneity_model_qc,
    linear_nested_model_comparison = as_tibble(het_comparison, rownames = "model_number"),
    spline_nested_model_comparison = as_tibble(het_spline_comparison, rownames = "model_number"),
    interpretation = paste0(
      "Linear interaction coefficients describe differences in annual log-linear slopes. ",
      "The spline interaction comparison evaluates demographic heterogeneity while allowing a nonlinear age pattern. ",
      "Interpret it alongside the current sample's linear-versus-spline comparison."
    )
  ),
  outputs = discordance_age_output_qc
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  discordance_age_checks,
  discordance_age_output_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)


# =============================================================================
# Cohort-relative tercile change versus continuous age
# =============================================================================
# Reuse the patient-level tertile assignments generated above rather than independently recalculating them. The age
# coefficient is an adjusted logit-linear association per year, conditional on sex and race.

reclassification_age_output_path <- file.path(
  table_dir, paste0("norm_reclassification_age_", site_name, ".csv")
)

reclassification_age_data <- purrr::imap(reclassification_results, function(result, result_name) {
  result$analysis_data %>%
    transmute(
      .source_row, hospitalization_id, reclassified = as.integer(reclassified),
      age, sex_category, race_category
    )
})

reclassification_age_labels <- c(
  elastance = "Normalized elastance (Goligher)",
  mechanical_power = "Mechanical power (Gattinoni)"
)

recl_age_fit <- purrr::map(reclassification_age_data, function(analysis_data) {
  glm(
    reclassified ~ age + sex_category + race_category,
    data = analysis_data, family = binomial, x = TRUE, y = TRUE, model = TRUE
  )
})

stopifnot(
  identical(names(reclassification_age_data), names(reclassification_age_labels)),
  identical(names(recl_age_fit), names(reclassification_age_labels))
)

# Conservative criteria for displaying model-based tercile-change curves. These are visualization heuristics rather
# than absolute model-validity thresholds.
MIN_RECLASSIFICATION_EVENTS_PER_PARAMETER <- 10
EXTREME_FITTED_PROBABILITY_THRESHOLD <- 1e-6

recl_age_model_qc <- purrr::imap_dfr(recl_age_fit, function(fitted_model, model_name) {
  outcome <- model.response(model.frame(fitted_model))
  model_matrix <- model.matrix(fitted_model)
  coefficient_table <- summary(fitted_model)$coefficients
  aliased_coefficients <- is.na(coef(fitted_model))
  n_events <- sum(outcome == 1L)
  n_nonevents <- sum(outcome == 0L)
  n_parameters_excluding_intercept <- fitted_model$rank - 1L
  
  tibble(
    metric_id = model_name, metric = unname(reclassification_age_labels[[model_name]]),
    formula = paste(deparse(formula(fitted_model)), collapse = " "),
    n_obs = nobs(fitted_model), n_reclassified = n_events, n_not_reclassified = n_nonevents,
    reclassification_percent = 100 * mean(outcome == 1L),
    n_parameters_excluding_intercept = n_parameters_excluding_intercept,
    events_per_parameter = n_events / n_parameters_excluding_intercept,
    converged = isTRUE(fitted_model$converged), boundary = isTRUE(fitted_model$boundary),
    iterations = fitted_model$iter, model_matrix_columns = ncol(model_matrix),
    model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    n_aliased_coefficients = sum(aliased_coefficients),
    coefficients_finite = all(is.finite(coef(fitted_model)[!aliased_coefficients])),
    standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
    minimum_fitted_probability = min(fitted(fitted_model)),
    maximum_fitted_probability = max(fitted(fitted_model)),
    n_extreme_fitted_probabilities = sum(
      fitted(fitted_model) < EXTREME_FITTED_PROBABILITY_THRESHOLD |
        fitted(fitted_model) > 1 - EXTREME_FITTED_PROBABILITY_THRESHOLD
    )
  )
})

recl_age_tbl <- purrr::imap_dfr(recl_age_fit, function(fitted_model, model_name) {
  age_coefficient <- summary(fitted_model)$coefficients["age", ]
  critical_value <- qnorm(0.975)
  metric_label <- unname(reclassification_age_labels[[model_name]])
  
  model_diagnostics <- recl_age_model_qc %>%
    filter(metric_id == model_name)
  
  stopifnot(nrow(model_diagnostics) == 1L)
  
  tibble(
    site = site_name,
    metric = metric_label,
    n_obs = nobs(fitted_model),
    n_reclassified = sum(model.response(model.frame(fitted_model)) == 1L),
    n_parameters_excluding_intercept = model_diagnostics$n_parameters_excluding_intercept,
    events_per_parameter = model_diagnostics$events_per_parameter,
    age_unit = "year",
    reclass_OR_per_year = exp(age_coefficient[["Estimate"]]),
    lo = exp(age_coefficient[["Estimate"]] - critical_value * age_coefficient[["Std. Error"]]),
    hi = exp(age_coefficient[["Estimate"]] + critical_value * age_coefficient[["Std. Error"]]),
    p_value = age_coefficient[["Pr(>|z|)"]],
    inference = "Model-based Wald"
  )
})

reclassification_age_cell_qc <- purrr::imap_dfr(
  reclassification_age_data,
  function(analysis_data, model_name) {
    bind_rows(
      analysis_data %>%
        count(group = as.character(sex_category), reclassified, name = "n") %>%
        complete(group, reclassified = 0:1, fill = list(n = 0L)) %>%
        mutate(grouping_variable = "sex_category"),
      analysis_data %>%
        count(group = as.character(race_category), reclassified, name = "n") %>%
        complete(group, reclassified = 0:1, fill = list(n = 0L)) %>%
        mutate(grouping_variable = "race_category")
    ) %>%
      mutate(
        metric_id = model_name,
        metric = unname(reclassification_age_labels[[model_name]]),
        .before = 1
      )
  }
)

reclassification_age_sample_qc <- purrr::imap_dfr(
  reclassification_age_data,
  function(analysis_data, model_name) {
    expected_result <- reclassification_results[[model_name]]
    
    tibble(
      metric_id = model_name,
      metric = unname(reclassification_age_labels[[model_name]]),
      n_reclassification_source = nrow(expected_result$analysis_data),
      n_age_model_source = nrow(analysis_data),
      n_model = nobs(recl_age_fit[[model_name]]),
      source_rows_match = identical(
        analysis_data$.source_row,
        expected_result$analysis_data$.source_row
      ),
      reclassification_count_source = sum(expected_result$analysis_data$reclassified),
      reclassification_count_model = sum(model.response(model.frame(recl_age_fit[[model_name]])) == 1L)
    )
  }
)

reclassification_age_checks <- list(
  both_reclassification_age_models_fitted = length(recl_age_fit) == 2L,
  reclassification_age_samples_nonempty = all(recl_age_model_qc$n_obs > 0L),
  reclassification_age_models_have_events = all(recl_age_model_qc$n_reclassified > 0L),
  reclassification_age_models_have_nonevents = all(recl_age_model_qc$n_not_reclassified > 0L),
  reclassification_age_models_use_expected_n = all(
    reclassification_age_sample_qc$n_model == reclassification_age_sample_qc$n_reclassification_source
  ),
  reclassification_age_source_rows_match = all(reclassification_age_sample_qc$source_rows_match),
  reclassification_age_event_counts_match = all(
    reclassification_age_sample_qc$reclassification_count_source ==
      reclassification_age_sample_qc$reclassification_count_model
  ),
  reclassification_age_models_converged = all(recl_age_model_qc$converged),
  no_reclassification_age_model_at_boundary = !any(recl_age_model_qc$boundary),
  reclassification_age_model_matrices_full_rank = all(recl_age_model_qc$full_rank),
  reclassification_age_models_have_no_aliased_coefficients = all(
    recl_age_model_qc$n_aliased_coefficients == 0L
  ),
  reclassification_age_coefficients_finite = all(recl_age_model_qc$coefficients_finite),
  reclassification_age_standard_errors_finite = all(recl_age_model_qc$standard_errors_finite),
  reclassification_age_fitted_probabilities_valid = all(
    is.finite(recl_age_model_qc$minimum_fitted_probability) &
      is.finite(recl_age_model_qc$maximum_fitted_probability) &
      recl_age_model_qc$minimum_fitted_probability >= 0 &
      recl_age_model_qc$maximum_fitted_probability <= 1
  ),
  reclassification_age_odds_ratios_finite = all(is.finite(recl_age_tbl$reclass_OR_per_year)),
  reclassification_age_confidence_intervals_finite = all(
    is.finite(recl_age_tbl$lo) & is.finite(recl_age_tbl$hi)
  ),
  reclassification_age_p_values_valid = all(
    recl_age_tbl$p_value >= 0 & recl_age_tbl$p_value <= 1
  )
)

if (!isTRUE(all(unlist(reclassification_age_checks, use.names = FALSE)))) {
  print(recl_age_model_qc, n = Inf, width = Inf)
  print(reclassification_age_sample_qc, n = Inf, width = Inf)
  print(recl_age_tbl, n = Inf, width = Inf)
  print(reclassification_age_checks)
}

stopifnot(isTRUE(all(unlist(reclassification_age_checks, use.names = FALSE))))

write_csv(recl_age_tbl, reclassification_age_output_path)

reclassification_age_output_qc <- tibble(
  output = "Cohort-relative tercile change versus continuous age",
  path = reclassification_age_output_path,
  exists = file.exists(reclassification_age_output_path),
  size_bytes = if (file.exists(reclassification_age_output_path)) {
    file.info(reclassification_age_output_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

reclassification_age_output_checks <- list(
  reclassification_age_output_exists = reclassification_age_output_qc$exists,
  reclassification_age_output_nonempty = reclassification_age_output_qc$nonempty
)

if (!isTRUE(all(unlist(reclassification_age_output_checks, use.names = FALSE)))) {
  print(reclassification_age_output_qc, n = Inf, width = Inf)
  print(reclassification_age_output_checks)
}

stopifnot(isTRUE(all(unlist(reclassification_age_output_checks, use.names = FALSE))))

reclassification_age_warnings <- recl_age_model_qc %>%
  mutate(
    extreme_fitted_probabilities_present = n_extreme_fitted_probabilities > 0L,
    curve_usable = converged & !boundary & full_rank & n_aliased_coefficients == 0L &
      coefficients_finite & standard_errors_finite &
      events_per_parameter >= MIN_RECLASSIFICATION_EVENTS_PER_PARAMETER &
      n_extreme_fitted_probabilities == 0L,
    interpretation = case_when(
      !converged ~ "Curve suppressed: model did not converge.",
      boundary ~ "Curve suppressed: model fit is on the boundary.",
      !full_rank | n_aliased_coefficients > 0L ~ "Curve suppressed: model matrix is rank deficient.",
      !coefficients_finite | !standard_errors_finite ~ "Curve suppressed: model inference is nonfinite.",
      n_extreme_fitted_probabilities > 0L ~
        "Curve suppressed: extreme fitted probabilities indicate possible separation.",
      events_per_parameter < MIN_RECLASSIFICATION_EVENTS_PER_PARAMETER ~ paste0(
        "Curve suppressed: fewer than ", MIN_RECLASSIFICATION_EVENTS_PER_PARAMETER,
        " tercile-change events per fitted parameter."
      ),
      TRUE ~ "Curve displayed: prespecified stability criteria passed."
    )
  ) %>%
  select(
    metric, n_obs, n_reclassified, n_parameters_excluding_intercept, events_per_parameter,
    converged, boundary, full_rank, n_aliased_coefficients, coefficients_finite,
    standard_errors_finite, n_extreme_fitted_probabilities, extreme_fitted_probabilities_present,
    curve_usable, interpretation
  )

required_curve_status_columns <- c(
  "metric", "events_per_parameter", "converged", "boundary", "full_rank",
  "n_aliased_coefficients", "coefficients_finite", "standard_errors_finite",
  "n_extreme_fitted_probabilities", "curve_usable", "interpretation"
)

missing_curve_status_columns <- setdiff(required_curve_status_columns, names(reclassification_age_warnings))

stopifnot(
  length(missing_curve_status_columns) == 0L,
  nrow(reclassification_age_warnings) == length(recl_age_fit),
  !any(is.na(reclassification_age_warnings$curve_usable))
)

if (any(!reclassification_age_warnings$curve_usable)) {
  message(
    "Script 05 | reclassification curves suppressed for: ",
    paste(
      reclassification_age_warnings$metric[
        !reclassification_age_warnings$curve_usable
      ],
      collapse = "; "
    ),
    ". Reasons and model diagnostics retained in the diagnostic RDS."
  )
}

normalization_analysis_diagnostics$reclassification$age_association <- list(
  interpretation = paste0(
    "Odds ratios describe the adjusted logit-linear association between continuous age and changing cohort-relative ",
    "tercile when switching normalization. These are descriptive, not causal, associations."
  ),
  curve_display_criteria = list(
    minimum_events_per_fitted_parameter = MIN_RECLASSIFICATION_EVENTS_PER_PARAMETER,
    requires_convergence = TRUE,
    permits_boundary_fit = FALSE,
    requires_full_rank = TRUE,
    permits_aliased_coefficients = FALSE,
    requires_finite_coefficients_and_standard_errors = TRUE,
    permits_extreme_fitted_probabilities = FALSE,
    extreme_fitted_probability_threshold = EXTREME_FITTED_PROBABILITY_THRESHOLD,
    note = paste0(
      "The events-per-parameter threshold is a conservative visualization heuristic, ",
      "not an absolute model-validity boundary."
    )
  ),
  model_diagnostics = recl_age_model_qc,
  sample_validation = reclassification_age_sample_qc,
  grouped_outcome_counts = reclassification_age_cell_qc,
  coefficient_summary = recl_age_tbl,
  warnings = reclassification_age_warnings,
  source_rows = purrr::map(reclassification_age_data, ~ .x %>% select(.source_row, hospitalization_id)),
  output = reclassification_age_output_qc
)

normalization_analysis_diagnostics$warnings$reclassification_age <- reclassification_age_warnings

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  reclassification_age_checks,
  reclassification_age_output_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# Part 1b figure: nonlinear age pattern and probability of cohort-relative tercile change
# =============================================================================

part1b_figure_path <- file.path(figure_dir, paste0("norm_discordance_part1b_", site_name, ".pdf"))

# Descriptive sex-specific age curves for each weight/PFVC ratio, conditional on the White race reference.
descriptive_age_model_data <- purrr::map(
  descriptive_method_data,
  ~ .x %>%
    select(.source_row, hospitalization_id, age, sex_category, race_category, .method_ratio) %>%
    drop_na(age, sex_category, race_category, .method_ratio)
)

descriptive_age_models <- purrr::map(descriptive_age_model_data, function(data) {
  stopifnot(nrow(data) > 0L, n_distinct(data$age) >= 4L)
  
  lm(
    log(.method_ratio) ~ ns(age, df = 3) * sex_category + ns(age, df = 3) * race_category,
    data = data, x = TRUE, y = TRUE, model = TRUE
  )
})

descriptive_age_model_qc <- purrr::imap_dfr(descriptive_age_models, function(model, method_id) {
  model_matrix <- model.matrix(model)
  coefficient_table <- summary(model)$coefficients
  
  tibble(
    weight_method = method_id,
    n_method_available = nrow(descriptive_method_data[[method_id]]),
    n_analyzed = nobs(model),
    n_excluded_demographics = n_method_available - n_analyzed,
    formula = paste(deparse(formula(model)), collapse = " "),
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    coefficients_finite = all(is.finite(coef(model))),
    standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
    residual_df = df.residual(model)
  )
})

descriptive_age_fit_checks <- list(
  descriptive_age_models_cover_all_methods = identical(
    names(descriptive_age_models),
    descriptive_weight_methods$weight_method
  ),
  descriptive_age_model_matrices_full_rank = all(descriptive_age_model_qc$full_rank),
  descriptive_age_model_coefficients_finite = all(descriptive_age_model_qc$coefficients_finite),
  descriptive_age_model_standard_errors_finite = all(descriptive_age_model_qc$standard_errors_finite),
  descriptive_age_model_residual_df_positive = all(descriptive_age_model_qc$residual_df > 0L),
  descriptive_age_model_samples_match = all(purrr::map2_lgl(
    descriptive_age_models,
    descriptive_age_model_data,
    ~ nobs(.x) == nrow(.y) && identical(rownames(model.frame(.x)), rownames(.y))
  ))
)

stopifnot(isTRUE(all(unlist(descriptive_age_fit_checks, use.names = FALSE))))

discordance_age_curve <- purrr::imap_dfr(descriptive_age_models, function(model, method_id) {
  analysis_data <- descriptive_age_model_data[[method_id]]
  registry_row <- descriptive_weight_methods %>%
    filter(.data$weight_method == .env$method_id)
  
  prediction_grid <- crossing(
    age = seq(min(analysis_data$age), max(analysis_data$age), length.out = 100),
    sex_category = factor(expected_sex_levels, levels = expected_sex_levels)
  ) %>%
    mutate(race_category = factor("WHITE", levels = expected_race_levels))
  
  prediction <- predict(model, newdata = prediction_grid, se.fit = TRUE)
  
  prediction_grid %>%
    mutate(
      weight_method = method_id,
      weight_label = registry_row$weight_label[[1]],
      fitted_log_ratio = as.numeric(prediction$fit),
      standard_error = as.numeric(prediction$se.fit),
      fitted_weight_pfvc = exp(fitted_log_ratio),
      conf_low = exp(fitted_log_ratio - 1.96 * standard_error),
      conf_high = exp(fitted_log_ratio + 1.96 * standard_error)
    )
}) %>%
  mutate(weight_label = factor(weight_label, levels = descriptive_weight_methods$weight_label))

# Retain the existing McLaren curve field names for the original diagnostics and checks.
descriptive_mclaren_age_curve <- discordance_age_curve %>%
  filter(weight_method == "mclaren")

original_mclaren_age_prediction <- predict(
  m_het_spline,
  newdata = descriptive_mclaren_age_curve,
  se.fit = TRUE
)

descriptive_age_prediction_checks <- list(
  descriptive_age_curve_has_expected_rows = nrow(discordance_age_curve) ==
    nrow(descriptive_weight_methods) * length(expected_sex_levels) * 100L,
  descriptive_age_curve_keys_unique = !anyDuplicated(
    discordance_age_curve[c("weight_method", "sex_category", "age")]
  ),
  descriptive_age_curve_values_finite = all(is.finite(c(
    discordance_age_curve$fitted_weight_pfvc,
    discordance_age_curve$conf_low,
    discordance_age_curve$conf_high
  ))),
  descriptive_age_curve_intervals_ordered = all(
    discordance_age_curve$conf_low <= discordance_age_curve$fitted_weight_pfvc &
      discordance_age_curve$fitted_weight_pfvc <= discordance_age_curve$conf_high
  ),
  descriptive_age_curve_reference_race_preserved = all(
    as.character(discordance_age_curve$race_category) == "WHITE"
  ),
  descriptive_age_mclaren_source_rows_preserved = identical(
    descriptive_age_model_data$mclaren$.source_row,
    disc_age_model_data$.source_row
  ),
  descriptive_age_mclaren_predictions_preserved = isTRUE(all.equal(
    descriptive_mclaren_age_curve$fitted_log_ratio,
    as.numeric(original_mclaren_age_prediction$fit),
    tolerance = 1e-8,
    check.attributes = FALSE
  ))
)

stopifnot(isTRUE(all(unlist(descriptive_age_prediction_checks, use.names = FALSE))))

discordance_age_curve <- discordance_age_curve %>%
  mutate(
    fitted_log_pbwpfvc = fitted_log_ratio,
    fitted_pbwpfvc = fitted_weight_pfvc
  )

p_het <- ggplot(
  discordance_age_curve,
  aes(x = age, y = fitted_weight_pfvc, color = weight_label, fill = weight_label)
) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ sex_category, nrow = 1, drop = FALSE) +
  scale_color_manual(values = descriptive_method_colors, name = NULL, drop = FALSE) +
  scale_fill_manual(values = descriptive_method_colors, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = seq(3, 17, by = 2)) +
  labs(
    title = "Weight/PFVC age patterns by method and sex",
    subtitle = str_wrap(
      paste0(
        "Fitted spline curves at the White race reference category. ",
        "Each method uses its available demographic-complete sample."
      ),
      width = 95
    ),
    x = "Age at admission (years)",
    y = "Fitted weight/PFVC conversion factor (kg/L)",
    caption = str_wrap(
      paste0(
        "Bands are approximate pointwise model-based 95% confidence intervals on the exponentiated fitted log-ratio. ",
        "These describe weight/PFVC ratios, not measured lung-volume effects. Actual body weight is not predicted weight."
      ),
      width = 105
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 8)),
    plot.caption = element_text(hjust = 0, margin = margin(t = 8)),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
  )

reclassification_curve_data <- purrr::imap_dfr(recl_age_fit, function(fitted_model, model_name) {
  status <- reclassification_age_warnings %>%
    filter(metric == unname(reclassification_age_labels[[model_name]]))
  
  required_status_columns <- c(
    "converged", "boundary", "full_rank", "n_aliased_coefficients",
    "coefficients_finite", "standard_errors_finite", "events_per_parameter",
    "n_extreme_fitted_probabilities", "curve_usable", "interpretation"
  )
  
  missing_status_columns <- setdiff(required_status_columns, names(status))
  
  if (length(missing_status_columns) > 0L) {
    stop(
      "Curve diagnostics for ", model_name, " are missing column(s): ",
      paste(missing_status_columns, collapse = ", ")
    )
  }
  
  stopifnot(nrow(status) == 1L)
  
  age_grid <- tibble(
    age = seq(min(model.frame(fitted_model)$age), max(model.frame(fitted_model)$age), length.out = 100),
    sex_category = factor("Male", levels = expected_sex_levels),
    race_category = factor("WHITE", levels = expected_race_levels)
  )
  
  model_usable <- isTRUE(status$converged[[1]]) &&
    !isTRUE(status$boundary[[1]]) &&
    isTRUE(status$full_rank[[1]]) &&
    status$n_aliased_coefficients[[1]] == 0L &&
    isTRUE(status$coefficients_finite[[1]]) &&
    isTRUE(status$standard_errors_finite[[1]]) &&
    is.finite(status$events_per_parameter[[1]]) &&
    status$events_per_parameter[[1]] >= MIN_RECLASSIFICATION_EVENTS_PER_PARAMETER &&
    status$n_extreme_fitted_probabilities[[1]] == 0L
  
  stopifnot(
    length(model_usable) == 1L,
    !is.na(model_usable),
    isTRUE(model_usable == status$curve_usable[[1]])
  )
  
  if (!model_usable) {
    return(
      age_grid %>%
        transmute(
          metric_id = model_name,
          metric = unname(reclassification_age_labels[[model_name]]),
          age, probability = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
          estimable = FALSE, plot_status = status$interpretation
        )
    )
  }
  
  prediction <- predict(fitted_model, newdata = age_grid, type = "link", se.fit = TRUE)
  
  age_grid %>%
    transmute(
      metric_id = model_name,
      metric = unname(reclassification_age_labels[[model_name]]),
      age,
      probability = plogis(as.numeric(prediction$fit)),
      conf_low = plogis(as.numeric(prediction$fit) - 1.96 * as.numeric(prediction$se.fit)),
      conf_high = plogis(as.numeric(prediction$fit) + 1.96 * as.numeric(prediction$se.fit)),
      estimable = TRUE,
      plot_status = status$interpretation
    )
})

reclassification_curve_labels <- reclassification_curve_data %>%
  filter(!estimable) %>%
  group_by(metric, plot_status) %>%
  summarise(age = median(age), .groups = "drop") %>%
  mutate(
    probability = 0.25,
    label = str_wrap(plot_status, width = 32)
  )

suppressed_reclassification_metrics <- reclassification_age_warnings %>%
  filter(!curve_usable) %>%
  pull(metric)

reclassification_curve_subtitle <- if (length(suppressed_reclassification_metrics) == 0L) {
  paste0(
    "Predictions are shown for the male and White reference categories. ",
    "All models passed the prespecified curve-display criteria."
  )
} else {
  paste0(
    "Predictions are shown for the male and White reference categories. ",
    "Suppressed curves: ", paste(suppressed_reclassification_metrics, collapse = "; "),
    ". See panel annotations and diagnostics."
  )
}

p_recl_age <- ggplot(reclassification_curve_data, aes(x = age, y = probability)) +
  geom_ribbon(
    data = reclassification_curve_data %>% filter(estimable),
    aes(ymin = conf_low, ymax = conf_high),
    fill = "#009E73", alpha = 0.15
  ) +
  geom_line(
    data = reclassification_curve_data %>% filter(estimable),
    color = "#009E73", linewidth = 0.9
  ) +
  geom_text(
    data = reclassification_curve_labels,
    aes(x = age, y = probability, label = label),
    inherit.aes = FALSE, color = "#B2182B", size = 3.2, lineheight = 1.05
  ) +
  facet_wrap(~ metric, nrow = 1, drop = FALSE) +
  scale_x_continuous(breaks = seq(3, 17, by = 2)) +
  scale_y_continuous(
    limits = c(0, 0.5),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Probability of changing cohort-relative tercile across age",
    subtitle = str_wrap(reclassification_curve_subtitle, width = 95),
    x = "Age at admission (years)",
    y = "Predicted probability of changing tercile",
    caption = str_wrap(
      paste0(
        "Models adjust for sex and race. The plotted age function is logit-linear. ",
        "Terciles are cohort-relative and are not validated clinical injury thresholds."
      ),
      width = 105
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 8)),
    plot.caption = element_text(hjust = 0, margin = margin(t = 8)),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
  )

part1b_figure <- p_het / p_recl_age + patchwork::plot_layout(heights = c(1, 1.15))

ggsave(filename = part1b_figure_path, plot = part1b_figure, width = 12, height = 9)

part1b_figure_qc <- tibble(
  output = "Part 1b age-pattern and tercile-change figure",
  path = part1b_figure_path,
  exists = file.exists(part1b_figure_path),
  size_bytes = if (file.exists(part1b_figure_path)) {
    file.info(part1b_figure_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

part1b_figure_checks <- list(
  discordance_age_curve_complete = all(complete.cases(
    discordance_age_curve[c("age", "fitted_pbwpfvc", "conf_low", "conf_high")]
  )),
  discordance_age_curve_finite = all(is.finite(c(
    discordance_age_curve$fitted_pbwpfvc,
    discordance_age_curve$conf_low,
    discordance_age_curve$conf_high
  ))),
  discordance_age_intervals_ordered = all(
    discordance_age_curve$conf_low <= discordance_age_curve$fitted_pbwpfvc &
      discordance_age_curve$fitted_pbwpfvc <= discordance_age_curve$conf_high
  ),
  reclassification_curve_status_matches_diagnostics = purrr::imap_lgl(
    recl_age_fit,
    function(.fitted_model, model_name) {
      expected_usable <- reclassification_age_warnings %>%
        filter(metric == unname(reclassification_age_labels[[model_name]])) %>%
        pull(curve_usable)
      
      observed_usable <- reclassification_curve_data %>%
        filter(metric_id == model_name) %>%
        distinct(estimable) %>%
        pull(estimable)
      
      length(expected_usable) == 1L && length(observed_usable) == 1L &&
        isTRUE(observed_usable[[1]] == expected_usable[[1]])
    }
  ) %>%
    all(),
  each_reclassification_model_has_curve_rows = all(
    names(recl_age_fit) %in% unique(reclassification_curve_data$metric_id)
  ),
  nonestimable_models_have_plot_labels = all(
    reclassification_age_warnings$metric[!reclassification_age_warnings$curve_usable] %in%
      reclassification_curve_labels$metric
  ),
  part1b_plot_objects_created = !is.null(p_het) && !is.null(p_recl_age) && !is.null(part1b_figure),
  part1b_figure_exists = part1b_figure_qc$exists,
  part1b_figure_nonempty = part1b_figure_qc$nonempty
)

if (!isTRUE(all(unlist(part1b_figure_checks, use.names = FALSE)))) {
  print(part1b_figure_qc, n = Inf, width = Inf)
  print(part1b_figure_checks)
}

stopifnot(isTRUE(all(unlist(part1b_figure_checks, use.names = FALSE))))

normalization_analysis_diagnostics$discordance$age_pattern$prediction_curve <-
  descriptive_mclaren_age_curve %>%
  transmute(
    age, sex_category, race_category,
    fitted_log_pbwpfvc = fitted_log_ratio,
    standard_error,
    fitted_pbwpfvc = fitted_weight_pfvc,
    conf_low, conf_high
  )

normalization_analysis_diagnostics$descriptive_weight_comparison$age_pattern <- list(
  formula = "log(method weight / PFVC) ~ ns(age, df = 3) * sex + ns(age, df = 3) * race",
  reference_race = "WHITE",
  model_diagnostics = descriptive_age_model_qc,
  coefficients = purrr::imap_dfr(descriptive_age_models, function(model, method_id) {
    broom::tidy(model, conf.int = TRUE) %>%
      mutate(weight_method = method_id, .before = 1)
  }),
  source_rows = purrr::map(
    descriptive_age_model_data,
    ~ .x %>% select(.source_row, hospitalization_id)
  ),
  prediction_curve = discordance_age_curve,
  method_colors = descriptive_method_colors,
  figure = part1b_figure_qc,
  interpretation = paste0(
    "Sex facets show three method-specific descriptive spline curves conditional on White race. ",
    "Method samples may differ; the original McLaren predictions are checked for preservation. ",
    "The bottom reclassification panel retains its original McLaren analysis."
  ),
  checks = c(descriptive_age_fit_checks, descriptive_age_prediction_checks)
)

normalization_analysis_diagnostics$reclassification$age_association$prediction_curve <- reclassification_curve_data
normalization_analysis_diagnostics$reclassification$age_association$curve_display_status <- reclassification_age_warnings %>%
  select(metric, curve_usable, interpretation)
normalization_analysis_diagnostics$reclassification$age_association$figure <- part1b_figure_qc

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  descriptive_age_fit_checks,
  descriptive_age_prediction_checks,
  part1b_figure_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)


# =============================================================================
# PART 2 -- Prognostic utility of normalization for 60-day mortality
# =============================================================================

# Each candidate model retains one alternative normalization/size term, illness severity, oxygenation, and BMI.
# Primary-adjusted models additionally include age, sex, and race. These are statistical adjustment specifications,
# not claims that the covariates form a sufficient causal adjustment set.
common_base_cov <- "sofa_total + sf10 + bmi"
demo_cov <- "age + sex_category + race_category"

# Retain the McLaren reference adjustment string for compatibility and preservation checks.
base_cov <- paste("vtpbw +", common_base_cov)

adjustment_specs <- c(
  minimal = base_cov,
  primary = paste(base_cov, "+", demo_cov)
)

adjustment_labels <- c(
  minimal = "Minimally adjusted",
  primary = "Primary adjusted"
)

get_normalization_adjustment <- function(normalization_spec, adjustment_id) {
  registry_row <- normalization_adjustment_specs %>%
    filter(.data$normalization_spec == .env$normalization_spec)
  
  if (nrow(registry_row) != 1L) {
    stop("Expected exactly one normalization specification: ", normalization_spec)
  }
  
  if (length(adjustment_id) != 1L || !adjustment_id %in% names(adjustment_labels)) {
    stop("Unknown adjustment specification: ", paste(adjustment_id, collapse = ", "))
  }
  
  covariates <- paste(registry_row$normalization_variable[[1]], "+", common_base_cov)
  
  if (adjustment_id == "primary") {
    covariates <- paste(covariates, "+", demo_cov)
  }
  
  covariates
}

normalization_covariate_grid <- crossing(
  normalization_adjustment_specs,
  adjustment = names(adjustment_labels)
) %>%
  mutate(
    adjusted = unname(adjustment_labels[adjustment]),
    covariates = purrr::map2_chr(normalization_spec, adjustment, get_normalization_adjustment)
  )

normalization_covariate_checks <- list(
  normalization_covariate_grid_has_fourteen_rows = nrow(normalization_covariate_grid) ==
    nrow(normalization_adjustment_specs) * length(adjustment_labels),
  normalization_covariate_keys_unique = !anyDuplicated(
    normalization_covariate_grid[c("normalization_spec", "adjustment")]
  ),
  each_covariate_spec_contains_exactly_its_requested_alternative = all(purrr::map2_lgl(
    normalization_covariate_grid$covariates,
    normalization_covariate_grid$normalization_variable,
    function(covariates, expected_variable) {
      formula_variables <- all.vars(as.formula(paste("~", covariates)))
      included_alternatives <- intersect(formula_variables, normalization_adjustment_vars)
      
      length(included_alternatives) == 1L && identical(included_alternatives, expected_variable)
    }
  )),
  original_mclaren_minimal_adjustment_preserved = identical(
    get_normalization_adjustment("vt_mclaren_pbw", "minimal"),
    adjustment_specs[["minimal"]]
  ),
  original_mclaren_primary_adjustment_preserved = identical(
    get_normalization_adjustment("vt_mclaren_pbw", "primary"),
    adjustment_specs[["primary"]]
  )
)

stopifnot(isTRUE(all(unlist(normalization_covariate_checks, use.names = FALSE))))

normalization_analysis_diagnostics$normalization_specifications$covariate_grid <- normalization_covariate_grid

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  normalization_covariate_checks
)

prog_specs <- tribble(
  ~family, ~spec, ~exposure,
  "Elastance", "Mechanic only", "log(ers)",
  "Elastance", "PBW-locked (Goligher)", "log(ers_pbw)",
  "Elastance", "PFVC-locked", "log(ers_pfvc)",
  "Elastance", "Separate (Ers + PBW)", "log(ers) + log(pbw)",
  "Elastance", "Separate (Ers + PFVC)", "log(ers) + log(pfvc)",
  "Mechanical power", "Mechanic only", "log(mechanical_power)",
  "Mechanical power", "PBW-locked (Gattinoni)", "log(mp_pbw)",
  "Mechanical power", "PFVC-locked", "log(mp_pfvc)",
  "Mechanical power", "Separate (MP + PBW)", "log(mechanical_power) + log(pbw)",
  "Mechanical power", "Separate (MP + PFVC)", "log(mechanical_power) + log(pfvc)"
)

family_data <- list(
  "Elastance" = ers_data,
  "Mechanical power" = mp_data
)

family_metric_vars <- list(
  "Elastance" = c("ers", "ers_pbw", "ers_pfvc"),
  "Mechanical power" = c("mechanical_power", "mp_pbw", "mp_pfvc")
)

prognostic_common_vars <- unique(c(
  "mortality_event_60", normalization_adjustment_vars,
  "sofa_total", "sf10", "bmi", "age", "sex_category", "race_category", "pbw", "pfvc"
))

prognostic_analysis_data <- purrr::imap(family_data, function(data, family_name) {
  required_vars <- unique(c(
    ".source_row", "hospitalization_id", prognostic_common_vars, family_metric_vars[[family_name]]
  ))
  
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0L) {
    stop("The ", family_name, " prognostic analysis is missing variable(s): ", paste(missing_vars, collapse = ", "))
  }
  
  data %>%
    select(all_of(required_vars)) %>%
    drop_na(all_of(setdiff(required_vars, c(".source_row", "hospitalization_id"))))
})

prognostic_sample_qc <- purrr::imap_dfr(prognostic_analysis_data, function(data, family_name) {
  tibble(
    family = family_name,
    n_metric_frame = nrow(family_data[[family_name]]),
    n_analyzed = nrow(data),
    n_excluded = n_metric_frame - n_analyzed,
    percent_analyzed = 100 * n_analyzed / n_metric_frame,
    n_events = sum(data$mortality_event_60 == 1L),
    n_nonevents = sum(data$mortality_event_60 == 0L),
    event_percent = 100 * mean(data$mortality_event_60 == 1L)
  )
})

prognostic_sample_checks <- list(
  both_family_samples_created = identical(names(prognostic_analysis_data), names(family_data)),
  family_samples_nonempty = all(purrr::map_int(prognostic_analysis_data, nrow) > 0L),
  family_sample_ids_unique = all(purrr::map_lgl(
    prognostic_analysis_data,
    ~ !anyDuplicated(.x$hospitalization_id)
  )),
  family_samples_complete = all(purrr::map_lgl(
    prognostic_analysis_data,
    ~ all(complete.cases(.x))
  )),
  family_samples_have_events = all(prognostic_sample_qc$n_events > 0L),
  family_samples_have_nonevents = all(prognostic_sample_qc$n_nonevents > 0L),
  family_sample_outcomes_binary = all(purrr::map_lgl(
    prognostic_analysis_data,
    ~ all(.x$mortality_event_60 %in% c(0L, 1L))
  ))
)

if (!isTRUE(all(unlist(prognostic_sample_checks, use.names = FALSE)))) {
  print(prognostic_sample_qc, n = Inf, width = Inf)
  print(prognostic_sample_checks)
}

stopifnot(isTRUE(all(unlist(prognostic_sample_checks, use.names = FALSE))))

prognostic_sample_missingness_qc <- purrr::imap_dfr(family_data, function(data, family_name) {
  required_vars <- unique(c(prognostic_common_vars, family_metric_vars[[family_name]]))
  
  purrr::map_dfr(required_vars, function(variable) {
    tibble(
      family = family_name,
      variable = variable,
      n_metric_frame = nrow(data),
      n_missing = sum(is.na(data[[variable]])),
      percent_missing = 100 * mean(is.na(data[[variable]]))
    )
  })
})

prognostic_excluded_rows <- purrr::imap_dfr(family_data, function(data, family_name) {
  required_vars <- unique(c(prognostic_common_vars, family_metric_vars[[family_name]]))
  excluded_indices <- which(!complete.cases(data[required_vars]))
  
  tibble(
    family = rep(family_name, length(excluded_indices)),
    source_row = data$.source_row[excluded_indices],
    hospitalization_id = data$hospitalization_id[excluded_indices],
    missing_variables = purrr::map_chr(excluded_indices, function(row_index) {
      missing_vars <- required_vars[purrr::map_lgl(
        required_vars,
        ~ is.na(data[[.x]][row_index])
      )]
      
      paste(missing_vars, collapse = "; ")
    })
  )
})

prognostic_sample_policy_checks <- list(
  prognostic_common_sample_covers_all_alternatives = all(
    normalization_adjustment_vars %in% prognostic_common_vars
  ),
  prognostic_exclusion_details_reconcile = all(purrr::map2_lgl(
    prognostic_sample_qc$family,
    prognostic_sample_qc$n_excluded,
    function(family_name, expected_n) {
      sum(prognostic_excluded_rows$family == family_name) == expected_n
    }
  )),
  prognostic_exclusion_reasons_nonempty = all(nzchar(prognostic_excluded_rows$missing_variables))
)

stopifnot(isTRUE(all(unlist(prognostic_sample_policy_checks, use.names = FALSE))))

normalization_analysis_diagnostics$prognostic_sample_configuration <- list(
  policy = paste0(
    "One common complete-case sample per mechanic family across all seven alternative normalization/size variables, ",
    "all candidate mechanic specifications, and both adjustment levels. Exclusions affect these outcome models only; ",
    "they do not change the upstream cohort or descriptive discordance/reclassification samples."
  ),
  required_common_variables = prognostic_common_vars,
  sample_summary = prognostic_sample_qc,
  variable_missingness = prognostic_sample_missingness_qc,
  excluded_rows = prognostic_excluded_rows,
  source_rows = purrr::map(
    prognostic_analysis_data,
    ~ .x %>% select(.source_row, hospitalization_id)
  )
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  prognostic_sample_policy_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

fit_prog <- function(family_name, spec, exposure, adjustment_id, normalization_spec) {
  analysis_data <- prognostic_analysis_data[[family_name]]
  
  registry_row <- normalization_adjustment_specs %>%
    filter(.data$normalization_spec == .env$normalization_spec)
  
  if (nrow(registry_row) != 1L) {
    stop("Expected exactly one normalization specification: ", normalization_spec)
  }
  
  covariates <- get_normalization_adjustment(normalization_spec, adjustment_id)
  model_formula <- as.formula(paste("mortality_event_60 ~", exposure, "+", covariates))
  
  fitted_model <- glm(
    model_formula, data = analysis_data, family = binomial, x = TRUE, y = TRUE, model = TRUE
  )
  
  model_frame_rows <- rownames(model.frame(fitted_model))
  analysis_data_rows <- rownames(analysis_data)
  
  if (!identical(model_frame_rows, analysis_data_rows)) {
    stop("The fitted prognostic model did not retain the expected analysis-data rows.")
  }
  
  list(
    meta = tibble(
      family = family_name,
      spec = spec,
      normalization_spec = normalization_spec,
      normalization_variable = registry_row$normalization_variable[[1]],
      normalization_label = registry_row$normalization_label[[1]],
      normalization_kind = registry_row$normalization_kind[[1]],
      adjustment = adjustment_id,
      adjusted = unname(adjustment_labels[[adjustment_id]]),
      formula = paste(deparse(formula(fitted_model)), collapse = " "),
      aic = AIC(fitted_model),
      auc = auc_from_model(fitted_model),
      n = nobs(fitted_model),
      n_events = sum(model.response(model.frame(fitted_model)) == 1L),
      n_nonevents = sum(model.response(model.frame(fitted_model)) == 0L)
    ),
    model = fitted_model
  )
}

# =============================================================================
# Fit and diagnose candidate prognostic models
# =============================================================================

prognostic_fit_output_path <- file.path(table_dir, paste0("norm_prognostic_fit_", site_name, ".csv"))

prognostic_fit_grid <- crossing(
  prog_specs,
  normalization_spec = normalization_adjustment_specs$normalization_spec,
  adjustment = names(adjustment_labels)
)

message("Script 05 | fitting candidate prognostic models")
prog_fits <- pmap(
  prognostic_fit_grid,
  function(family, spec, exposure, normalization_spec, adjustment) {
    fit_prog(
      family_name = family,
      spec = spec,
      exposure = exposure,
      adjustment_id = adjustment,
      normalization_spec = normalization_spec
    )
  }
)

names(prog_fits) <- paste(
  prognostic_fit_grid$family,
  prognostic_fit_grid$normalization_spec,
  prognostic_fit_grid$adjustment,
  prognostic_fit_grid$spec,
  sep = " | "
)

prognostic_fit_grid_checks <- list(
  prognostic_grid_has_expected_number_of_models = nrow(prognostic_fit_grid) ==
    nrow(prog_specs) * nrow(normalization_adjustment_specs) * length(adjustment_labels),
  prognostic_grid_keys_unique = !anyDuplicated(
    prognostic_fit_grid[c("family", "normalization_spec", "adjustment", "spec")]
  ),
  prognostic_fit_names_unique = !anyDuplicated(names(prog_fits)),
  prognostic_fit_count_matches_grid = length(prog_fits) == nrow(prognostic_fit_grid),
  prognostic_fit_metadata_matches_grid = all(purrr::map2_lgl(
    prog_fits,
    seq_len(nrow(prognostic_fit_grid)),
    function(fit, fit_index) {
      identical(fit$meta$family[[1]], prognostic_fit_grid$family[[fit_index]]) &&
        identical(fit$meta$spec[[1]], prognostic_fit_grid$spec[[fit_index]]) &&
        identical(fit$meta$normalization_spec[[1]], prognostic_fit_grid$normalization_spec[[fit_index]]) &&
        identical(fit$meta$adjustment[[1]], prognostic_fit_grid$adjustment[[fit_index]])
    }
  ))
)

stopifnot(isTRUE(all(unlist(prognostic_fit_grid_checks, use.names = FALSE))))

normalization_analysis_diagnostics$prognostic_sample_configuration$fit_grid <- prognostic_fit_grid

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  prognostic_fit_grid_checks
)

prognostic_model_qc <- purrr::map_dfr(prog_fits, function(fit) {
  fitted_model <- fit$model
  model_matrix <- model.matrix(fitted_model)
  coefficient_table <- summary(fitted_model)$coefficients
  aliased_coefficients <- is.na(coef(fitted_model))
  outcome <- model.response(model.frame(fitted_model))
  n_events <- sum(outcome == 1L)
  n_nonevents <- sum(outcome == 0L)
  n_parameters_excluding_intercept <- fitted_model$rank - 1L
  
  fit$meta %>%
    transmute(
      family, spec,
      normalization_spec, normalization_variable, normalization_label, normalization_kind,
      adjustment, adjusted, formula, n_obs = n,
      n_events, n_nonevents,
      n_parameters_excluding_intercept = n_parameters_excluding_intercept,
      events_per_parameter = n_events / n_parameters_excluding_intercept,
      converged = isTRUE(fitted_model$converged),
      boundary = isTRUE(fitted_model$boundary),
      iterations = fitted_model$iter,
      model_matrix_columns = ncol(model_matrix),
      model_matrix_rank = qr(model_matrix)$rank,
      full_rank = qr(model_matrix)$rank == ncol(model_matrix),
      n_aliased_coefficients = sum(aliased_coefficients),
      coefficients_finite = all(is.finite(coef(fitted_model)[!aliased_coefficients])),
      standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
      minimum_fitted_probability = min(fitted(fitted_model)),
      maximum_fitted_probability = max(fitted(fitted_model)),
      n_extreme_fitted_probabilities = sum(
        fitted(fitted_model) < EXTREME_FITTED_PROBABILITY_THRESHOLD |
          fitted(fitted_model) > 1 - EXTREME_FITTED_PROBABILITY_THRESHOLD
      ),
      aic, auc
    )
})

MIN_PROGNOSTIC_EVENTS_PER_PARAMETER <- 10
HIGHLY_SPARSE_PROGNOSTIC_EVENTS_PER_PARAMETER <- 5

prognostic_model_qc <- prognostic_model_qc %>%
  mutate(
    computationally_usable = converged & !boundary & full_rank & n_aliased_coefficients == 0L &
      coefficients_finite & standard_errors_finite & n_extreme_fitted_probabilities == 0L,
    sparse_events = events_per_parameter < MIN_PROGNOSTIC_EVENTS_PER_PARAMETER,
    highly_sparse_events = events_per_parameter < HIGHLY_SPARSE_PROGNOSTIC_EVENTS_PER_PARAMETER,
    comparison_status = case_when(
      !computationally_usable ~ "Suppressed: computational instability or separation",
      highly_sparse_events ~ "Exploratory: highly sparse events",
      sparse_events ~ "Exploratory: sparse events",
      TRUE ~ "Eligible for conventional comparison"
    ),
    comparison_usable = computationally_usable,
    interpretation = case_when(
      !converged ~ "Do not interpret: model did not converge.",
      boundary ~ "Do not interpret: model fit is on the boundary.",
      !full_rank | n_aliased_coefficients > 0L ~ "Do not interpret: model matrix is rank deficient.",
      !coefficients_finite | !standard_errors_finite ~ "Do not interpret: model inference is nonfinite.",
      n_extreme_fitted_probabilities > 0L ~
        "Do not interpret: extreme fitted probabilities indicate possible separation.",
      events_per_parameter < HIGHLY_SPARSE_PROGNOSTIC_EVENTS_PER_PARAMETER ~
        "Exploratory only: fewer than 5 deaths per fitted parameter.",
      events_per_parameter < MIN_PROGNOSTIC_EVENTS_PER_PARAMETER ~
        "Interpret cautiously: fewer than 10 deaths per fitted parameter.",
      TRUE ~ "No prespecified sparse-event or computational warning."
    )
  )

prognostic_same_rows_qc <- prognostic_model_qc %>%
  distinct(family, normalization_spec, adjustment) %>%
  mutate(
    same_rows = purrr::pmap_lgl(
      list(family, normalization_spec, adjustment),
      function(family_name, normalization_id, adjustment_id) {
        indices <- which(
          prognostic_fit_grid$family == family_name &
            prognostic_fit_grid$normalization_spec == normalization_id &
            prognostic_fit_grid$adjustment == adjustment_id
        )
        
        if (length(indices) != 5L) {
          stop(
            "Expected five prognostic models for ", family_name, " | ",
            normalization_id, " | ", adjustment_id, "; found ", length(indices), "."
          )
        }
        
        reference_rows <- rownames(model.frame(prog_fits[[indices[[1]]]]$model))
        
        all(purrr::map_lgl(
          prog_fits[indices],
          ~ identical(rownames(model.frame(.x$model)), reference_rows)
        ))
      }
    )
  )

# The common-sample policy also requires identical retained rows across alternatives and adjustment levels.
prognostic_cross_specification_rows_qc <- tibble(family = names(prognostic_analysis_data)) %>%
  mutate(
    same_rows_across_normalizations_and_adjustments = purrr::map_lgl(family, function(family_name) {
      indices <- which(prognostic_fit_grid$family == family_name)
      expected_rows <- rownames(prognostic_analysis_data[[family_name]])
      
      if (length(indices) != 5L * nrow(normalization_adjustment_specs) * length(adjustment_labels)) {
        stop("Unexpected number of prognostic models for mechanic family: ", family_name)
      }
      
      all(purrr::map_lgl(
        prog_fits[indices],
        ~ identical(rownames(model.frame(.x$model)), expected_rows)
      ))
    })
  )

prognostic_family_status <- prognostic_model_qc %>%
  group_by(
    family, normalization_spec, normalization_variable, normalization_label, normalization_kind,
    adjustment, adjusted
  ) %>%  
  summarise(
    n_models = n(), n_obs = first(n_obs), n_events = first(n_events),
    minimum_events_per_parameter = min(events_per_parameter),
    all_models_converged = all(converged), any_model_at_boundary = any(boundary),
    all_model_matrices_full_rank = all(full_rank),
    all_coefficients_finite = all(coefficients_finite),
    all_standard_errors_finite = all(standard_errors_finite),
    maximum_extreme_fitted_probabilities = max(n_extreme_fitted_probabilities),
    all_models_computationally_usable = all(computationally_usable),
    all_models_conventionally_eligible = all(
      comparison_status == "Eligible for conventional comparison"
    ),
    family_comparison_status = case_when(
      !all(computationally_usable) ~ "Suppressed: at least one model is computationally unstable",
      any(highly_sparse_events) ~ "Exploratory: highly sparse events",
      any(sparse_events) ~ "Exploratory: sparse events",
      TRUE ~ "Eligible for conventional comparison"
    ),
    .groups = "drop"
  )

prognostic_model_integrity_checks <- list(
  expected_number_of_models_fitted = length(prog_fits) ==
    nrow(prog_specs) * nrow(normalization_adjustment_specs) * length(adjustment_labels),
  prognostic_qc_row_count_matches_fitted_models = nrow(prognostic_model_qc) == length(prog_fits),
  prognostic_qc_model_keys_unique = !anyDuplicated(
    prognostic_model_qc[c("family", "normalization_spec", "adjustment", "spec")]
  ),
  all_requested_normalizations_present = setequal(
    unique(prognostic_model_qc$normalization_spec),
    normalization_adjustment_specs$normalization_spec
  ),
  five_models_per_family_normalization_and_adjustment = prognostic_model_qc %>%
    count(family, normalization_spec, adjustment) %>%
    summarise(
      valid = n() == length(family_data) * nrow(normalization_adjustment_specs) * length(adjustment_labels) &&
        all(n == 5L)
    ) %>%
    pull(valid),
  prognostic_same_row_checks_have_expected_comparison_sets = nrow(prognostic_same_rows_qc) ==
    length(family_data) * nrow(normalization_adjustment_specs) * length(adjustment_labels),
  all_models_use_expected_family_n = purrr::map2_lgl(
    prognostic_model_qc$family,
    prognostic_model_qc$n_obs,
    ~ .y == nrow(prognostic_analysis_data[[.x]])
  ) %>%
    all(),
  all_models_have_expected_event_counts = purrr::map2_lgl(
    prognostic_model_qc$family,
    prognostic_model_qc$n_events,
    ~ .y == sum(prognostic_analysis_data[[.x]]$mortality_event_60 == 1L)
  ) %>%
    all(),
  all_models_have_expected_nonevent_counts = purrr::map2_lgl(
    prognostic_model_qc$family,
    prognostic_model_qc$n_nonevents,
    ~ .y == sum(prognostic_analysis_data[[.x]]$mortality_event_60 == 0L)
  ) %>%
    all(),
  all_comparison_sets_use_same_rows = all(prognostic_same_rows_qc$same_rows),
  all_normalizations_and_adjustments_use_common_family_rows = all(
    prognostic_cross_specification_rows_qc$same_rows_across_normalizations_and_adjustments
  ),
  model_formulas_contain_the_expected_normalization_term = purrr::map2_lgl(
    prognostic_model_qc$formula,
    prognostic_model_qc$normalization_variable,
    function(model_formula, expected_variable) {
      expected_variable %in% all.vars(as.formula(model_formula))
    }
  ) %>%
    all(),
  all_model_matrices_full_rank = all(prognostic_model_qc$full_rank),
  no_aliased_coefficients = all(prognostic_model_qc$n_aliased_coefficients == 0L),
  all_coefficients_finite = all(prognostic_model_qc$coefficients_finite),
  all_standard_errors_finite = all(prognostic_model_qc$standard_errors_finite),
  all_fitted_probability_ranges_valid = all(
    is.finite(prognostic_model_qc$minimum_fitted_probability) &
      is.finite(prognostic_model_qc$maximum_fitted_probability) &
      prognostic_model_qc$minimum_fitted_probability >= 0 &
      prognostic_model_qc$maximum_fitted_probability <= 1
  ),
  all_aic_values_finite = all(is.finite(prognostic_model_qc$aic)),
  all_auc_values_valid = all(
    is.finite(prognostic_model_qc$auc) &
      prognostic_model_qc$auc >= 0 &
      prognostic_model_qc$auc <= 1
  )
)

if (!isTRUE(all(unlist(prognostic_model_integrity_checks, use.names = FALSE)))) {
  print(prognostic_model_qc, n = Inf, width = Inf)
  print(prognostic_same_rows_qc, n = Inf, width = Inf)
  print(prognostic_cross_specification_rows_qc, n = Inf, width = Inf)
  print(prognostic_family_status, n = Inf, width = Inf)
  print(prognostic_model_integrity_checks)
}

stopifnot(isTRUE(all(unlist(prognostic_model_integrity_checks, use.names = FALSE))))

normalization_analysis_diagnostics$prognostic_fit_validation <- list(
  fit_grid = prognostic_fit_grid,
  model_diagnostics = prognostic_model_qc,
  within_comparison_set_row_validation = prognostic_same_rows_qc,
  across_normalizations_and_adjustments_row_validation = prognostic_cross_specification_rows_qc,
  comparison_set_status = prognostic_family_status,
  checks = prognostic_model_integrity_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

prognostic_influence_qc <- purrr::map2_dfr(
  prog_fits,
  seq_along(prog_fits),
  function(fit, fit_index) {
    fitted_model <- fit$model
    family_name <- prognostic_fit_grid$family[[fit_index]]
    spec_name <- prognostic_fit_grid$spec[[fit_index]]
    normalization_id <- prognostic_fit_grid$normalization_spec[[fit_index]]
    adjustment_id <- prognostic_fit_grid$adjustment[[fit_index]]
    analysis_data <- prognostic_analysis_data[[family_name]]
    
    registry_row <- normalization_adjustment_specs %>%
      filter(.data$normalization_spec == .env$normalization_id)
    
    stopifnot(nrow(registry_row) == 1L)
    model_frame_rows <- rownames(model.frame(fitted_model))
    analysis_data_rows <- rownames(analysis_data)
    model_row_index <- match(model_frame_rows, analysis_data_rows)
    
    if (any(is.na(model_row_index))) {
      stop("A prognostic influence row could not be matched to its source analysis row.")
    }
    
    n <- nobs(fitted_model)
    p <- fitted_model$rank
    cooks_values <- cooks.distance(fitted_model)
    leverage_values <- hatvalues(fitted_model)
    standardized_residuals <- rstandard(fitted_model)
    dfbeta_values <- dfbetas(fitted_model)
    maximum_absolute_dfbeta <- apply(abs(dfbeta_values), 1L, max, na.rm = TRUE)
    
    tibble(
      family = family_name,
      spec = spec_name,
      normalization_spec = normalization_id,
      normalization_variable = registry_row$normalization_variable[[1]],
      normalization_label = registry_row$normalization_label[[1]],
      normalization_kind = registry_row$normalization_kind[[1]],
      adjustment = adjustment_id,
      adjusted = unname(adjustment_labels[[adjustment_id]]),
      analysis_row = seq_len(n),
      source_row = analysis_data$.source_row[model_row_index],
      hospitalization_id = analysis_data$hospitalization_id[model_row_index],
      cooks_distance = unname(cooks_values),
      cooks_threshold = 4 / n,
      high_cooks_distance = cooks_distance > cooks_threshold,
      leverage = unname(leverage_values),
      leverage_threshold = 2 * p / n,
      high_leverage = leverage > leverage_threshold,
      standardized_residual = unname(standardized_residuals),
      large_standardized_residual = abs(standardized_residual) > 3,
      maximum_absolute_dfbeta = maximum_absolute_dfbeta,
      dfbeta_threshold = 2 / sqrt(n),
      high_dfbeta = maximum_absolute_dfbeta > dfbeta_threshold,
      flagged_by_any_measure = high_cooks_distance | high_leverage |
        large_standardized_residual | high_dfbeta
    )
  }
)

prognostic_influence_summary <- prognostic_influence_qc %>%
  group_by(
    family, spec, normalization_spec, normalization_variable, normalization_label, normalization_kind,
    adjustment, adjusted
  ) %>%
  summarise(
    n_obs = n(),
    n_high_cooks_distance = sum(high_cooks_distance),
    n_high_leverage = sum(high_leverage),
    n_large_standardized_residual = sum(large_standardized_residual),
    n_high_dfbeta = sum(high_dfbeta),
    n_flagged_by_any_measure = sum(flagged_by_any_measure),
    maximum_cooks_distance = max(cooks_distance),
    maximum_leverage = max(leverage),
    maximum_absolute_standardized_residual = max(abs(standardized_residual)),
    maximum_absolute_dfbeta = max(maximum_absolute_dfbeta),
    .groups = "drop"
  )

prognostic_influence_checks <- list(
  influence_diagnostics_have_expected_rows = nrow(prognostic_influence_qc) ==
    sum(prognostic_model_qc$n_obs),
  influence_source_rows_complete = !any(is.na(prognostic_influence_qc$source_row)),
  influence_hospitalization_ids_complete = !any(is.na(prognostic_influence_qc$hospitalization_id)),
  influence_analysis_rows_sequential_within_model = prognostic_influence_qc %>%
    group_by(family, spec, normalization_spec, adjustment) %>%
    summarise(
      rows_sequential = identical(analysis_row, seq_len(n())),
      .groups = "drop"
    ) %>%
    pull(rows_sequential) %>%
    all(),
  influence_source_rows_valid = all(
    prognostic_influence_qc$source_row %in% seq_len(nrow(cross_sectional))
  ),
  influence_rows_unique_within_model = prognostic_influence_qc %>%
    group_by(family, spec, normalization_spec, adjustment) %>%
    summarise(rows_unique = n_distinct(source_row) == n(), .groups = "drop") %>%
    pull(rows_unique) %>%
    all(),
  influence_values_finite = all(is.finite(c(
    prognostic_influence_qc$cooks_distance,
    prognostic_influence_qc$leverage,
    prognostic_influence_qc$standardized_residual,
    prognostic_influence_qc$maximum_absolute_dfbeta
  )))
)

if (!isTRUE(all(unlist(prognostic_influence_checks, use.names = FALSE)))) {
  print(prognostic_influence_summary, n = Inf, width = Inf)
  print(prognostic_influence_checks)
}

stopifnot(isTRUE(all(unlist(prognostic_influence_checks, use.names = FALSE))))

prog_tbl <- prognostic_model_qc %>%
  group_by(family, normalization_spec, adjustment, adjusted) %>%
  mutate(
    family_computationally_usable = all(computationally_usable),
    aic_raw = aic,
    auc_raw = auc,
    aic = if_else(family_computationally_usable, aic_raw, NA_real_),
    auc = if_else(family_computationally_usable, auc_raw, NA_real_),
    delta_aic = if_else(family_computationally_usable, aic_raw - min(aic_raw), NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    site = site_name,
    spec = factor(spec, levels = unique(prog_specs$spec)),
    .before = 1
  ) %>%
  select(
    site, family, spec,
    normalization_spec, normalization_variable, normalization_label, normalization_kind,
    adjustment, adjusted, formula, aic, auc, delta_aic,
    n = n_obs, n_events, n_nonevents, events_per_parameter,
    computationally_usable, family_computationally_usable,
    sparse_events, highly_sparse_events, comparison_status, interpretation,
    aic_raw, auc_raw
  )

prognostic_suppression_checks <- list(
  computational_status_consistent_with_convergence = all(
    !prognostic_model_qc$computationally_usable | prognostic_model_qc$converged
  ),
  unstable_models_have_suppressed_comparison_values = all(
    prog_tbl$computationally_usable |
      (is.na(prog_tbl$aic) & is.na(prog_tbl$auc) & is.na(prog_tbl$delta_aic))
  ),
  stable_families_retain_comparison_values = all(
    !prog_tbl$family_computationally_usable |
      (is.finite(prog_tbl$aic) & is.finite(prog_tbl$auc) & is.finite(prog_tbl$delta_aic))
  ),
  family_status_agrees_with_model_status = prog_tbl %>%
    group_by(family, normalization_spec, adjustment) %>%
    summarise(
      stored_family_status = first(family_computationally_usable),
      calculated_family_status = all(computationally_usable),
      status_matches = stored_family_status == calculated_family_status,
      .groups = "drop"
    ) %>%
    pull(status_matches) %>%
    all()
)

prognostic_suppression_qc <- prog_tbl %>%
  group_by(
    family, normalization_spec, normalization_variable, normalization_label, normalization_kind,
    adjustment, adjusted
  ) %>%
  summarise(
    n_models = n(),
    family_computationally_usable = first(family_computationally_usable),
    n_computationally_usable_models = sum(computationally_usable),
    n_aic_values_exported = sum(!is.na(aic)),
    n_auc_values_exported = sum(!is.na(auc)),
    n_delta_aic_values_exported = sum(!is.na(delta_aic)),
    all_raw_aic_values_finite = all(is.finite(aic_raw)),
    all_raw_auc_values_valid = all(is.finite(auc_raw) & auc_raw >= 0 & auc_raw <= 1),
    .groups = "drop"
  )

if (!isTRUE(all(unlist(prognostic_suppression_checks, use.names = FALSE)))) {
  print(prognostic_suppression_qc, n = Inf, width = Inf)
  print(prognostic_suppression_checks)
}

stopifnot(isTRUE(all(unlist(prognostic_suppression_checks, use.names = FALSE))))

write_csv(prog_tbl, prognostic_fit_output_path)

prognostic_fit_output_qc <- tibble(
  output = "Prognostic fit comparison",
  path = prognostic_fit_output_path,
  exists = file.exists(prognostic_fit_output_path),
  size_bytes = if (file.exists(prognostic_fit_output_path)) {
    file.info(prognostic_fit_output_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

prognostic_fit_output_checks <- list(
  prognostic_fit_output_exists = prognostic_fit_output_qc$exists,
  prognostic_fit_output_nonempty = prognostic_fit_output_qc$nonempty
)

stopifnot(isTRUE(all(unlist(prognostic_fit_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$prognostic_utility <- list(
  outcome = "60-day mortality",
  statistical_limitation = paste0(
    "The cross-sectional cohort contains ",
    sum(cross_sectional$mortality_event_60 == 1L),
    " deaths. Events per fitted parameter are reported for every model. Computational convergence and full rank do not ",
    "remove sparse-event bias or imprecision."
  ),
  adjustment_specs = adjustment_specs,
  adjustment_labels = adjustment_labels,
  normalization_specifications = normalization_adjustment_specs,
  normalization_covariate_grid = normalization_covariate_grid,
  adjustment_specs_note = paste0(
    "adjustment_specs retains the original McLaren reference strings. ",
    "normalization_covariate_grid records the actual adjustment strings used for all seven alternatives."
  ),
  common_sample_configuration = normalization_analysis_diagnostics$prognostic_sample_configuration,
  candidate_specifications = prog_specs,
  comparison_thresholds = list(
    minimum_events_per_parameter_for_conventional_interpretation =
      MIN_PROGNOSTIC_EVENTS_PER_PARAMETER,
    highly_sparse_events_per_parameter_threshold =
      HIGHLY_SPARSE_PROGNOSTIC_EVENTS_PER_PARAMETER,
    note = paste0(
      "These thresholds are interpretation heuristics rather than absolute validity boundaries. ",
      "Computationally stable models below the thresholds are retained as exploratory. ",
      "AIC, AUC, and delta AIC are suppressed within a five-model comparison set if any model in that set fails ",
      "the existing computational-use criteria. Each set is defined by mechanic family, normalization/size ",
      "specification, and adjustment level."
    )
  ),
  sample_summary = prognostic_sample_qc,
  source_rows = purrr::map(
    prognostic_analysis_data,
    ~ .x %>% select(.source_row, hospitalization_id)
  ),
  model_diagnostics = prognostic_model_qc,
  influence_diagnostics = list(
    summary = prognostic_influence_summary,
    flagged_observations = prognostic_influence_qc %>%
      filter(flagged_by_any_measure),
    thresholds = tibble(
      diagnostic = c("Cook's distance", "Leverage", "Absolute standardized residual", "Absolute DFBETA"),
      threshold = c("4 / n", "2 × model rank / n", "3", "2 / sqrt(n)")
    )
  ),
  comparison_row_validation = prognostic_same_rows_qc,
  cross_specification_row_validation = prognostic_cross_specification_rows_qc,
  comparison_set_definition = "Mechanic family × normalization/size specification × adjustment level",
  family_status = prognostic_family_status,
  fit_comparison = prog_tbl,
  suppression_diagnostics = prognostic_suppression_qc,
  output = prognostic_fit_output_qc
)

normalization_analysis_diagnostics$warnings$prognostic_models <- prognostic_model_qc %>%
  filter(interpretation != "No prespecified sparse-event or computational warning.")

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  prognostic_sample_checks,
  prognostic_model_integrity_checks,
  prognostic_influence_checks,
  prognostic_suppression_checks,
  prognostic_fit_output_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# Prognostic coefficient estimates
# =============================================================================
# Exposure coefficients are retained on the raw log-odds scale in the local diagnostic RDS. Odds ratios are exported
# only for computationally usable models. Computationally stable but sparse models remain explicitly exploratory.

prognostic_coefficient_output_path <- file.path(
  table_dir, paste0("norm_prognostic_coefs_", site_name, ".csv")
)

z975 <- qnorm(0.975)

prognostic_coefficient_diagnostics <- purrr::map2_dfr(
  prog_fits,
  seq_along(prog_fits),
  function(fit, fit_index) {
    fitted_model <- fit$model
    family_name <- prognostic_fit_grid$family[[fit_index]]
    spec_name <- prognostic_fit_grid$spec[[fit_index]]
    normalization_id <- prognostic_fit_grid$normalization_spec[[fit_index]]
    adjustment_id <- prognostic_fit_grid$adjustment[[fit_index]]
    analysis_data <- prognostic_analysis_data[[family_name]]
    
    model_status <- prognostic_model_qc %>%
      filter(
        .data$family == .env$family_name,
        .data$spec == .env$spec_name,
        .data$normalization_spec == .env$normalization_id,
        .data$adjustment == .env$adjustment_id
      )
    
    stopifnot(nrow(model_status) == 1L)
    
    broom::tidy(fitted_model) %>%
      filter(str_detect(term, "^log\\(")) %>%
      mutate(
        exposure_variable = str_match(term, "^log\\(([^)]+)\\)$")[, 2],
        exposure_sd = purrr::map_dbl(exposure_variable, function(variable) {
          if (is.na(variable) || !variable %in% names(analysis_data)) {
            stop("Could not identify the source variable for prognostic term: ", term[[1]])
          }
          
          sd(log(analysis_data[[variable]]))
        }),
        family = family_name,
        spec = spec_name,
        normalization_spec = normalization_id,
        normalization_variable = model_status$normalization_variable[[1]],
        normalization_label = model_status$normalization_label[[1]],
        normalization_kind = model_status$normalization_kind[[1]],
        adjustment = adjustment_id,
        adjusted = unname(adjustment_labels[[adjustment_id]]),
        n = nobs(fitted_model),
        n_events = sum(model.response(model.frame(fitted_model)) == 1L),
        computationally_usable = model_status$computationally_usable[[1]],
        comparison_status = model_status$comparison_status[[1]],
        interpretation = model_status$interpretation[[1]],
        .before = 1
      ) %>%
      transmute(
        family, spec,
        normalization_spec, normalization_variable, normalization_label, normalization_kind,
        adjustment, adjusted, term, exposure_variable, exposure_sd,
        log_odds_coefficient = estimate,
        std_error = std.error,
        statistic,
        p_value = p.value,
        n,
        n_events,
        computationally_usable,
        comparison_status,
        interpretation
      )
  }
)

prognostic_coefficient_diagnostics <- prognostic_coefficient_diagnostics %>%
  mutate(
    inferential_values_exported = computationally_usable,
    or_per_log = if_else(computationally_usable, exp(log_odds_coefficient), NA_real_),
    or_per_log_lo = if_else(
      computationally_usable,
      exp(log_odds_coefficient - z975 * std_error),
      NA_real_
    ),
    or_per_log_hi = if_else(
      computationally_usable,
      exp(log_odds_coefficient + z975 * std_error),
      NA_real_
    ),
    or_per_sd = if_else(
      computationally_usable,
      exp(log_odds_coefficient * exposure_sd),
      NA_real_
    ),
    or_per_sd_lo = if_else(
      computationally_usable,
      exp((log_odds_coefficient - z975 * std_error) * exposure_sd),
      NA_real_
    ),
    or_per_sd_hi = if_else(
      computationally_usable,
      exp((log_odds_coefficient + z975 * std_error) * exposure_sd),
      NA_real_
    )
  )

prog_coefs <- prognostic_coefficient_diagnostics %>%
  transmute(
    site = site_name,
    family,
    spec,
    normalization_spec,
    normalization_variable,
    normalization_label,
    normalization_kind,
    adjustment,
    adjusted,
    term,
    exposure_variable,
    sd_log = exposure_sd,
    or_per_log,
    or_per_log_lo,
    or_per_log_hi,
    or_per_sd,
    or_per_sd_lo,
    or_per_sd_hi,
    std_error = if_else(computationally_usable, std_error, NA_real_),
    p_value = if_else(computationally_usable, p_value, NA_real_),
    n,
    n_events,
    computationally_usable,
    comparison_status,
    interpretation
  )

exposure_terms_per_spec <- purrr::map_int(prog_specs$exposure, function(exposure_spec) {
  length(all.vars(as.formula(paste("~", exposure_spec))))
})

expected_exposure_coefficient_rows <- n_distinct(prognostic_fit_grid$adjustment) *
  n_distinct(prognostic_fit_grid$normalization_spec) * sum(exposure_terms_per_spec)

prognostic_coefficient_checks <- list(
  coefficient_results_nonempty = nrow(prognostic_coefficient_diagnostics) > 0L,
  coefficient_model_term_keys_unique = !anyDuplicated(
    prognostic_coefficient_diagnostics[c("family", "normalization_spec", "adjustment", "spec", "term")]
  ),
  coefficient_export_keys_unique = !anyDuplicated(
    prog_coefs[c("site", "family", "normalization_spec", "adjustment", "spec", "term")]
  ),
  all_requested_normalizations_have_coefficient_results = setequal(
    unique(prognostic_coefficient_diagnostics$normalization_spec),
    normalization_adjustment_specs$normalization_spec
  ),
  coefficient_results_cover_every_fitted_model = nrow(distinct(
    prognostic_coefficient_diagnostics,
    family, normalization_spec, adjustment, spec
  )) == nrow(prognostic_fit_grid),
  coefficient_counts_match_each_mechanic_specification = all(purrr::pmap_lgl(
    prognostic_fit_grid,
    function(family, spec, exposure, normalization_spec, adjustment) {
      expected_terms <- all.vars(as.formula(paste("~", exposure)))
      
      observed_variables <- prognostic_coefficient_diagnostics %>%
        filter(
          .data$family == .env$family,
          .data$spec == .env$spec,
          .data$normalization_spec == .env$normalization_spec,
          .data$adjustment == .env$adjustment
        ) %>%
        pull(exposure_variable)
      
      length(observed_variables) == length(expected_terms) &&
        setequal(observed_variables, expected_terms)
    }
  )),
  coefficient_sample_sizes_match_common_family_samples = all(purrr::map2_lgl(
    prognostic_coefficient_diagnostics$family,
    prognostic_coefficient_diagnostics$n,
    ~ .y == nrow(prognostic_analysis_data[[.x]])
  )),
  coefficient_event_counts_match_common_family_samples = all(purrr::map2_lgl(
    prognostic_coefficient_diagnostics$family,
    prognostic_coefficient_diagnostics$n_events,
    ~ .y == sum(prognostic_analysis_data[[.x]]$mortality_event_60 == 1L)
  )),
  exposure_terms_per_spec_valid = identical(exposure_terms_per_spec, c(1L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L)),
  expected_number_of_exposure_coefficient_rows = nrow(prognostic_coefficient_diagnostics) ==
    expected_exposure_coefficient_rows,
  exposure_variables_identified = all(
    !is.na(prognostic_coefficient_diagnostics$exposure_variable) &
      nzchar(prognostic_coefficient_diagnostics$exposure_variable)
  ),
  exposure_standard_deviations_finite = all(
    is.finite(prognostic_coefficient_diagnostics$exposure_sd)
  ),
  exposure_standard_deviations_positive = all(
    prognostic_coefficient_diagnostics$exposure_sd > 0
  ),
  raw_log_odds_coefficients_finite = all(
    is.finite(prognostic_coefficient_diagnostics$log_odds_coefficient)
  ),
  raw_standard_errors_finite = all(
    is.finite(prognostic_coefficient_diagnostics$std_error)
  ),
  computationally_usable_estimates_exported = check_subset(
    prog_coefs,
    prog_coefs$computationally_usable,
    function(data) {
      all(
        is.finite(data$or_per_log) &
          is.finite(data$or_per_log_lo) &
          is.finite(data$or_per_log_hi) &
          is.finite(data$or_per_sd) &
          is.finite(data$or_per_sd_lo) &
          is.finite(data$or_per_sd_hi) &
          is.finite(data$std_error) &
          is.finite(data$p_value)
      )
    }
  ),
  computationally_unusable_estimates_suppressed = check_subset(
    prog_coefs,
    !prog_coefs$computationally_usable,
    function(data) {
      all(
        is.na(data$or_per_log) &
          is.na(data$or_per_log_lo) &
          is.na(data$or_per_log_hi) &
          is.na(data$or_per_sd) &
          is.na(data$or_per_sd_lo) &
          is.na(data$or_per_sd_hi) &
          is.na(data$std_error) &
          is.na(data$p_value)
      )
    }
  ),
  exported_p_values_valid = all(
    is.na(prog_coefs$p_value) |
      (prog_coefs$p_value >= 0 & prog_coefs$p_value <= 1)
  ),
  coefficient_sample_sizes_positive = all(prog_coefs$n > 0L),
  coefficient_event_counts_positive = all(prog_coefs$n_events > 0L)
)

if (!isTRUE(all(unlist(prognostic_coefficient_checks, use.names = FALSE)))) {
  print(prognostic_coefficient_diagnostics, n = Inf, width = Inf)
  print(prog_coefs, n = Inf, width = Inf)
  print(prognostic_coefficient_checks)
}

stopifnot(isTRUE(all(unlist(prognostic_coefficient_checks, use.names = FALSE))))

write_csv(prog_coefs, prognostic_coefficient_output_path)

prognostic_coefficient_output_qc <- tibble(
  output = "Prognostic coefficient estimates",
  path = prognostic_coefficient_output_path,
  exists = file.exists(prognostic_coefficient_output_path),
  size_bytes = if (file.exists(prognostic_coefficient_output_path)) {
    file.info(prognostic_coefficient_output_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

prognostic_coefficient_output_checks <- list(
  prognostic_coefficient_output_exists = prognostic_coefficient_output_qc$exists,
  prognostic_coefficient_output_nonempty = prognostic_coefficient_output_qc$nonempty
)

stopifnot(isTRUE(all(unlist(prognostic_coefficient_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$prognostic_utility$coefficient_diagnostics <-
  prognostic_coefficient_diagnostics

normalization_analysis_diagnostics$prognostic_utility$coefficient_export <- list(
  data = prog_coefs,
  output = prognostic_coefficient_output_qc,
  expected_rows = expected_exposure_coefficient_rows,
  key_columns = c("site", "family", "normalization_spec", "adjustment", "spec", "term"),
  interpretation = paste0(
    "Exported coefficients describe the log-transformed mechanic exposure conditional on the selected alternative ",
    "normalization/size adjustment variable and the remaining model covariates. ",
    "OR per log-unit corresponds to an e-fold increase in the mechanic exposure; OR per SD uses the sample SD of its ",
    "log-transformed values in the common mechanic-family sample. ",
    "For separate mechanic-plus-size models, each coefficient is conditional on the other component. ",
    "The alternative adjustment-variable coefficients are not included in this mechanic-exposure export."
  )
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  prognostic_coefficient_checks,
  prognostic_coefficient_output_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# Prognostic-fit summary
# =============================================================================

# =============================================================================
# Form-versus-physiology decomposition
# =============================================================================
# Within each size reference, the locked model is nested in the corresponding separate model. For example,
# log(Ers x PBW) is the constrained form of log(Ers) + log(PBW), where both component coefficients are equal.
# The one-degree-of-freedom LRT tests whether freeing those coefficients improves fit.
#
# PBW-versus-PFVC comparisons within the same form are nonnested and are compared using AIC. Positive delta AIC values
# below indicate that the PFVC specification has lower AIC. All comparisons use the same family-specific patient sample.
# Computationally unstable families are retained as suppressed rows rather than assigned misleading comparison results.

form_physiology_output_path <- file.path(table_dir, paste0("norm_form_vs_physiology_", site_name, ".csv"))

get_prog_model <- function(family_name, spec_name, adjustment_id, normalization_spec) {
  model_index <- which(
    prognostic_fit_grid$family == family_name &
      prognostic_fit_grid$spec == spec_name &
      prognostic_fit_grid$adjustment == adjustment_id &
      prognostic_fit_grid$normalization_spec == normalization_spec
  )
  
  if (length(model_index) != 1L) {
    stop(
      "Expected exactly one prognostic model for ", family_name, " | ", normalization_spec,
      " | ", adjustment_id, " | ", spec_name, "; found ", length(model_index), "."
    )
  }
  
  prog_fits[[model_index]]$model
}

decomposition_grid <- crossing(
  family = names(family_data),
  normalization_spec = normalization_adjustment_specs$normalization_spec,
  adjustment = names(adjustment_labels)
)

decomp_tbl <- pmap_dfr(decomposition_grid, function(family, normalization_spec, adjustment) {
  specification_names <- if (family == "Elastance") {
    c(
      pbw_lock = "PBW-locked (Goligher)",
      pfvc_lock = "PFVC-locked",
      pbw_sep = "Separate (Ers + PBW)",
      pfvc_sep = "Separate (Ers + PFVC)"
    )
  } else {
    c(
      pbw_lock = "PBW-locked (Gattinoni)",
      pfvc_lock = "PFVC-locked",
      pbw_sep = "Separate (MP + PBW)",
      pfvc_sep = "Separate (MP + PFVC)"
    )
  }
  
  models <- purrr::map(
    specification_names,
    ~ get_prog_model(family, .x, adjustment, normalization_spec)
  )
  
  registry_row <- normalization_adjustment_specs %>%
    filter(.data$normalization_spec == .env$normalization_spec)
  
  model_status <- prognostic_model_qc %>%
    filter(
      .data$family == .env$family,
      .data$normalization_spec == .env$normalization_spec,
      .data$adjustment == .env$adjustment,
      .data$spec %in% unname(specification_names)
    )
  
  family_status <- prognostic_family_status %>%
    filter(
      .data$family == .env$family,
      .data$normalization_spec == .env$normalization_spec,
      .data$adjustment == .env$adjustment
    )
  
  stopifnot(nrow(registry_row) == 1L)
  
  stopifnot(
    length(models) == 4L,
    nrow(model_status) == 4L,
    nrow(family_status) == 1L
  )
  
  same_n <- n_distinct(purrr::map_int(models, nobs)) == 1L
  same_rows <- all(purrr::map_lgl(
    models,
    ~ identical(rownames(model.frame(.x)), rownames(model.frame(models[[1]])))
  ))
  
  computationally_usable <- all(model_status$computationally_usable) && same_n && same_rows
  
  if (!computationally_usable) {
    return(
      tibble(
        site = site_name, family = family,
        normalization_spec = normalization_spec,
        normalization_variable = registry_row$normalization_variable[[1]],
        normalization_label = registry_row$normalization_label[[1]],
        normalization_kind = registry_row$normalization_kind[[1]],
        adjustment = adjustment,
        adjusted = unname(adjustment_labels[[adjustment]]),
        n = nobs(models[[1]]),
        n_events = sum(model.response(model.frame(models[[1]])) == 1L),
        computationally_usable = FALSE,
        comparison_status = "Suppressed: computational instability or separation",
        form_pbw_lrt_df = NA_real_, form_pbw_lrt_p = NA_real_,
        form_pfvc_lrt_df = NA_real_, form_pfvc_lrt_p = NA_real_,
        phys_locked_dAIC = NA_real_, phys_separate_dAIC = NA_real_,
        aic_pbw_sep = NA_real_, aic_pfvc_sep = NA_real_
      )
    )
  }
  
  form_pbw <- anova(models$pbw_lock, models$pbw_sep, test = "LRT")
  form_pfvc <- anova(models$pfvc_lock, models$pfvc_sep, test = "LRT")
  
  comparison_status <- if (any(model_status$highly_sparse_events)) {
    "Exploratory: highly sparse events"
  } else if (any(model_status$sparse_events)) {
    "Exploratory: sparse events"
  } else {
    "Eligible for conventional comparison"
  }
  
  tibble(
    site = site_name, family = family,
    normalization_spec = normalization_spec,
    normalization_variable = registry_row$normalization_variable[[1]],
    normalization_label = registry_row$normalization_label[[1]],
    normalization_kind = registry_row$normalization_kind[[1]],
    adjustment = adjustment,
    adjusted = unname(adjustment_labels[[adjustment]]),
    n = nobs(models[[1]]),
    n_events = sum(model.response(model.frame(models[[1]])) == 1L),
    computationally_usable = TRUE,
    comparison_status = comparison_status,
    form_pbw_lrt_df = form_pbw$Df[2],
    form_pbw_lrt_p = form_pbw$`Pr(>Chi)`[2],
    form_pfvc_lrt_df = form_pfvc$Df[2],
    form_pfvc_lrt_p = form_pfvc$`Pr(>Chi)`[2],
    phys_locked_dAIC = AIC(models$pbw_lock) - AIC(models$pfvc_lock),
    phys_separate_dAIC = AIC(models$pbw_sep) - AIC(models$pfvc_sep),
    aic_pbw_sep = AIC(models$pbw_sep),
    aic_pfvc_sep = AIC(models$pfvc_sep)
  )
})

decomposition_model_row_qc <- pmap_dfr(decomposition_grid, function(family, normalization_spec, adjustment) {
  specification_names <- if (family == "Elastance") {
    c(
      pbw_lock = "PBW-locked (Goligher)",
      pfvc_lock = "PFVC-locked",
      pbw_sep = "Separate (Ers + PBW)",
      pfvc_sep = "Separate (Ers + PFVC)"
    )
  } else {
    c(
      pbw_lock = "PBW-locked (Gattinoni)",
      pfvc_lock = "PFVC-locked",
      pbw_sep = "Separate (MP + PBW)",
      pfvc_sep = "Separate (MP + PFVC)"
    )
  }
  
  models <- purrr::map(
    specification_names,
    ~ get_prog_model(family, .x, adjustment, normalization_spec)
  )
  reference_rows <- rownames(model.frame(models[[1]]))
  
  tibble(
    family = family, normalization_spec = normalization_spec, adjustment = adjustment,
    n_distinct_sample_sizes = n_distinct(purrr::map_int(models, nobs)),
    all_models_use_same_rows = all(purrr::map_lgl(
      models,
      ~ identical(rownames(model.frame(.x)), reference_rows)
    )),
    locked_vs_separate_pbw_df_difference = attr(logLik(models$pbw_sep), "df") -
      attr(logLik(models$pbw_lock), "df"),
    locked_vs_separate_pfvc_df_difference = attr(logLik(models$pfvc_sep), "df") -
      attr(logLik(models$pfvc_lock), "df")
  )
})

decomposition_checks <- list(
  expected_decomposition_rows_created = nrow(decomp_tbl) == nrow(decomposition_grid),
  decomposition_result_keys_unique = !anyDuplicated(
    decomp_tbl[c("site", "family", "normalization_spec", "adjustment")]
  ),
  decomposition_validation_rows_match_grid = nrow(decomposition_model_row_qc) == nrow(decomposition_grid),
  decomposition_validation_keys_unique = !anyDuplicated(
    decomposition_model_row_qc[c("family", "normalization_spec", "adjustment")]
  ),
  decomposition_covers_all_requested_normalizations = setequal(
    unique(decomp_tbl$normalization_spec),
    normalization_adjustment_specs$normalization_spec
  ),
  all_decomposition_model_sets_use_one_sample_size = all(
    decomposition_model_row_qc$n_distinct_sample_sizes == 1L
  ),
  all_decomposition_model_sets_use_same_rows = all(
    decomposition_model_row_qc$all_models_use_same_rows
  ),
  locked_vs_separate_pbw_models_differ_by_one_parameter = all(
    decomposition_model_row_qc$locked_vs_separate_pbw_df_difference == 1L
  ),
  locked_vs_separate_pfvc_models_differ_by_one_parameter = all(
    decomposition_model_row_qc$locked_vs_separate_pfvc_df_difference == 1L
  ),
  usable_lrt_p_values_valid = check_subset(
    decomp_tbl,
    decomp_tbl$computationally_usable,
    function(data) {
      all(
        data$form_pbw_lrt_p >= 0 & data$form_pbw_lrt_p <= 1 &
          data$form_pfvc_lrt_p >= 0 & data$form_pfvc_lrt_p <= 1
      )
    }
  ),
  usable_aic_comparisons_finite = check_subset(
    decomp_tbl,
    decomp_tbl$computationally_usable,
    function(data) {
      all(
        is.finite(data$phys_locked_dAIC) &
          is.finite(data$phys_separate_dAIC) &
          is.finite(data$aic_pbw_sep) &
          is.finite(data$aic_pfvc_sep)
      )
    }
  ),
  unusable_comparisons_suppressed = check_subset(
    decomp_tbl,
    !decomp_tbl$computationally_usable,
    function(data) {
      all(
        is.na(data$form_pbw_lrt_p) &
          is.na(data$form_pfvc_lrt_p) &
          is.na(data$phys_locked_dAIC) &
          is.na(data$phys_separate_dAIC) &
          is.na(data$aic_pbw_sep) &
          is.na(data$aic_pfvc_sep)
      )
    }
  )
)

if (!isTRUE(all(unlist(decomposition_checks, use.names = FALSE)))) {
  print(decomp_tbl, n = Inf, width = Inf)
  print(decomposition_model_row_qc, n = Inf, width = Inf)
  print(decomposition_checks)
}

stopifnot(isTRUE(all(unlist(decomposition_checks, use.names = FALSE))))

write_csv(decomp_tbl, form_physiology_output_path)

form_physiology_output_qc <- tibble(
  output = "Form-versus-physiology decomposition",
  path = form_physiology_output_path,
  exists = file.exists(form_physiology_output_path),
  size_bytes = if (file.exists(form_physiology_output_path)) {
    file.info(form_physiology_output_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

form_physiology_output_checks <- list(
  form_physiology_output_exists = form_physiology_output_qc$exists,
  form_physiology_output_nonempty = form_physiology_output_qc$nonempty
)

stopifnot(isTRUE(all(unlist(form_physiology_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$prognostic_utility$form_vs_physiology <- list(
  comparison_grid = decomposition_grid,
  normalization_specifications = normalization_adjustment_specs,
  comparison_set_definition = "Mechanic family × normalization/size specification × adjustment level",
  interpretation = paste0(
    "Locked-versus-separate comparisons use one-degree-of-freedom nested likelihood-ratio tests. ",
    "PBW-versus-PFVC comparisons within a form are nonnested and use AIC differences. ",
    "Positive physiology delta AIC values favor PFVC. Computationally unstable comparisons are suppressed."
  ),
  results = decomp_tbl,
  sample_and_nesting_validation = decomposition_model_row_qc,
  output = form_physiology_output_qc
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  decomposition_checks,
  form_physiology_output_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# Demographic effect modification of locked normalization metrics
# =============================================================================
# Each reduced model includes the locked metric and primary adjustment covariates. The expanded model adds the metric's
# interactions with continuous age, sex, and race. A nonsignificant joint test means that this analysis did not detect
# effect modification; it does not establish demographic invariance or a universal normalization target.

demographic_modification_output_path <- file.path(
  table_dir, paste0("norm_demographic_invariance_", site_name, ".csv")
)

demographic_modification_analysis_label <- "Demographic effect-modification test"

invar_specs <- tribble(
  ~family, ~exposure, ~normalization, ~locked_spec,
  "Elastance", "log(ers_pbw)", "Ers x PBW (Goligher)", "PBW-locked (Goligher)",
  "Elastance", "log(ers_pfvc)", "Ers x PFVC", "PFVC-locked",
  "Mechanical power", "log(mp_pbw)", "MP/PBW (Gattinoni)", "PBW-locked (Gattinoni)",
  "Mechanical power", "log(mp_pfvc)", "MP/PFVC", "PFVC-locked"
)

demographic_modification_fit_grid <- crossing(
  invar_specs,
  normalization_spec = normalization_adjustment_specs$normalization_spec
)

demographic_modification_grid_checks <- list(
  demographic_modification_grid_has_expected_rows = nrow(demographic_modification_fit_grid) ==
    nrow(invar_specs) * nrow(normalization_adjustment_specs),
  demographic_modification_grid_keys_unique = !anyDuplicated(
    demographic_modification_fit_grid[c("family", "normalization", "normalization_spec")]
  )
)

stopifnot(isTRUE(all(unlist(demographic_modification_grid_checks, use.names = FALSE))))

summarise_demographic_modification_model <- function(fitted_model, model_role) {
  outcome <- model.response(model.frame(fitted_model))
  model_matrix <- model.matrix(fitted_model)
  coefficient_table <- summary(fitted_model)$coefficients
  aliased_coefficients <- is.na(coef(fitted_model))
  n_events <- sum(outcome == 1L)
  n_parameters_excluding_intercept <- fitted_model$rank - 1L
  
  tibble(
    model_role = model_role, formula = paste(deparse(formula(fitted_model)), collapse = " "),
    n_obs = nobs(fitted_model), n_events = n_events,
    n_nonevents = sum(outcome == 0L),
    n_parameters_excluding_intercept = n_parameters_excluding_intercept,
    events_per_parameter = n_events / n_parameters_excluding_intercept,
    converged = isTRUE(fitted_model$converged), boundary = isTRUE(fitted_model$boundary),
    iterations = fitted_model$iter, model_matrix_columns = ncol(model_matrix),
    model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    n_aliased_coefficients = sum(aliased_coefficients),
    coefficients_finite = all(is.finite(coef(fitted_model)[!aliased_coefficients])),
    standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
    minimum_fitted_probability = min(fitted(fitted_model)),
    maximum_fitted_probability = max(fitted(fitted_model)),
    n_extreme_fitted_probabilities = sum(
      fitted(fitted_model) < EXTREME_FITTED_PROBABILITY_THRESHOLD |
        fitted(fitted_model) > 1 - EXTREME_FITTED_PROBABILITY_THRESHOLD
    )
  ) %>%
    mutate(
      computationally_usable = converged & !boundary & full_rank & n_aliased_coefficients == 0L &
        coefficients_finite & standard_errors_finite & n_extreme_fitted_probabilities == 0L
    )
}

invar_fits <- pmap(
  demographic_modification_fit_grid,
  function(family, exposure, normalization, locked_spec, normalization_spec) {
    analysis_data <- prognostic_analysis_data[[family]]
    
    registry_row <- normalization_adjustment_specs %>%
      filter(.data$normalization_spec == .env$normalization_spec)
    
    stopifnot(nrow(registry_row) == 1L)
    
    reduced_formula <- as.formula(
      paste(
        "mortality_event_60 ~", exposure, "+",
        get_normalization_adjustment(normalization_spec, "primary")
      )
    )
    
    expanded_formula <- as.formula(
      paste0(
        "mortality_event_60 ~ ", exposure, " * (age + sex_category + race_category) + ",
        get_normalization_adjustment(normalization_spec, "minimal")
      )
    )
  
  reduced_model <- glm(
    reduced_formula, data = analysis_data, family = binomial, x = TRUE, y = TRUE, model = TRUE
  )
  
  expanded_model <- glm(
    expanded_formula, data = analysis_data, family = binomial, x = TRUE, y = TRUE, model = TRUE
  )
  
  reduced_qc <- summarise_demographic_modification_model(reduced_model, "Reduced primary-adjusted model")
  expanded_qc <- summarise_demographic_modification_model(expanded_model, "Expanded interaction model")
  same_rows <- identical(rownames(model.frame(reduced_model)), rownames(model.frame(expanded_model)))
  df_difference <- attr(logLik(expanded_model), "df") - attr(logLik(reduced_model), "df")
  comparison_usable <- reduced_qc$computationally_usable && expanded_qc$computationally_usable && same_rows
  
  interaction_comparison <- if (comparison_usable) {
    anova(reduced_model, expanded_model, test = "LRT")
  } else {
    NULL
  }
  
  comparison_status <- if (!comparison_usable) {
    "Suppressed: computational instability or separation"
  } else if (expanded_qc$events_per_parameter < HIGHLY_SPARSE_PROGNOSTIC_EVENTS_PER_PARAMETER) {
    "Exploratory: highly sparse events"
  } else if (expanded_qc$events_per_parameter < MIN_PROGNOSTIC_EVENTS_PER_PARAMETER) {
    "Exploratory: sparse events"
  } else {
    "Eligible for conventional comparison"
  }
  
  demo_interaction_df <- if (comparison_usable) {
    as.numeric(interaction_comparison$Df[2])
  } else {
    NA_real_
  }
  
  demo_interaction_p <- if (comparison_usable) {
    as.numeric(interaction_comparison$`Pr(>Chi)`[2])
  } else {
    NA_real_
  }
  
  interaction_interpretation <- if (!comparison_usable) {
    "Joint interaction test suppressed because at least one model was computationally unstable."
  } else if (demo_interaction_p < 0.05) {
    "Evidence of demographic effect modification in the joint interaction test."
  } else {
    "No demographic effect modification detected; this does not establish invariance or universality."
  }
  
  stopifnot(
    length(demo_interaction_df) == 1L,
    length(demo_interaction_p) == 1L,
    length(interaction_interpretation) == 1L,
    comparison_usable || (is.na(demo_interaction_df) && is.na(demo_interaction_p))
  )
  
  result <- tibble(
    site = site_name,
    analysis_label = demographic_modification_analysis_label,
    family = family,
    normalization = normalization,
    normalization_spec = normalization_spec,
    normalization_variable = registry_row$normalization_variable[[1]],
    normalization_label = registry_row$normalization_label[[1]],
    normalization_kind = registry_row$normalization_kind[[1]],
    adjustment = "primary", adjusted = adjustment_labels[["primary"]],
    n = nobs(reduced_model),
    n_events = sum(model.response(model.frame(reduced_model)) == 1L),
    computationally_usable = comparison_usable,
    comparison_status = comparison_status,
    demo_interaction_df = demo_interaction_df,
    demo_interaction_p = demo_interaction_p,
    interpretation = interaction_interpretation
  )
  
  list(
    metadata = tibble(
      family = family, normalization = normalization, exposure = exposure, locked_spec = locked_spec,
      normalization_spec = normalization_spec,
      normalization_variable = registry_row$normalization_variable[[1]],
      normalization_label = registry_row$normalization_label[[1]],
      normalization_kind = registry_row$normalization_kind[[1]],
      same_rows = same_rows, df_difference = df_difference
    ),
    reduced_model = reduced_model,
    expanded_model = expanded_model,
    model_diagnostics = bind_rows(reduced_qc, expanded_qc) %>%
      mutate(
        family = family,
        normalization = normalization,
        normalization_spec = normalization_spec,
        normalization_variable = registry_row$normalization_variable[[1]],
        normalization_label = registry_row$normalization_label[[1]],
        normalization_kind = registry_row$normalization_kind[[1]],
        .before = 1
      ),
    comparison = if (is.null(interaction_comparison)) tibble() else
      as_tibble(interaction_comparison, rownames = "model_number"),
    result = result
  )
})

names(invar_fits) <- paste(
  demographic_modification_fit_grid$normalization_spec,
  demographic_modification_fit_grid$normalization,
  sep = " | "
)

invar_result_structure_qc <- purrr::imap_dfr(invar_fits, function(fit, fit_name) {
  tibble(
    fit_id = fit_name,
    normalization = fit$metadata$normalization[[1]],
    normalization_spec = fit$metadata$normalization_spec[[1]],
    result_is_tibble = tibble::is_tibble(fit$result),
    result_nrow = nrow(fit$result),
    result_ncol = ncol(fit$result),
    comparison_nrow = nrow(fit$comparison)
  )
})

stopifnot(
  nrow(invar_result_structure_qc) == nrow(demographic_modification_fit_grid),
  !anyDuplicated(names(invar_fits)),
  all(invar_result_structure_qc$result_is_tibble),
  all(invar_result_structure_qc$result_nrow == 1L)
)

invar_tbl <- purrr::map_dfr(invar_fits, "result")
demographic_modification_model_qc <- purrr::map_dfr(invar_fits, "model_diagnostics")
demographic_modification_comparison_qc <- purrr::map_dfr(invar_fits, "metadata")

demographic_modification_checks <- list(
  demographic_modification_expected_comparisons_evaluated =
    length(invar_fits) == nrow(demographic_modification_fit_grid) &&
    nrow(invar_tbl) == nrow(demographic_modification_fit_grid),
  demographic_modification_expected_model_count =
    nrow(demographic_modification_model_qc) == 2L * nrow(demographic_modification_fit_grid),
  demographic_modification_result_keys_unique = !anyDuplicated(
    invar_tbl[c("site", "family", "normalization", "normalization_spec", "adjustment")]
  ),
  demographic_modification_model_keys_unique = !anyDuplicated(
    demographic_modification_model_qc[c("family", "normalization", "normalization_spec", "model_role")]
  ),
  demographic_modification_all_requested_alternatives_present = setequal(
    unique(invar_tbl$normalization_spec),
    normalization_adjustment_specs$normalization_spec
  ),
  demographic_modification_four_locked_metrics_per_alternative = invar_tbl %>%
    count(normalization_spec) %>%
    summarise(
      valid = n() == nrow(normalization_adjustment_specs) && all(n == nrow(invar_specs))
    ) %>%
    pull(valid),
  demographic_modification_models_use_same_rows =
    all(demographic_modification_comparison_qc$same_rows),
  demographic_modification_interaction_blocks_add_four_parameters =
    all(demographic_modification_comparison_qc$df_difference == 4L),
  demographic_modification_models_use_expected_family_n = purrr::map2_lgl(
    demographic_modification_model_qc$family,
    demographic_modification_model_qc$n_obs,
    ~ .y == nrow(prognostic_analysis_data[[.x]])
  ) %>%
    all(),
  demographic_modification_model_matrices_full_rank =
    all(demographic_modification_model_qc$full_rank),
  demographic_modification_models_have_no_aliased_coefficients =
    all(demographic_modification_model_qc$n_aliased_coefficients == 0L),
  demographic_modification_coefficients_finite =
    all(demographic_modification_model_qc$coefficients_finite),
  demographic_modification_standard_errors_finite =
    all(demographic_modification_model_qc$standard_errors_finite),
  demographic_modification_usable_interaction_p_values_valid = check_subset(
    invar_tbl,
    invar_tbl$computationally_usable,
    function(data) {
      all(data$demo_interaction_p >= 0 & data$demo_interaction_p <= 1)
    }
  ),
  demographic_modification_unusable_interaction_tests_suppressed = check_subset(
    invar_tbl,
    !invar_tbl$computationally_usable,
    function(data) {
      all(is.na(data$demo_interaction_df) & is.na(data$demo_interaction_p))
    }
  )
)

if (!isTRUE(all(unlist(demographic_modification_checks, use.names = FALSE)))) {
  print(invar_tbl, n = Inf, width = Inf)
  print(demographic_modification_model_qc, n = Inf, width = Inf)
  print(demographic_modification_comparison_qc, n = Inf, width = Inf)
  print(demographic_modification_checks)
}

stopifnot(isTRUE(all(unlist(demographic_modification_checks, use.names = FALSE))))

write_csv(invar_tbl, demographic_modification_output_path)

demographic_modification_output_qc <- tibble(
  output = "Demographic effect-modification tests",
  path = demographic_modification_output_path,
  exists = file.exists(demographic_modification_output_path),
  size_bytes = if (file.exists(demographic_modification_output_path)) {
    file.info(demographic_modification_output_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

demographic_modification_output_checks <- list(
  demographic_modification_output_exists = demographic_modification_output_qc$exists,
  demographic_modification_output_nonempty = demographic_modification_output_qc$nonempty
)

stopifnot(isTRUE(all(unlist(demographic_modification_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$prognostic_utility$demographic_effect_modification <- list(
  fit_grid = demographic_modification_fit_grid,
  normalization_specifications = normalization_adjustment_specs,
  grid_checks = demographic_modification_grid_checks,
  interpretation = paste0(
    "The joint LRT evaluates interactions of the locked normalization metric with continuous age, sex, and race. ",
    "A nonsignificant test does not prove demographic invariance or establish a universal normalization target."
  ),
  results = invar_tbl,
  result_structure_validation = invar_result_structure_qc,
  comparison_validation = demographic_modification_comparison_qc,
  model_diagnostics = demographic_modification_model_qc,
  nested_model_comparisons = purrr::map(invar_fits, "comparison"),
  output = demographic_modification_output_qc
)

normalization_analysis_diagnostics$warnings$demographic_effect_modification <- invar_tbl %>%
  filter(comparison_status != "Eligible for conventional comparison")

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  demographic_modification_grid_checks,
  demographic_modification_checks,
  demographic_modification_output_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# Part 2 figure: prognostic fit and discrimination
# =============================================================================

prognostic_figure_path <- file.path(figure_dir, paste0("norm_prognostic_fit_", site_name, ".pdf"))

prognostic_spec_levels <- unique(prog_specs$spec)

prognostic_plot_data <- prog_tbl %>%
  filter(
    family_computationally_usable,
    is.finite(delta_aic),
    is.finite(auc)
  ) %>%
  mutate(
    spec = factor(as.character(spec), levels = rev(prognostic_spec_levels)),
    normalization_label = factor(
      normalization_label,
      levels = normalization_adjustment_specs$normalization_label
    ),
    adjusted = factor(adjusted, levels = unname(adjustment_labels))
  )

# Keep every mechanic-family/alternative panel, including panels with no exported comparison values.
prognostic_facet_scaffold <- crossing(
  family = unique(prog_specs$family),
  normalization_spec = normalization_adjustment_specs$normalization_spec
) %>%
  left_join(
    normalization_adjustment_specs %>% select(normalization_spec, normalization_label),
    by = "normalization_spec"
  ) %>%
  mutate(
    normalization_label = factor(
      normalization_label,
      levels = normalization_adjustment_specs$normalization_label
    ),
    spec = factor("PFVC-locked", levels = rev(prognostic_spec_levels)),
    aic_x = 0,
    auc_x = if (nrow(prognostic_plot_data) > 0L) median(prognostic_plot_data$auc) else 0.5
  )

prognostic_suppressed_panels <- prognostic_family_status %>%
  filter(!all_models_computationally_usable) %>%
  group_by(family, normalization_spec, normalization_label) %>%
  summarise(
    suppressed_adjustments = paste(adjusted, collapse = " and "),
    status = paste(unique(family_comparison_status), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    normalization_label = factor(
      normalization_label,
      levels = normalization_adjustment_specs$normalization_label
    ),
    spec = factor("PFVC-locked", levels = rev(prognostic_spec_levels)),
    aic_x = 0,
    auc_x = if (nrow(prognostic_plot_data) > 0L) median(prognostic_plot_data$auc) else 0.5,
    label = str_wrap(
      paste0(
        suppressed_adjustments, " comparisons suppressed. ",
        "Models failed the existing computational-use criteria."
      ),
      width = 28
    )
  )

prognostic_exploratory_groups <- prognostic_family_status %>%
  filter(
    all_models_computationally_usable,
    family_comparison_status != "Eligible for conventional comparison"
  ) %>%
  transmute(
    family,
    normalization_spec,
    normalization_label,
    adjusted,
    status = family_comparison_status
  )

prognostic_figure_caption <- if (nrow(prognostic_exploratory_groups) > 0L) {
  str_wrap(
    paste0(
      "Displayed comparisons are exploratory where event counts are sparse. ",
      "Suppressed panels represent computationally unstable model families. ",
      "Mortality models are event-limited; computational stability does not establish statistical reliability. ",
      "Apparent C-statistics are in-sample descriptive comparisons."
    ),
    width = 115
  )
} else {
  str_wrap(
    paste0(
      "Mortality models are event-limited; computational stability does not establish statistical reliability. ",
      "Apparent C-statistics are in-sample descriptive comparisons. ",
      "All displayed model families passed the prespecified computational checks."
    ),
    width = 115
  )
}

p_aic <- ggplot() +
  geom_blank(
    data = prognostic_facet_scaffold,
    aes(x = aic_x, y = spec)
  ) +
  geom_point(
    data = prognostic_plot_data,
    aes(x = delta_aic, y = spec, color = adjusted),
    size = 2.7,
    position = position_dodge(width = 0.55)
  ) +
  geom_text(
    data = prognostic_suppressed_panels,
    aes(x = aic_x, y = spec, label = label),
    inherit.aes = FALSE,
    color = "#B2182B",
    size = 3.1,
    lineheight = 1.05,
    hjust = 0
  ) +
  facet_grid(normalization_label ~ family, scales = "free_y", drop = FALSE) +
  scale_color_manual(
    values = c("Minimally adjusted" = "#0072B2", "Primary adjusted" = "#D55E00"),
    name = NULL,
    drop = FALSE
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(
    title = "Exploratory relative prognostic fit of normalization specifications",
    subtitle = str_wrap(
      paste0(
        site_name,
        " — delta AIC is calculated relative to the lowest-AIC specification within each mechanic family ",
        "and adjustment set for 60-day mortality."
      ),
      width = 105
    ),
    x = "Delta AIC within family and adjustment set",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 8)),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
  )

p_auc <- ggplot() +
  geom_blank(
    data = prognostic_facet_scaffold,
    aes(x = auc_x, y = spec)
  ) +
  geom_point(
    data = prognostic_plot_data,
    aes(x = auc, y = spec, color = adjusted),
    size = 2.7,
    position = position_dodge(width = 0.55)
  ) +
  geom_text(
    data = prognostic_suppressed_panels,
    aes(x = auc_x, y = spec, label = label),
    inherit.aes = FALSE,
    color = "#B2182B",
    size = 3.1,
    lineheight = 1.05,
    hjust = 0.5
  ) +
  facet_grid(normalization_label ~ family, scales = "free_y", drop = FALSE) +
  scale_color_manual(
    values = c("Minimally adjusted" = "#0072B2", "Primary adjusted" = "#D55E00"),
    guide = "none",
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = if (nrow(prognostic_plot_data) > 0L) {
      range(prognostic_plot_data$auc) + c(-0.03, 0.03)
    } else {
      c(0, 1)
    },
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(
    title = "Apparent in-sample discrimination",
    subtitle = str_wrap(
      "Apparent C-statistics are descriptive and optimistic because they are evaluated in the model-development sample.",
      width = 100
    ),
    x = "Apparent C-statistic",
    y = NULL,
    caption = prognostic_figure_caption
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 8)),
    plot.caption = element_text(hjust = 0, margin = margin(t = 8)),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
  )

prognostic_figure <- p_aic / p_auc + patchwork::plot_layout(heights = c(1, 1))

ggsave(
  filename = prognostic_figure_path, plot = prognostic_figure,
  width = 14, height = 32, units = "in", limitsize = TRUE
)

prognostic_figure_qc <- tibble(
  output = "Prognostic fit and discrimination figure",
  path = prognostic_figure_path,
  exists = file.exists(prognostic_figure_path),
  size_bytes = if (file.exists(prognostic_figure_path)) {
    file.info(prognostic_figure_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

expected_prognostic_panel_keys <- prognostic_fit_grid %>%
  distinct(family, normalization_spec)

suppressed_prognostic_panel_keys <- prognostic_family_status %>%
  filter(!all_models_computationally_usable) %>%
  distinct(family, normalization_spec)

observed_suppression_annotation_keys <- prognostic_suppressed_panels %>%
  distinct(family, normalization_spec)

prognostic_figure_checks <- list(
  plotted_rows_are_computationally_usable = all(
    prognostic_plot_data$family_computationally_usable
  ),
  plotted_delta_aic_values_finite = all(is.finite(prognostic_plot_data$delta_aic)),
  plotted_auc_values_finite = all(is.finite(prognostic_plot_data$auc)),
  plotted_alternative_labels_complete = !any(is.na(prognostic_plot_data$normalization_label)),
  scaffold_alternative_labels_complete = !any(is.na(prognostic_facet_scaffold$normalization_label)),
  prognostic_scaffold_keys_unique = !anyDuplicated(
    prognostic_facet_scaffold[c("family", "normalization_spec")]
  ),
  all_expected_panels_retained_in_scaffold = nrow(anti_join(
    expected_prognostic_panel_keys,
    prognostic_facet_scaffold,
    by = c("family", "normalization_spec")
  )) == 0L && nrow(prognostic_facet_scaffold) == nrow(expected_prognostic_panel_keys),
  all_suppressed_panels_annotated = nrow(anti_join(
    suppressed_prognostic_panel_keys,
    observed_suppression_annotation_keys,
    by = c("family", "normalization_spec")
  )) == 0L,
  no_unexpected_suppression_annotations = nrow(anti_join(
    observed_suppression_annotation_keys,
    suppressed_prognostic_panel_keys,
    by = c("family", "normalization_spec")
  )) == 0L,
  prognostic_plot_objects_created = !is.null(p_aic) && !is.null(p_auc) && !is.null(prognostic_figure),
  prognostic_figure_exists = prognostic_figure_qc$exists,
  prognostic_figure_nonempty = prognostic_figure_qc$nonempty
)

if (!isTRUE(all(unlist(prognostic_figure_checks, use.names = FALSE)))) {
  print(prognostic_plot_data, n = Inf, width = Inf)
  print(prognostic_suppressed_panels, n = Inf, width = Inf)
  print(prognostic_figure_qc, n = Inf, width = Inf)
  print(prognostic_figure_checks)
}

stopifnot(isTRUE(all(unlist(prognostic_figure_checks, use.names = FALSE))))

normalization_analysis_diagnostics$prognostic_utility$figure <- list(
  plotted_data = prognostic_plot_data,
  panel_scaffold = prognostic_facet_scaffold,
  suppressed_panels = prognostic_suppressed_panels,
  exploratory_groups = prognostic_exploratory_groups,
  layout = list(
    row_variable = "normalization_label",
    column_variable = "family",
    plot_blocks = c("Delta AIC", "Apparent C-statistic"),
    width_inches = 14,
    height_inches = 32,
    visual_review_status = "Not established by runtime checks"
  ),
  output = prognostic_figure_qc
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  prognostic_figure_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)


# =============================================================================
# PART 2b -- Modification of the mechanic association by age and predicted FVC
# =============================================================================
# This analysis evaluates whether the association between each respiratory mechanic and 60-day mortality varies across
# continuous age and predicted FVC. Age may partly reflect developmental and recoil-related differences, whereas PFVC
# represents predicted lung volume. These modifiers are correlated, especially in pediatric patients, so simultaneous
# interaction estimates require cautious interpretation.
#
# Every model retains VT/PBW, pSOFA, SF ratio, and BMI. The second adjustment level additionally includes sex and race.
# Age is not included in the adjustment string because it is introduced explicitly through age_c in the ladder rungs.

age_ref <- 10

# Retain the McLaren reference strings for existing diagnostic fields.
# Specification-aware ladder fits use get_ladder_adjustment().
ladder_adjustment_specs <- c(
  minimal = base_cov,
  demographic = paste(base_cov, "+ sex_category + race_category")
)

ladder_adjustment_labels <- c(
  minimal = "Minimally adjusted",
  demographic = "Minimally adjusted + sex/race"
)

get_ladder_adjustment <- function(normalization_spec, adjustment_id) {
  if (length(adjustment_id) != 1L || !adjustment_id %in% names(ladder_adjustment_labels)) {
    stop("Unknown ladder adjustment specification: ", paste(adjustment_id, collapse = ", "))
  }
  
  covariates <- get_normalization_adjustment(normalization_spec, "minimal")
  
  if (adjustment_id == "demographic") {
    covariates <- paste(covariates, "+ sex_category + race_category")
  }
  
  covariates
}

ladder_covariate_grid <- crossing(
  normalization_adjustment_specs,
  adjustment = names(ladder_adjustment_labels)
) %>%
  mutate(
    adjusted = unname(ladder_adjustment_labels[adjustment]),
    covariates = purrr::map2_chr(normalization_spec, adjustment, get_ladder_adjustment)
  )

ladder_covariate_checks <- list(
  ladder_covariate_grid_has_expected_rows = nrow(ladder_covariate_grid) ==
    nrow(normalization_adjustment_specs) * length(ladder_adjustment_labels),
  ladder_covariate_grid_keys_unique = !anyDuplicated(
    ladder_covariate_grid[c("normalization_spec", "adjustment")]
  ),
  ladder_covariates_do_not_duplicate_age = all(purrr::map_lgl(
    ladder_covariate_grid$covariates,
    ~ !any(c("age", "age_c", "age_at_admission") %in% all.vars(as.formula(paste("~", .x))))
  )),
  ladder_covariates_contain_exactly_the_selected_alternative = all(purrr::map2_lgl(
    ladder_covariate_grid$covariates,
    ladder_covariate_grid$normalization_variable,
    function(covariates, expected_variable) {
      included_alternatives <- intersect(
        all.vars(as.formula(paste("~", covariates))),
        normalization_adjustment_vars
      )
      
      length(included_alternatives) == 1L && identical(included_alternatives, expected_variable)
    }
  ))
)

stopifnot(isTRUE(all(unlist(ladder_covariate_checks, use.names = FALSE))))


# Centering makes the mechanic main effect correspond to age 10 and mean log-PFVC. It changes coefficient
# interpretation but does not alter model fit, likelihood-ratio tests, or fitted values.
ladder_rhs <- c(
  a = "mech_c",
  b = "mech_c + age_c",
  c = "mech_c * age_c",
  sep = "mech_c + pfvc_c",
  sep_int = "mech_c + pfvc_c + mech_c:pfvc_c",
  d = "mech_c * age_c + pfvc_c",
  d_int = "mech_c * age_c + pfvc_c + mech_c:pfvc_c",
  e = "mech_c * age_c + pfvc_c * age_c",
  full = "mech_c * age_c + pfvc_c * age_c + mech_c:pfvc_c"
)

exps <- tribble(
  ~family, ~mx, ~data_name,
  "Estimated driving pressure", "dp", "dp_data",
  "Elastance", "ers", "ers_data",
  "Mechanical power", "mechanical_power", "mp_data"
)

data_map <- list(
  dp_data = dp_data,
  ers_data = ers_data,
  mp_data = mp_data
)

ladder_common_vars <- unique(c(
  ".source_row", "hospitalization_id", "mortality_event_60", "pfvc",
  normalization_adjustment_vars,
  "sofa_total", "sf10", "bmi", "age", "sex_category", "race_category"
))

prepare_ladder_data <- function(family, mx, data_name) {
  source_data <- data_map[[data_name]]
  required_vars <- unique(c(ladder_common_vars, mx))
  missing_vars <- setdiff(required_vars, names(source_data))
  
  if (length(missing_vars) > 0L) {
    stop("The ", family, " ladder analysis is missing variable(s): ", paste(missing_vars, collapse = ", "))
  }
  
  analysis_data <- source_data %>%
    select(all_of(required_vars)) %>%
    drop_na(all_of(required_vars)) %>%
    filter(
      is.finite(.data[[mx]]), .data[[mx]] > 0,
      is.finite(pfvc), pfvc > 0,
      is.finite(age)
    )
  
  if (nrow(analysis_data) == 0L) {
    stop("No complete observations are available for the ", family, " interaction ladder.")
  }
  
  mean_log_mechanic <- mean(log(analysis_data[[mx]]))
  mean_log_pfvc <- mean(log(analysis_data$pfvc))
  
  analysis_data <- analysis_data %>%
    mutate(
      mech_c = log(.data[[mx]]) - mean_log_mechanic,
      pfvc_c = log(pfvc) - mean_log_pfvc,
      age_c = age - age_ref
    )
  
  list(
    data = analysis_data,
    family = family,
    mechanic_variable = mx,
    source_data_name = data_name,
    mean_log_mechanic = mean_log_mechanic,
    mean_log_pfvc = mean_log_pfvc,
    age_reference = age_ref
  )
}

ladder_prepared_data <- pmap(
  exps,
  function(family, mx, data_name) {
    prepare_ladder_data(family = family, mx = mx, data_name = data_name)
  }
)

names(ladder_prepared_data) <- exps$family

ladder_sample_missingness_qc <- purrr::imap_dfr(ladder_prepared_data, function(prepared, family_name) {
  source_data <- data_map[[prepared$source_data_name]]
  required_vars <- unique(c(ladder_common_vars, prepared$mechanic_variable))
  
  purrr::map_dfr(required_vars, function(variable) {
    tibble(
      family = family_name,
      variable = variable,
      n_metric_frame = nrow(source_data),
      n_missing = sum(is.na(source_data[[variable]])),
      percent_missing = 100 * mean(is.na(source_data[[variable]]))
    )
  })
})

ladder_excluded_rows <- purrr::imap_dfr(ladder_prepared_data, function(prepared, family_name) {
  source_data <- data_map[[prepared$source_data_name]]
  required_vars <- unique(c(ladder_common_vars, prepared$mechanic_variable))
  mechanic_values <- source_data[[prepared$mechanic_variable]]
  
  eligible_rows <- complete.cases(source_data[required_vars]) &
    is.finite(mechanic_values) & mechanic_values > 0 &
    is.finite(source_data$pfvc) & source_data$pfvc > 0 &
    is.finite(source_data$age)
  
  excluded_indices <- which(!eligible_rows)
  
  tibble(
    family = rep(family_name, length(excluded_indices)),
    source_row = source_data$.source_row[excluded_indices],
    hospitalization_id = source_data$hospitalization_id[excluded_indices],
    exclusion_reason = purrr::map_chr(excluded_indices, function(row_index) {
      missing_vars <- required_vars[purrr::map_lgl(
        required_vars,
        ~ is.na(source_data[[.x]][row_index])
      )]
      
      reasons <- character()
      
      if (length(missing_vars) > 0L) {
        reasons <- c(reasons, paste0("Missing: ", paste(missing_vars, collapse = ", ")))
      }
      
      if (!is.finite(mechanic_values[[row_index]]) || mechanic_values[[row_index]] <= 0) {
        reasons <- c(reasons, "Mechanic not finite and positive")
      }
      
      if (!is.finite(source_data$pfvc[[row_index]]) || source_data$pfvc[[row_index]] <= 0) {
        reasons <- c(reasons, "PFVC not finite and positive")
      }
      
      if (!is.finite(source_data$age[[row_index]])) {
        reasons <- c(reasons, "Age not finite")
      }
      
      paste(reasons, collapse = "; ")
    })
  )
})

ladder_sample_qc <- purrr::imap_dfr(ladder_prepared_data, function(prepared, family_name) {
  analysis_data <- prepared$data
  source_data <- data_map[[prepared$source_data_name]]
  
  tibble(
    family = family_name,
    mechanic_variable = prepared$mechanic_variable,
    n_metric_frame = nrow(source_data),
    n_analyzed = nrow(analysis_data),
    n_excluded = n_metric_frame - n_analyzed,
    percent_analyzed = 100 * n_analyzed / n_metric_frame,
    n_events = sum(analysis_data$mortality_event_60 == 1L),
    n_nonevents = sum(analysis_data$mortality_event_60 == 0L),
    event_percent = 100 * mean(analysis_data$mortality_event_60 == 1L),
    mean_log_mechanic = prepared$mean_log_mechanic,
    mean_log_pfvc = prepared$mean_log_pfvc,
    age_reference = prepared$age_reference,
    minimum_age = min(analysis_data$age),
    maximum_age = max(analysis_data$age),
    minimum_pfvc = min(analysis_data$pfvc),
    maximum_pfvc = max(analysis_data$pfvc)
  )
})

ladder_centering_qc <- purrr::imap_dfr(ladder_prepared_data, function(prepared, family_name) {
  analysis_data <- prepared$data
  
  tibble(
    family = family_name,
    mean_mech_c = mean(analysis_data$mech_c),
    mean_pfvc_c = mean(analysis_data$pfvc_c),
    age_c_at_reference = age_ref - prepared$age_reference,
    maximum_mech_reconstruction_difference = max(abs(
      exp(analysis_data$mech_c + prepared$mean_log_mechanic) -
        analysis_data[[prepared$mechanic_variable]]
    )),
    maximum_pfvc_reconstruction_difference = max(abs(
      exp(analysis_data$pfvc_c + prepared$mean_log_pfvc) - analysis_data$pfvc
    ))
  )
})

ladder_sample_checks <- list(
  all_three_mechanic_samples_created = identical(names(ladder_prepared_data), exps$family),
  ladder_common_sample_covers_all_alternatives = all(
    normalization_adjustment_vars %in% ladder_common_vars
  ),
  ladder_exclusion_details_reconcile = all(purrr::map2_lgl(
    ladder_sample_qc$family,
    ladder_sample_qc$n_excluded,
    function(family_name, expected_n) {
      sum(ladder_excluded_rows$family == family_name) == expected_n
    }
  )),
  ladder_exclusion_keys_unique = !anyDuplicated(
    ladder_excluded_rows[c("family", "source_row")]
  ),
  ladder_exclusion_reasons_nonempty = all(nzchar(ladder_excluded_rows$exclusion_reason)),
  ladder_samples_nonempty = all(ladder_sample_qc$n_analyzed > 0L),
  ladder_samples_have_events = all(ladder_sample_qc$n_events > 0L),
  ladder_samples_have_nonevents = all(ladder_sample_qc$n_nonevents > 0L),
  ladder_sample_ids_unique = all(purrr::map_lgl(
    ladder_prepared_data,
    ~ !anyDuplicated(.x$data$hospitalization_id)
  )),
  ladder_samples_complete = all(purrr::map_lgl(
    ladder_prepared_data,
    ~ all(complete.cases(.x$data))
  )),
  ladder_outcomes_binary = all(purrr::map_lgl(
    ladder_prepared_data,
    ~ all(.x$data$mortality_event_60 %in% c(0L, 1L))
  )),
  mechanics_positive = all(purrr::map_lgl(
    ladder_prepared_data,
    ~ all(.x$data[[.x$mechanic_variable]] > 0)
  )),
  pfvc_positive = all(purrr::map_lgl(
    ladder_prepared_data,
    ~ all(.x$data$pfvc > 0)
  )),
  centered_mechanic_means_near_zero = all(abs(ladder_centering_qc$mean_mech_c) <= 1e-12),
  centered_pfvc_means_near_zero = all(abs(ladder_centering_qc$mean_pfvc_c) <= 1e-12),
  age_reference_consistent = all(abs(ladder_centering_qc$age_c_at_reference) <= 1e-12),
  mechanic_reconstruction_valid = all(
    ladder_centering_qc$maximum_mech_reconstruction_difference <= 1e-10
  ),
  pfvc_reconstruction_valid = all(
    ladder_centering_qc$maximum_pfvc_reconstruction_difference <= 1e-10
  )
)

if (!isTRUE(all(unlist(ladder_sample_checks, use.names = FALSE)))) {
  print(ladder_sample_qc, n = Inf, width = Inf)
  print(ladder_centering_qc, n = Inf, width = Inf)
  print(ladder_sample_checks)
}

stopifnot(isTRUE(all(unlist(ladder_sample_checks, use.names = FALSE))))


normalization_analysis_diagnostics$mechanic_modification <- list(
  outcome = "60-day mortality",
  age_reference_years = age_ref,
  adjustment_specs = ladder_adjustment_specs,
  adjustment_labels = ladder_adjustment_labels,
  normalization_specifications = normalization_adjustment_specs,
  normalization_covariate_grid = ladder_covariate_grid,
  adjustment_specs_note = paste0(
    "adjustment_specs retains the original McLaren reference strings. ",
    "normalization_covariate_grid records the adjustment strings for the expanded ladder fits."
  ),
  common_sample_policy = paste0(
    "One common complete-case sample per mechanic family across all seven alternatives and both ladder adjustment ",
    "levels. All ladder rungs share this sample and its centering constants. ",
    "These exclusions do not change the upstream cohort or descriptive discordance/reclassification analyses."
  ),
  required_common_variables = ladder_common_vars,
  variable_missingness = ladder_sample_missingness_qc,
  excluded_rows = ladder_excluded_rows,
  ladder_formulas = ladder_rhs,
  mechanic_specifications = exps,
  sample_summary = ladder_sample_qc,
  centering_validation = ladder_centering_qc,
  source_rows = purrr::map(
    ladder_prepared_data,
    ~ .x$data %>% select(.source_row, hospitalization_id)
  )
)


normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  ladder_covariate_checks,
  ladder_sample_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

message("Script 05 | saving validated interaction-ladder inputs")
save_normalization_diagnostics(normalization_analysis_diagnostics)


# =============================================================================
# Fit and diagnose the interaction ladders
# =============================================================================

ladder_output_path <- file.path(table_dir, paste0("norm_age_interaction_ladder_", site_name, ".csv"))

summarise_ladder_model <- function(
    fitted_model, family_name, mechanic_variable, adjustment_id, rung, normalization_spec
) {
  registry_row <- normalization_adjustment_specs %>%
    filter(.data$normalization_spec == .env$normalization_spec)
  
  stopifnot(nrow(registry_row) == 1L)
  
  outcome <- model.response(model.frame(fitted_model))
  model_matrix <- model.matrix(fitted_model)
  coefficient_table <- summary(fitted_model)$coefficients
  aliased_coefficients <- is.na(coef(fitted_model))
  n_events <- sum(outcome == 1L)
  n_nonevents <- sum(outcome == 0L)
  n_parameters_excluding_intercept <- fitted_model$rank - 1L
  
  tibble(
    family = family_name, mechanic_variable = mechanic_variable,
    normalization_spec = normalization_spec,
    normalization_variable = registry_row$normalization_variable[[1]],
    normalization_label = registry_row$normalization_label[[1]],
    normalization_kind = registry_row$normalization_kind[[1]],
    adjustment = adjustment_id, adjusted = unname(ladder_adjustment_labels[[adjustment_id]]),
    rung = rung, rhs = unname(ladder_rhs[[rung]]),
    formula = paste(deparse(formula(fitted_model)), collapse = " "),
    n_obs = nobs(fitted_model), n_events = n_events, n_nonevents = n_nonevents,
    n_parameters_excluding_intercept = n_parameters_excluding_intercept,
    events_per_parameter = n_events / n_parameters_excluding_intercept,
    converged = isTRUE(fitted_model$converged), boundary = isTRUE(fitted_model$boundary),
    iterations = fitted_model$iter, model_matrix_columns = ncol(model_matrix),
    model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    n_aliased_coefficients = sum(aliased_coefficients),
    coefficients_finite = all(is.finite(coef(fitted_model)[!aliased_coefficients])),
    standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
    minimum_fitted_probability = min(fitted(fitted_model)),
    maximum_fitted_probability = max(fitted(fitted_model)),
    n_extreme_fitted_probabilities = sum(
      fitted(fitted_model) < EXTREME_FITTED_PROBABILITY_THRESHOLD |
        fitted(fitted_model) > 1 - EXTREME_FITTED_PROBABILITY_THRESHOLD
    ),
    log_likelihood = as.numeric(logLik(fitted_model)),
    aic = AIC(fitted_model)
  ) %>%
    mutate(
      computationally_usable = converged & !boundary & full_rank & n_aliased_coefficients == 0L &
        coefficients_finite & standard_errors_finite & n_extreme_fitted_probabilities == 0L,
      sparse_events = events_per_parameter < MIN_PROGNOSTIC_EVENTS_PER_PARAMETER,
      highly_sparse_events = events_per_parameter < HIGHLY_SPARSE_PROGNOSTIC_EVENTS_PER_PARAMETER,
      model_status = case_when(
        !converged ~ "Suppressed: model did not converge",
        boundary ~ "Suppressed: boundary fit",
        !full_rank | n_aliased_coefficients > 0L ~ "Suppressed: rank-deficient model",
        !coefficients_finite | !standard_errors_finite ~ "Suppressed: nonfinite inference",
        n_extreme_fitted_probabilities > 0L ~ "Suppressed: possible separation",
        highly_sparse_events ~ "Exploratory: highly sparse events",
        sparse_events ~ "Exploratory: sparse events",
        TRUE ~ "Eligible for conventional interpretation"
      )
    )
}

ladder_comparison_specs <- tribble(
  ~test_id, ~reduced_rung, ~expanded_rung, ~interpretation,
  "age_interaction", "b", "c", "Mechanic-by-age interaction",
  "pfvc_adds", "c", "d", "PFVC main effect beyond age-modified mechanic",
  "pfvc_x_age", "d", "e", "PFVC-by-age interaction",
  "size_int_alone", "sep", "sep_int", "Mechanic-by-PFVC interaction without age modification",
  "size_int_joint", "d", "d_int", "Mechanic-by-PFVC interaction beyond mechanic-by-age",
  "size_int_full", "e", "full", "Mechanic-by-PFVC interaction beyond both age interactions",
  "pfvc_age_in_full", "d_int", "full", "PFVC-by-age interaction beyond mechanic-by-PFVC"
)

fit_ladder <- function(family_name, adjustment_id, normalization_spec) {
  prepared <- ladder_prepared_data[[family_name]]
  analysis_data <- prepared$data
  mechanic_variable <- prepared$mechanic_variable
  
  registry_row <- normalization_adjustment_specs %>%
    filter(.data$normalization_spec == .env$normalization_spec)
  
  stopifnot(nrow(registry_row) == 1L)
  
  covariates <- get_ladder_adjustment(normalization_spec, adjustment_id)
  
  models <- purrr::imap(ladder_rhs, function(rhs, rung) {
    model_formula <- as.formula(paste("mortality_event_60 ~", rhs, "+", covariates))
    
    glm(
      model_formula, data = analysis_data, family = binomial,
      x = TRUE, y = TRUE, model = TRUE
    )
  })
  
  model_diagnostics <- purrr::imap_dfr(models, function(fitted_model, rung) {
    summarise_ladder_model(
      fitted_model = fitted_model,
      family_name = family_name,
      mechanic_variable = mechanic_variable,
      adjustment_id = adjustment_id,
      rung = rung,
      normalization_spec = normalization_spec
    )
  })
  
  reference_rows <- rownames(model.frame(models[[1]]))
  
  row_validation <- tibble(
    family = family_name, normalization_spec = normalization_spec, adjustment = adjustment_id,
    expected_n = nrow(analysis_data),
    n_distinct_model_sample_sizes = n_distinct(purrr::map_int(models, nobs)),
    all_models_use_expected_n = all(purrr::map_int(models, nobs) == nrow(analysis_data)),
    all_models_use_same_rows = all(purrr::map_lgl(
      models,
      ~ identical(rownames(model.frame(.x)), reference_rows)
    ))
  )
  
  comparisons <- pmap_dfr(
    ladder_comparison_specs,
    function(test_id, reduced_rung, expanded_rung, interpretation) {
      reduced_model <- models[[reduced_rung]]
      expanded_model <- models[[expanded_rung]]
      
      reduced_status <- model_diagnostics %>%
        filter(rung == reduced_rung)
      
      expanded_status <- model_diagnostics %>%
        filter(rung == expanded_rung)
      
      stopifnot(nrow(reduced_status) == 1L, nrow(expanded_status) == 1L)
      
      same_rows <- identical(
        rownames(model.frame(reduced_model)),
        rownames(model.frame(expanded_model))
      )
      
      df_difference <- attr(logLik(expanded_model), "df") - attr(logLik(reduced_model), "df")
      
      computationally_usable <- reduced_status$computationally_usable[[1]] &&
        expanded_status$computationally_usable[[1]] &&
        same_rows &&
        df_difference > 0L
      
      comparison <- if (computationally_usable) {
        anova(reduced_model, expanded_model, test = "LRT")
      } else {
        NULL
      }
      
      lrt_statistic <- if (computationally_usable) {
        as.numeric(comparison$Deviance[2])
      } else {
        NA_real_
      }
      
      lrt_df <- if (computationally_usable) {
        as.numeric(comparison$Df[2])
      } else {
        NA_real_
      }
      
      lrt_p_value <- if (computationally_usable) {
        as.numeric(comparison$`Pr(>Chi)`[2])
      } else {
        NA_real_
      }
      
      comparison_status <- if (!computationally_usable) {
        "Suppressed: at least one model is computationally unstable"
      } else if (
        reduced_status$highly_sparse_events[[1]] ||
        expanded_status$highly_sparse_events[[1]]
      ) {
        "Exploratory: highly sparse events"
      } else if (
        reduced_status$sparse_events[[1]] ||
        expanded_status$sparse_events[[1]]
      ) {
        "Exploratory: sparse events"
      } else {
        "Eligible for conventional interpretation"
      }
      
      tibble(
        family = family_name, mechanic_variable = mechanic_variable,
        normalization_spec = normalization_spec,
        normalization_variable = registry_row$normalization_variable[[1]],
        normalization_label = registry_row$normalization_label[[1]],
        normalization_kind = registry_row$normalization_kind[[1]],
        adjustment = adjustment_id, adjusted = unname(ladder_adjustment_labels[[adjustment_id]]),
        test_id = test_id, interpretation = interpretation,
        reduced_rung = reduced_rung, expanded_rung = expanded_rung,
        n_obs = nobs(reduced_model),
        n_events = sum(model.response(model.frame(reduced_model)) == 1L),
        same_rows = same_rows, df_difference = df_difference,
        computationally_usable = computationally_usable,
        comparison_status = comparison_status,
        lrt_statistic = lrt_statistic, lrt_df = lrt_df, lrt_p_value = lrt_p_value
      )
    }
  )
  
  list(
    family = family_name,
    mechanic_variable = mechanic_variable,
    normalization_spec = normalization_spec,
    normalization_variable = registry_row$normalization_variable[[1]],
    normalization_label = registry_row$normalization_label[[1]],
    normalization_kind = registry_row$normalization_kind[[1]],
    adjustment = adjustment_id,
    adjusted = unname(ladder_adjustment_labels[[adjustment_id]]),
    data = analysis_data,
    models = models,
    model_diagnostics = model_diagnostics,
    row_validation = row_validation,
    comparisons = comparisons,
    mean_log_mechanic = prepared$mean_log_mechanic,
    mean_log_pfvc = prepared$mean_log_pfvc,
    age_reference = prepared$age_reference
  )
}

ladder_fit_grid <- crossing(
  family = exps$family,
  normalization_spec = normalization_adjustment_specs$normalization_spec,
  adjustment = names(ladder_adjustment_labels)
)

message("Script 05 | fitting age/PFVC interaction ladders")
ladder_fits <- pmap(
  ladder_fit_grid,
  function(family, normalization_spec, adjustment) {
    fit_ladder(
      family_name = family,
      adjustment_id = adjustment,
      normalization_spec = normalization_spec
    )
  }
)

names(ladder_fits) <- paste(
  ladder_fit_grid$family,
  ladder_fit_grid$normalization_spec,
  ladder_fit_grid$adjustment,
  sep = " | "
)

ladder_fit_grid_checks <- list(
  ladder_grid_has_expected_number_of_sets = nrow(ladder_fit_grid) ==
    nrow(exps) * nrow(normalization_adjustment_specs) * length(ladder_adjustment_labels),
  ladder_grid_keys_unique = !anyDuplicated(
    ladder_fit_grid[c("family", "normalization_spec", "adjustment")]
  ),
  ladder_fit_names_unique = !anyDuplicated(names(ladder_fits)),
  ladder_fit_count_matches_grid = length(ladder_fits) == nrow(ladder_fit_grid),
  ladder_fit_metadata_matches_grid = all(purrr::map2_lgl(
    ladder_fits,
    seq_len(nrow(ladder_fit_grid)),
    function(fit, fit_index) {
      identical(fit$family, ladder_fit_grid$family[[fit_index]]) &&
        identical(fit$normalization_spec, ladder_fit_grid$normalization_spec[[fit_index]]) &&
        identical(fit$adjustment, ladder_fit_grid$adjustment[[fit_index]])
    }
  ))
)

stopifnot(isTRUE(all(unlist(ladder_fit_grid_checks, use.names = FALSE))))

normalization_analysis_diagnostics$mechanic_modification$fit_grid <- ladder_fit_grid

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  ladder_fit_grid_checks
)

ladder_model_qc <- purrr::map_dfr(ladder_fits, "model_diagnostics")
ladder_comparison_qc <- purrr::map_dfr(ladder_fits, "comparisons")
ladder_row_validation_qc <- purrr::map_dfr(ladder_fits, "row_validation")

ladder_cross_specification_rows_qc <- tibble(family = names(ladder_prepared_data)) %>%
  mutate(
    same_rows_across_alternatives_and_adjustments = purrr::map_lgl(family, function(family_name) {
      fit_indices <- which(ladder_fit_grid$family == family_name)
      expected_rows <- rownames(ladder_prepared_data[[family_name]]$data)
      
      all(purrr::map_lgl(ladder_fits[fit_indices], function(fit) {
        all(purrr::map_lgl(
          fit$models,
          ~ identical(rownames(model.frame(.x)), expected_rows)
        ))
      }))
    })
  )

ladder_common_row_checks <- list(
  ladder_cross_specification_validation_covers_all_families = setequal(
    ladder_cross_specification_rows_qc$family,
    exps$family
  ),
  ladder_alternatives_and_adjustments_share_common_family_rows = all(
    ladder_cross_specification_rows_qc$same_rows_across_alternatives_and_adjustments
  )
)

stopifnot(isTRUE(all(unlist(ladder_common_row_checks, use.names = FALSE))))

normalization_analysis_diagnostics$mechanic_modification$cross_specification_row_validation <-
  ladder_cross_specification_rows_qc

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  ladder_common_row_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

ladder_fit_integrity_checks <- list(
  ladder_expected_sets_fitted = length(ladder_fits) == nrow(ladder_fit_grid),
  ladder_model_keys_unique = !anyDuplicated(
    ladder_model_qc[c("family", "normalization_spec", "adjustment", "rung")]
  ),
  ladder_comparison_keys_unique = !anyDuplicated(
    ladder_comparison_qc[c("family", "normalization_spec", "adjustment", "test_id")]
  ),
  ladder_row_validation_count_matches_grid = nrow(ladder_row_validation_qc) == nrow(ladder_fit_grid),
  ladder_nine_models_per_set = ladder_model_qc %>%
    count(family, normalization_spec, adjustment) %>%
    summarise(
      passed = n() == nrow(ladder_fit_grid) && all(n == length(ladder_rhs))
    ) %>%
    pull(passed),
  ladder_expected_number_of_models_fitted = nrow(ladder_model_qc) ==
    nrow(ladder_fit_grid) * length(ladder_rhs),
  ladder_sets_use_one_sample_size = all(
    ladder_row_validation_qc$n_distinct_model_sample_sizes == 1L
  ),
  ladder_models_use_expected_sample_size = all(
    ladder_row_validation_qc$all_models_use_expected_n
  ),
  ladder_models_within_set_use_same_rows = all(
    ladder_row_validation_qc$all_models_use_same_rows
  ),
  ladder_model_matrices_full_rank = all(ladder_model_qc$full_rank),
  ladder_models_have_no_aliased_coefficients = all(
    ladder_model_qc$n_aliased_coefficients == 0L
  ),
  ladder_coefficients_finite = all(ladder_model_qc$coefficients_finite),
  ladder_standard_errors_finite = all(ladder_model_qc$standard_errors_finite),
  ladder_fitted_probability_ranges_valid = all(
    is.finite(ladder_model_qc$minimum_fitted_probability) &
      is.finite(ladder_model_qc$maximum_fitted_probability) &
      ladder_model_qc$minimum_fitted_probability >= 0 &
      ladder_model_qc$maximum_fitted_probability <= 1
  ),
  ladder_log_likelihoods_finite = all(is.finite(ladder_model_qc$log_likelihood)),
  ladder_aic_values_finite = all(is.finite(ladder_model_qc$aic)),
  ladder_seven_comparisons_per_set = ladder_comparison_qc %>%
    count(family, normalization_spec, adjustment) %>%
    summarise(
      passed = n() == nrow(ladder_fit_grid) && all(n == nrow(ladder_comparison_specs))
    ) %>%
    pull(passed),
  ladder_expected_number_of_comparisons = nrow(ladder_comparison_qc) ==
    nrow(ladder_fit_grid) * nrow(ladder_comparison_specs),
  ladder_nested_comparisons_use_same_rows = all(ladder_comparison_qc$same_rows),
  ladder_nested_comparisons_add_one_parameter = all(ladder_comparison_qc$df_difference == 1L),
  ladder_usable_lrt_results_valid = check_subset(
    ladder_comparison_qc,
    ladder_comparison_qc$computationally_usable,
    function(data) {
      all(
        is.finite(data$lrt_statistic) &
          data$lrt_statistic >= 0 &
          is.finite(data$lrt_df) &
          data$lrt_df == 1 &
          is.finite(data$lrt_p_value) &
          data$lrt_p_value >= 0 &
          data$lrt_p_value <= 1
      )
    }
  ),
  ladder_unusable_lrt_results_suppressed = check_subset(
    ladder_comparison_qc,
    !ladder_comparison_qc$computationally_usable,
    function(data) {
      all(is.na(data$lrt_statistic) & is.na(data$lrt_df) & is.na(data$lrt_p_value))
    }
  )
)

if (!isTRUE(all(unlist(ladder_fit_integrity_checks, use.names = FALSE)))) {
  print(
    ladder_model_qc %>%
      select(
        family, normalization_spec, normalization_label, adjustment, adjusted, rung, n_obs, n_events,
        n_parameters_excluding_intercept, events_per_parameter,
        converged, boundary, full_rank, n_extreme_fitted_probabilities,
        computationally_usable, model_status
      ),
    n = Inf,
    width = Inf
  )
  print(ladder_comparison_qc, n = Inf, width = Inf)
  print(ladder_row_validation_qc, n = Inf, width = Inf)
  print(ladder_fit_integrity_checks)
}

stopifnot(isTRUE(all(unlist(ladder_fit_integrity_checks, use.names = FALSE))))

# Create the federated ladder summary while preserving the original principal column names.
ladder_tbl <- ladder_fit_grid %>%
  pmap_dfr(function(family, normalization_spec, adjustment) {
    fit_name <- paste(family, normalization_spec, adjustment, sep = " | ")
    fit <- ladder_fits[[fit_name]]
    
    if (is.null(fit)) {
      stop("The expected ladder fit was not found: ", fit_name)
    }
    
    model_qc <- fit$model_diagnostics
    comparison_qc <- fit$comparisons
    
    get_model_value <- function(rung, variable) {
      value <- model_qc %>%
        filter(.data$rung == !!rung) %>%
        pull(all_of(variable))
      
      if (length(value) != 1L) {
        stop("Expected one `", variable, "` value for rung `", rung, "` in ", fit_name, ".")
      }
      
      value[[1]]
    }
    
    get_comparison_value <- function(test_id, variable) {
      value <- comparison_qc %>%
        filter(.data$test_id == !!test_id) %>%
        pull(all_of(variable))
      
      if (length(value) != 1L) {
        stop("Expected one `", variable, "` value for comparison `", test_id, "` in ", fit_name, ".")
      }
      
      value[[1]]
    }
    
    full_usable <- get_model_value("full", "computationally_usable")
    full_model <- fit$models$full
    
    tibble(
      site = site_name, family = family, mx = fit$mechanic_variable,
      normalization_spec = normalization_spec,
      normalization_variable = fit$normalization_variable,
      normalization_label = fit$normalization_label,
      normalization_kind = fit$normalization_kind,
      adjustment = adjustment, adjusted = fit$adjusted,
      n = nrow(fit$data),
      n_events = sum(fit$data$mortality_event_60 == 1L),
      aic_mech = if (get_model_value("a", "computationally_usable")) AIC(fit$models$a) else NA_real_,
      aic_mech_x_age = if (get_model_value("c", "computationally_usable")) AIC(fit$models$c) else NA_real_,
      aic_age_pfvc = if (get_model_value("d", "computationally_usable")) AIC(fit$models$d) else NA_real_,
      aic_joint = if (get_model_value("d_int", "computationally_usable")) AIC(fit$models$d_int) else NA_real_,
      aic_full = if (full_usable) AIC(full_model) else NA_real_,
      age_interaction_p = get_comparison_value("age_interaction", "lrt_p_value"),
      pfvc_adds_p = get_comparison_value("pfvc_adds", "lrt_p_value"),
      pfvc_x_age_p = get_comparison_value("pfvc_x_age", "lrt_p_value"),
      size_int_alone_p = get_comparison_value("size_int_alone", "lrt_p_value"),
      size_int_joint_p = get_comparison_value("size_int_joint", "lrt_p_value"),
      size_int_full_p = get_comparison_value("size_int_full", "lrt_p_value"),
      pfvc_age_in_full_p = get_comparison_value("pfvc_age_in_full", "lrt_p_value"),
      max_vif_full = if (full_usable) vif_max(full_model) else NA_real_,
      full_model_computationally_usable = full_usable,
      full_model_status = get_model_value("full", "model_status")
    )
  })

ladder_output_checks <- list(
  ladder_summary_has_expected_rows = nrow(ladder_tbl) == nrow(ladder_fit_grid),
  ladder_summary_keys_unique = !anyDuplicated(
    ladder_tbl[c("site", "family", "normalization_spec", "adjustment")]
  ),
  ladder_summary_covers_all_alternatives = setequal(
    unique(ladder_tbl$normalization_spec),
    normalization_adjustment_specs$normalization_spec
  ),
  ladder_summary_sample_sizes_match_prepared_families = all(purrr::map2_lgl(
    ladder_tbl$family,
    ladder_tbl$n,
    ~ .y == nrow(ladder_prepared_data[[.x]]$data)
  )),
  ladder_summary_event_counts_match_prepared_families = all(purrr::map2_lgl(
    ladder_tbl$family,
    ladder_tbl$n_events,
    ~ .y == sum(ladder_prepared_data[[.x]]$data$mortality_event_60 == 1L)
  )),
  ladder_summary_sample_sizes_positive = all(ladder_tbl$n > 0L),
  ladder_summary_event_counts_positive = all(ladder_tbl$n_events > 0L),
  ladder_unusable_full_models_have_suppressed_aic_and_vif = check_subset(
    ladder_tbl,
    !ladder_tbl$full_model_computationally_usable,
    function(data) {
      all(is.na(data$aic_full) & is.na(data$max_vif_full))
    }
  ),
  ladder_usable_full_models_have_finite_aic_and_vif = check_subset(
    ladder_tbl,
    ladder_tbl$full_model_computationally_usable,
    function(data) {
      all(is.finite(data$aic_full) & is.finite(data$max_vif_full))
    }
  )
)

if (!isTRUE(all(unlist(ladder_output_checks, use.names = FALSE)))) {
  print(ladder_tbl, n = Inf, width = Inf)
  print(ladder_output_checks)
}

stopifnot(isTRUE(all(unlist(ladder_output_checks, use.names = FALSE))))

write_csv(ladder_tbl, ladder_output_path)

ladder_file_qc <- tibble(
  output = "Mechanic age/PFVC interaction ladder",
  path = ladder_output_path,
  exists = file.exists(ladder_output_path),
  size_bytes = if (file.exists(ladder_output_path)) {
    file.info(ladder_output_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

ladder_file_checks <- list(
  ladder_output_exists = ladder_file_qc$exists,
  ladder_output_nonempty = ladder_file_qc$nonempty
)

stopifnot(isTRUE(all(unlist(ladder_file_checks, use.names = FALSE))))

normalization_analysis_diagnostics$mechanic_modification$model_diagnostics <- ladder_model_qc
normalization_analysis_diagnostics$mechanic_modification$comparison_diagnostics <- ladder_comparison_qc
normalization_analysis_diagnostics$mechanic_modification$row_validation <- ladder_row_validation_qc
normalization_analysis_diagnostics$mechanic_modification$summary <- ladder_tbl
normalization_analysis_diagnostics$mechanic_modification$output <- ladder_file_qc

normalization_analysis_diagnostics$warnings$mechanic_modification <- ladder_model_qc %>%
  filter(model_status != "Eligible for conventional interpretation")

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  ladder_fit_integrity_checks,
  ladder_output_checks,
  ladder_file_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

format_ladder_p <- function(p_value) {
  if (length(p_value) != 1L || !is.finite(p_value)) {
    return("suppressed")
  }
  
  format.pval(p_value, digits = 3, eps = 0.001)
}

# =============================================================================
# Balanced single-modifier models
# =============================================================================
# Each mechanic uses the previously validated common sample and centering constants. Both interaction models contain
# identical main effects:
#
#   mechanic + age + PFVC + VT/PBW + pSOFA + SF ratio + BMI + sex + race
#
# The age-modifier model adds mechanic × age. The PFVC-modifier model instead adds mechanic × PFVC. Each interaction is
# therefore evaluated separately while retaining the other modifier as an adjustment variable.

single_modifier_output_path <- file.path(
  table_dir, paste0("norm_single_modifier_diagnostics_", site_name, ".csv")
)

single_modifier_main_rhs <- paste(
  "mech_c + age_c + pfvc_c +",
  ladder_adjustment_specs[["demographic"]]
)

get_single_modifier_formulas <- function(normalization_spec) {
  main_rhs <- paste(
    "mech_c + age_c + pfvc_c +",
    get_ladder_adjustment(normalization_spec, "demographic")
  )
  
  list(
    main = as.formula(paste("mortality_event_60 ~", main_rhs)),
    age_modifier = as.formula(paste("mortality_event_60 ~", main_rhs, "+ mech_c:age_c")),
    pfvc_modifier = as.formula(paste("mortality_event_60 ~", main_rhs, "+ mech_c:pfvc_c"))
  )
}

# Retain the McLaren reference formulas; specification-aware fits obtain their own formulas.
single_modifier_formulas <- get_single_modifier_formulas("vt_mclaren_pbw")

summarise_single_modifier_model <- function(fitted_model, family_name, model_name, normalization_spec) {
  registry_row <- normalization_adjustment_specs %>%
    filter(.data$normalization_spec == .env$normalization_spec)
  
  stopifnot(nrow(registry_row) == 1L)
  
  outcome <- model.response(model.frame(fitted_model))
  model_matrix <- model.matrix(fitted_model)
  coefficient_table <- summary(fitted_model)$coefficients
  aliased_coefficients <- is.na(coef(fitted_model))
  n_events <- sum(outcome == 1L)
  n_nonevents <- sum(outcome == 0L)
  n_parameters_excluding_intercept <- fitted_model$rank - 1L
  
  tibble(
    family = family_name, model = model_name,
    normalization_spec = normalization_spec,
    normalization_variable = registry_row$normalization_variable[[1]],
    normalization_label = registry_row$normalization_label[[1]],
    normalization_kind = registry_row$normalization_kind[[1]],
    formula = paste(deparse(formula(fitted_model)), collapse = " "),
    n_obs = nobs(fitted_model), n_events = n_events, n_nonevents = n_nonevents,
    n_parameters_excluding_intercept = n_parameters_excluding_intercept,
    events_per_parameter = n_events / n_parameters_excluding_intercept,
    converged = isTRUE(fitted_model$converged), boundary = isTRUE(fitted_model$boundary),
    iterations = fitted_model$iter, model_matrix_columns = ncol(model_matrix),
    model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    n_aliased_coefficients = sum(aliased_coefficients),
    coefficients_finite = all(is.finite(coef(fitted_model)[!aliased_coefficients])),
    standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
    minimum_fitted_probability = min(fitted(fitted_model)),
    maximum_fitted_probability = max(fitted(fitted_model)),
    n_extreme_fitted_probabilities = sum(
      fitted(fitted_model) < EXTREME_FITTED_PROBABILITY_THRESHOLD |
        fitted(fitted_model) > 1 - EXTREME_FITTED_PROBABILITY_THRESHOLD
    ),
    log_likelihood = as.numeric(logLik(fitted_model)),
    aic = AIC(fitted_model)
  ) %>%
    mutate(
      computationally_usable = converged & !boundary & full_rank & n_aliased_coefficients == 0L &
        coefficients_finite & standard_errors_finite & n_extreme_fitted_probabilities == 0L,
      sparse_events = events_per_parameter < MIN_PROGNOSTIC_EVENTS_PER_PARAMETER,
      highly_sparse_events = events_per_parameter < HIGHLY_SPARSE_PROGNOSTIC_EVENTS_PER_PARAMETER,
      model_status = case_when(
        !converged ~ "Suppressed: model did not converge",
        boundary ~ "Suppressed: boundary fit",
        !full_rank | n_aliased_coefficients > 0L ~ "Suppressed: rank-deficient model",
        !coefficients_finite | !standard_errors_finite ~ "Suppressed: nonfinite inference",
        n_extreme_fitted_probabilities > 0L ~ "Suppressed: possible separation",
        highly_sparse_events ~ "Exploratory: highly sparse events",
        sparse_events ~ "Exploratory: sparse events",
        TRUE ~ "Eligible for conventional interpretation"
      )
    )
}

fit_single_modifier_models <- function(family_name, normalization_spec) {
  prepared <- ladder_prepared_data[[family_name]]
  analysis_data <- prepared$data
  
  registry_row <- normalization_adjustment_specs %>%
    filter(.data$normalization_spec == .env$normalization_spec)
  
  stopifnot(nrow(registry_row) == 1L)
  
  model_formulas <- get_single_modifier_formulas(normalization_spec)
  
  models <- purrr::map(model_formulas, function(model_formula) {
    glm(
      model_formula, data = analysis_data, family = binomial,
      x = TRUE, y = TRUE, model = TRUE
    )
  })
  
  model_diagnostics <- purrr::imap_dfr(models, function(fitted_model, model_name) {
    summarise_single_modifier_model(
      fitted_model = fitted_model,
      family_name = family_name,
      model_name = model_name,
      normalization_spec = normalization_spec
    )
  })
  
  reference_rows <- rownames(model.frame(models$main))
  
  row_validation <- tibble(
    family = family_name,
    normalization_spec = normalization_spec,
    expected_n = nrow(analysis_data),
    n_distinct_model_sample_sizes = n_distinct(purrr::map_int(models, nobs)),
    all_models_use_expected_n = all(purrr::map_int(models, nobs) == nrow(analysis_data)),
    all_models_use_same_rows = all(purrr::map_lgl(
      models,
      ~ identical(rownames(model.frame(.x)), reference_rows)
    ))
  )
  
  interaction_specs <- tribble(
    ~modifier, ~expanded_model,
    "Age", "age_modifier",
    "PFVC", "pfvc_modifier"
  )
  
  interaction_results <- pmap_dfr(interaction_specs, function(modifier, expanded_model) {
    reduced_model <- models$main
    interaction_model <- models[[expanded_model]]
    
    reduced_status <- model_diagnostics %>%
      filter(model == "main")
    
    expanded_status <- model_diagnostics %>%
      filter(model == expanded_model)
    
    stopifnot(nrow(reduced_status) == 1L, nrow(expanded_status) == 1L)
    
    same_rows <- identical(
      rownames(model.frame(reduced_model)),
      rownames(model.frame(interaction_model))
    )
    
    df_difference <- attr(logLik(interaction_model), "df") - attr(logLik(reduced_model), "df")
    
    computationally_usable <- reduced_status$computationally_usable[[1]] &&
      expanded_status$computationally_usable[[1]] &&
      same_rows &&
      df_difference == 1L
    
    comparison <- if (computationally_usable) {
      anova(reduced_model, interaction_model, test = "LRT")
    } else {
      NULL
    }
    
    interaction_term_candidates <- if (modifier == "Age") {
      c("mech_c:age_c", "age_c:mech_c")
    } else {
      c("mech_c:pfvc_c", "pfvc_c:mech_c")
    }
    
    interaction_term <- intersect(interaction_term_candidates, names(coef(interaction_model)))
    
    if (length(interaction_term) != 1L) {
      stop("Expected exactly one ", modifier, " interaction term for ", family_name, ".")
    }
    
    coefficient_table <- summary(interaction_model)$coefficients
    interaction_estimate <- coefficient_table[interaction_term, "Estimate"]
    interaction_standard_error <- coefficient_table[interaction_term, "Std. Error"]
    interaction_wald_p <- coefficient_table[interaction_term, "Pr(>|z|)"]
    
    comparison_status <- if (!computationally_usable) {
      "Suppressed: at least one model is computationally unstable"
    } else if (expanded_status$highly_sparse_events[[1]]) {
      "Exploratory: highly sparse events"
    } else if (expanded_status$sparse_events[[1]]) {
      "Exploratory: sparse events"
    } else {
      "Eligible for conventional interpretation"
    }
    
    tibble(
      family = family_name, modifier = modifier,
      normalization_spec = normalization_spec,
      normalization_variable = registry_row$normalization_variable[[1]],
      normalization_label = registry_row$normalization_label[[1]],
      normalization_kind = registry_row$normalization_kind[[1]],
      reduced_model = "main", expanded_model = expanded_model,
      interaction_term = interaction_term,
      n_obs = nobs(interaction_model),
      n_events = sum(model.response(model.frame(interaction_model)) == 1L),
      n_parameters_excluding_intercept = expanded_status$n_parameters_excluding_intercept[[1]],
      events_per_parameter = expanded_status$events_per_parameter[[1]],
      same_rows = same_rows, df_difference = df_difference,
      computationally_usable = computationally_usable,
      comparison_status = comparison_status,
      interaction_estimate = if (computationally_usable) interaction_estimate else NA_real_,
      interaction_standard_error = if (computationally_usable) interaction_standard_error else NA_real_,
      interaction_wald_p = if (computationally_usable) interaction_wald_p else NA_real_,
      interaction_lrt_statistic = if (computationally_usable) as.numeric(comparison$Deviance[2]) else NA_real_,
      interaction_lrt_p = if (computationally_usable) as.numeric(comparison$`Pr(>Chi)`[2]) else NA_real_
    )
  })
  
  list(
    family = family_name,
    normalization_spec = normalization_spec,
    normalization_variable = registry_row$normalization_variable[[1]],
    normalization_label = registry_row$normalization_label[[1]],
    normalization_kind = registry_row$normalization_kind[[1]],
    data = analysis_data,
    formulas = model_formulas,
    models = models,
    model_diagnostics = model_diagnostics,
    row_validation = row_validation,
    interaction_results = interaction_results,
    mean_log_mechanic = prepared$mean_log_mechanic,
    mean_log_pfvc = prepared$mean_log_pfvc,
    age_reference = prepared$age_reference
  )
}

single_modifier_fit_grid <- crossing(
  family = exps$family,
  normalization_spec = normalization_adjustment_specs$normalization_spec
)

message("Script 05 | fitting balanced single-modifier models")
single_modifier_fits <- pmap(
  single_modifier_fit_grid,
  function(family, normalization_spec) {
    fit_single_modifier_models(
      family_name = family,
      normalization_spec = normalization_spec
    )
  }
)

names(single_modifier_fits) <- paste(
  single_modifier_fit_grid$family,
  single_modifier_fit_grid$normalization_spec,
  sep = " | "
)

single_modifier_model_qc <- purrr::map_dfr(single_modifier_fits, "model_diagnostics")
single_modifier_row_qc <- purrr::map_dfr(single_modifier_fits, "row_validation")
single_modifier_diagnostics <- purrr::map_dfr(single_modifier_fits, "interaction_results") %>%
  mutate(site = site_name, .before = 1)

single_modifier_cross_specification_row_qc <- tibble(family = exps$family) %>%
  mutate(
    common_rows_preserved = purrr::map_lgl(family, function(family_name) {
      fit_indices <- which(single_modifier_fit_grid$family == family_name)
      expected_rows <- rownames(ladder_prepared_data[[family_name]]$data)
      
      length(fit_indices) == nrow(normalization_adjustment_specs) &&
        all(purrr::map_lgl(single_modifier_fits[fit_indices], function(fit) {
          all(purrr::map_lgl(
            fit$models,
            ~ identical(rownames(model.frame(.x)), expected_rows)
          ))
        }))
    })
  )

single_modifier_checks <- list(
  single_modifier_expected_sets_fitted =
    length(single_modifier_fits) == nrow(single_modifier_fit_grid),
  single_modifier_fit_grid_keys_unique = !anyDuplicated(
    single_modifier_fit_grid[c("family", "normalization_spec")]
  ),
  single_modifier_fit_names_unique = !anyDuplicated(names(single_modifier_fits)),
  single_modifier_fit_metadata_matches_grid = all(purrr::map2_lgl(
    single_modifier_fits,
    seq_len(nrow(single_modifier_fit_grid)),
    function(fit, fit_index) {
      identical(fit$family, single_modifier_fit_grid$family[[fit_index]]) &&
        identical(fit$normalization_spec, single_modifier_fit_grid$normalization_spec[[fit_index]])
    }
  )),
  single_modifier_model_keys_unique = !anyDuplicated(
    single_modifier_model_qc[c("family", "normalization_spec", "model")]
  ),
  single_modifier_interaction_keys_unique = !anyDuplicated(
    single_modifier_diagnostics[c("site", "family", "normalization_spec", "modifier")]
  ),
  single_modifier_three_models_per_set = single_modifier_model_qc %>%
    count(family, normalization_spec) %>%
    summarise(passed = n() == nrow(single_modifier_fit_grid) && all(n == 3L)) %>%
    pull(passed),
  single_modifier_expected_number_of_models = nrow(single_modifier_model_qc) ==
    nrow(single_modifier_fit_grid) * length(single_modifier_formulas),
  single_modifier_row_validation_matches_grid = nrow(single_modifier_row_qc) ==
    nrow(single_modifier_fit_grid),
  single_modifier_all_alternatives_present = setequal(
    unique(single_modifier_model_qc$normalization_spec),
    normalization_adjustment_specs$normalization_spec
  ),
  single_modifier_common_rows_across_alternatives = all(
    single_modifier_cross_specification_row_qc$common_rows_preserved
  ),
  single_modifier_model_sets_use_one_sample_size = all(
    single_modifier_row_qc$n_distinct_model_sample_sizes == 1L
  ),
  single_modifier_models_use_expected_n = all(single_modifier_row_qc$all_models_use_expected_n),
  single_modifier_models_within_mechanic_use_same_rows = all(
    single_modifier_row_qc$all_models_use_same_rows
  ),
  single_modifier_model_matrices_full_rank = all(single_modifier_model_qc$full_rank),
  single_modifier_models_have_no_aliased_coefficients = all(
    single_modifier_model_qc$n_aliased_coefficients == 0L
  ),
  single_modifier_coefficients_finite = all(single_modifier_model_qc$coefficients_finite),
  single_modifier_standard_errors_finite = all(single_modifier_model_qc$standard_errors_finite),
  single_modifier_fitted_probability_ranges_valid = all(
    is.finite(single_modifier_model_qc$minimum_fitted_probability) &
      is.finite(single_modifier_model_qc$maximum_fitted_probability) &
      single_modifier_model_qc$minimum_fitted_probability >= 0 &
      single_modifier_model_qc$maximum_fitted_probability <= 1
  ),
  single_modifier_two_interaction_tests_per_set = single_modifier_diagnostics %>%
    count(family, normalization_spec) %>%
    summarise(passed = n() == nrow(single_modifier_fit_grid) && all(n == 2L)) %>%
    pull(passed),
  single_modifier_expected_number_of_interaction_tests = nrow(single_modifier_diagnostics) ==
    2L * nrow(single_modifier_fit_grid),
  single_modifier_interaction_models_add_one_parameter = all(
    single_modifier_diagnostics$df_difference == 1L
  ),
  single_modifier_interaction_comparisons_use_same_rows = all(
    single_modifier_diagnostics$same_rows
  ),
  single_modifier_usable_interaction_results_valid = check_subset(
    single_modifier_diagnostics,
    single_modifier_diagnostics$computationally_usable,
    function(data) {
      all(
        is.finite(data$interaction_estimate) &
          is.finite(data$interaction_standard_error) &
          data$interaction_standard_error > 0 &
          data$interaction_wald_p >= 0 &
          data$interaction_wald_p <= 1 &
          data$interaction_lrt_p >= 0 &
          data$interaction_lrt_p <= 1
      )
    }
  ),
  single_modifier_unusable_interaction_results_suppressed = check_subset(
    single_modifier_diagnostics,
    !single_modifier_diagnostics$computationally_usable,
    function(data) {
      all(
        is.na(data$interaction_estimate) &
          is.na(data$interaction_standard_error) &
          is.na(data$interaction_wald_p) &
          is.na(data$interaction_lrt_statistic) &
          is.na(data$interaction_lrt_p)
      )
    }
  )
)

if (!isTRUE(all(unlist(single_modifier_checks, use.names = FALSE)))) {
  print(
    single_modifier_model_qc %>%
      select(
        family, normalization_spec, normalization_label, model, n_obs, n_events, n_parameters_excluding_intercept,
        events_per_parameter, converged, boundary, full_rank,
        n_extreme_fitted_probabilities, computationally_usable, model_status
      ),
    n = Inf,
    width = Inf
  )
  print(single_modifier_diagnostics, n = Inf, width = Inf)
  print(single_modifier_row_qc, n = Inf, width = Inf)
  print(single_modifier_cross_specification_row_qc, n = Inf, width = Inf)
  print(single_modifier_checks)
}

stopifnot(isTRUE(all(unlist(single_modifier_checks, use.names = FALSE))))

write_csv(single_modifier_diagnostics, single_modifier_output_path)

single_modifier_output_qc <- tibble(
  output = "Balanced single-modifier diagnostics",
  path = single_modifier_output_path,
  exists = file.exists(single_modifier_output_path),
  size_bytes = if (file.exists(single_modifier_output_path)) {
    file.info(single_modifier_output_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
) %>%
  mutate(nonempty = exists & is.finite(size_bytes) & size_bytes > 0)

single_modifier_output_checks <- list(
  single_modifier_output_exists = single_modifier_output_qc$exists,
  single_modifier_output_nonempty = single_modifier_output_qc$nonempty
)

stopifnot(isTRUE(all(unlist(single_modifier_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$mechanic_modification$single_modifier_models <- list(
  interpretation = paste0(
    "Age and PFVC modification are evaluated in separate models with identical main effects within each alternative ",
    "adjustment specification. All alternatives share the prepared mechanic-family sample and centering constants. ",
    "Computationally stable sparse-event models are retained as exploratory; unstable interaction tests are suppressed."
  ),
  fit_grid = single_modifier_fit_grid,
  normalization_specifications = normalization_adjustment_specs,
  formulas = purrr::map(single_modifier_fits, function(fit) {
    purrr::map_chr(fit$formulas, ~ paste(deparse(.x), collapse = " "))
  }),
  model_diagnostics = single_modifier_model_qc,
  row_validation = single_modifier_row_qc,
  cross_specification_row_validation = single_modifier_cross_specification_row_qc,
  interaction_results = single_modifier_diagnostics,
  output = single_modifier_output_qc
)

normalization_analysis_diagnostics$warnings$single_modifier_models <- single_modifier_model_qc %>%
  filter(model_status != "Eligible for conventional interpretation")

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  single_modifier_checks,
  single_modifier_output_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# Balanced single-modifier effect curves
# =============================================================================
# Curves show the 60-day mortality OR per 10% higher mechanic across continuous age or PFVC. A curve is displayed only
# when the interaction model is computationally stable and has at least the prespecified number of deaths per fitted
# parameter. This threshold is a conservative visualization heuristic, not an absolute model-validity boundary.

MIN_MODIFIER_CURVE_EVENTS_PER_PARAMETER <- 10
MECHANIC_CURVE_LOG_CONTRAST <- log(1.10)

modifier_curve_csv_path <- file.path(table_dir, paste0("norm_age_interaction_slopes_", site_name, ".csv"))
modifier_curve_pdf_path <- file.path(figure_dir, paste0("norm_age_interaction_slopes_", site_name, ".pdf"))

make_single_modifier_curve <- function(fit, family_name, modifier) {
  model_name <- if (modifier == "Age") "age_modifier" else "pfvc_modifier"
  fitted_model <- fit$models[[model_name]]
  
  model_status <- fit$model_diagnostics %>%
    filter(model == model_name)
  
  interaction_status <- fit$interaction_results %>%
    filter(.data$modifier == !!modifier)
  
  stopifnot(nrow(model_status) == 1L, nrow(interaction_status) == 1L)
  
  curve_usable <- interaction_status$computationally_usable[[1]] &&
    is.finite(model_status$events_per_parameter[[1]]) &&
    model_status$events_per_parameter[[1]] >= MIN_MODIFIER_CURVE_EVENTS_PER_PARAMETER
  
  curve_status <- if (!interaction_status$computationally_usable[[1]]) {
    "Curve suppressed: interaction model is computationally unstable."
  } else if (model_status$events_per_parameter[[1]] < MIN_MODIFIER_CURVE_EVENTS_PER_PARAMETER) {
    paste0(
      "Curve suppressed: ", round(model_status$events_per_parameter[[1]], 2),
      " deaths per fitted parameter; threshold = ", MIN_MODIFIER_CURVE_EVENTS_PER_PARAMETER, "."
    )
  } else {
    "Curve displayed: prespecified stability criteria passed."
  }
  
  if (modifier == "Age") {
    centered_grid <- seq(
      quantile(fit$data$age_c, 0.05, names = FALSE),
      quantile(fit$data$age_c, 0.95, names = FALSE),
      length.out = 100
    )
    
    display_grid <- centered_grid + fit$age_reference
    interaction_candidates <- c("mech_c:age_c", "age_c:mech_c")
    axis_label <- "Age at admission (years)"
  } else {
    centered_grid <- seq(
      quantile(fit$data$pfvc_c, 0.05, names = FALSE),
      quantile(fit$data$pfvc_c, 0.95, names = FALSE),
      length.out = 100
    )
    
    display_grid <- exp(centered_grid + fit$mean_log_pfvc)
    interaction_candidates <- c("mech_c:pfvc_c", "pfvc_c:mech_c")
    axis_label <- "Predicted FVC (L)"
  }
  
  if (!curve_usable) {
    return(
      tibble(
        site = site_name, family = family_name, modifier = modifier,
        normalization_spec = fit$normalization_spec,
        normalization_variable = fit$normalization_variable,
        normalization_label = fit$normalization_label,
        normalization_kind = fit$normalization_kind,
        x = display_grid, axis = axis_label,
        ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
        n_obs = nobs(fitted_model),
        n_events = model_status$n_events[[1]],
        n_parameters_excluding_intercept = model_status$n_parameters_excluding_intercept[[1]],
        events_per_parameter = model_status$events_per_parameter[[1]],
        estimable = FALSE, curve_status = curve_status
      )
    )
  }
  
  coefficient_estimates <- coef(fitted_model)
  covariance_matrix <- vcov(fitted_model)
  interaction_term <- intersect(interaction_candidates, names(coefficient_estimates))
  
  if (length(interaction_term) != 1L) {
    stop("Expected one ", modifier, " interaction coefficient for ", family_name, ".")
  }
  
  coefficient_terms <- c("mech_c", interaction_term)
  coefficient_subset <- coefficient_estimates[coefficient_terms]
  covariance_subset <- covariance_matrix[coefficient_terms, coefficient_terms, drop = FALSE]
  
  curve_results <- purrr::map_dfr(centered_grid, function(modifier_value) {
    contrast_weights <- c(1, modifier_value) * MECHANIC_CURVE_LOG_CONTRAST
    log_odds_ratio <- sum(contrast_weights * coefficient_subset)
    log_odds_ratio_variance <- as.numeric(
      crossprod(contrast_weights, covariance_subset %*% contrast_weights)
    )
    
    log_odds_ratio_standard_error <- sqrt(max(log_odds_ratio_variance, 0))
    
    tibble(
      ratio = exp(log_odds_ratio),
      conf_low = exp(log_odds_ratio - 1.96 * log_odds_ratio_standard_error),
      conf_high = exp(log_odds_ratio + 1.96 * log_odds_ratio_standard_error)
    )
  })
  
  bind_cols(
    tibble(
      site = site_name, family = family_name, modifier = modifier,
      normalization_spec = fit$normalization_spec,
      normalization_variable = fit$normalization_variable,
      normalization_label = fit$normalization_label,
      normalization_kind = fit$normalization_kind,
      x = display_grid, axis = axis_label
    ),
    curve_results
  ) %>%
    mutate(
      n_obs = nobs(fitted_model),
      n_events = model_status$n_events[[1]],
      n_parameters_excluding_intercept = model_status$n_parameters_excluding_intercept[[1]],
      events_per_parameter = model_status$events_per_parameter[[1]],
      estimable = TRUE,
      curve_status = curve_status
    )
}

modifier_curve_grid <- crossing(
  single_modifier_fit_grid,
  modifier = c("Age", "PFVC")
)

modifier_curve_data <- modifier_curve_grid %>%
  pmap_dfr(function(family, normalization_spec, modifier) {
    fit_name <- paste(family, normalization_spec, sep = " | ")
    fit <- single_modifier_fits[[fit_name]]
    
    if (is.null(fit)) {
      stop("The expected single-modifier fit was not found: ", fit_name)
    }
    
    stopifnot(
      identical(fit$family, family),
      identical(fit$normalization_spec, normalization_spec)
    )
    
    make_single_modifier_curve(
      fit = fit,
      family_name = family,
      modifier = modifier
    )
  }) %>%
  mutate(
    family = factor(
      family,
      levels = c("Estimated driving pressure", "Elastance", "Mechanical power")
    ),
    normalization_label = factor(
      normalization_label,
      levels = normalization_adjustment_specs$normalization_label
    ),
    modifier = factor(modifier, levels = c("Age", "PFVC"))
  )

modifier_curve_status <- modifier_curve_data %>%
  distinct(
    family, normalization_spec, normalization_variable, normalization_label, normalization_kind,
    modifier, n_obs, n_events, n_parameters_excluding_intercept,
    events_per_parameter, estimable, curve_status
  )

write_csv(modifier_curve_data, modifier_curve_csv_path)

modifier_curve_labels <- modifier_curve_data %>%
  filter(!estimable) %>%
  group_by(family, normalization_spec, normalization_label, modifier, curve_status) %>%
  summarise(x = median(x), .groups = "drop") %>%
  mutate(
    y = 1,
    label = str_wrap(curve_status, width = 34)
  )

age_curve_data <- modifier_curve_data %>%
  filter(modifier == "Age")

pfvc_curve_data <- modifier_curve_data %>%
  filter(modifier == "PFVC")

p_age_modifier <- ggplot(age_curve_data, aes(x = x, y = ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_ribbon(
    data = age_curve_data %>% filter(estimable),
    aes(ymin = conf_low, ymax = conf_high),
    fill = "#0072B2", alpha = 0.15
  ) +
  geom_line(
    data = age_curve_data %>% filter(estimable),
    color = "#0072B2", linewidth = 0.9
  ) +
  geom_text(
    data = modifier_curve_labels %>% filter(modifier == "Age"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, color = "#B2182B", size = 3.1, lineheight = 1.05
  ) +
  facet_grid(normalization_label ~ family, drop = FALSE) +
  scale_x_continuous(breaks = seq(3, 17, by = 2)) +
  labs(
    title = "Mechanic association with mortality across age",
    subtitle = str_wrap(
      "Age-modifier models include PFVC as a main effect. Curves span the 5th–95th percentile of observed age.",
      width = 100
    ),
    x = "Age at admission (years)",
    y = "60-day mortality OR per 10% higher mechanic"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 8)),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
  )

p_pfvc_modifier <- ggplot(pfvc_curve_data, aes(x = x, y = ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_ribbon(
    data = pfvc_curve_data %>% filter(estimable),
    aes(ymin = conf_low, ymax = conf_high),
    fill = "#E69F00", alpha = 0.15
  ) +
  geom_line(
    data = pfvc_curve_data %>% filter(estimable),
    color = "#E69F00", linewidth = 0.9
  ) +
  geom_text(
    data = modifier_curve_labels %>% filter(modifier == "PFVC"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, color = "#B2182B", size = 3.1, lineheight = 1.05
  ) +
  facet_grid(normalization_label ~ family, drop = FALSE) +
  labs(
    title = "Mechanic association with mortality across predicted FVC",
    subtitle = str_wrap(
      "PFVC-modifier models include continuous age as a main effect. Curves span the 5th–95th percentile of PFVC.",
      width = 100
    ),
    x = "Predicted FVC (L)",
    y = "60-day mortality OR per 10% higher mechanic",
    caption = str_wrap(
      paste0(
        "Curves are displayed only when the interaction model is computationally stable and has at least ",
        MIN_MODIFIER_CURVE_EVENTS_PER_PARAMETER,
        " deaths per fitted parameter. This is a conservative visualization heuristic, not an absolute validity boundary."
      ),
      width = 115
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 8)),
    plot.caption = element_text(hjust = 0, margin = margin(t = 8)),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
  )

modifier_curve_figure <- p_age_modifier / p_pfvc_modifier +
  patchwork::plot_layout(heights = c(1, 1.1)) +
  patchwork::plot_annotation(
    title = "Modification of respiratory-mechanic associations with 60-day mortality",
    subtitle = str_wrap(
      paste0(
        site_name,
        " — balanced single-modifier models. Suppressed panels remain visible so the output structure is consistent ",
        "across sites."
      ),
      width = 110
    )
  )

ggsave(
  filename = modifier_curve_pdf_path,
  plot = modifier_curve_figure,
  width = 16,
  height = 32,
  units = "in",
  limitsize = TRUE
)

modifier_curve_output_qc <- tibble(
  output = c("Single-modifier curve data", "Single-modifier curve figure"),
  path = c(modifier_curve_csv_path, modifier_curve_pdf_path)
) %>%
  mutate(
    exists = file.exists(path),
    size_bytes = purrr::map_dbl(path, function(output_path) {
      if (file.exists(output_path)) file.info(output_path, extra_cols = FALSE)$size[[1]] else NA_real_
    }),
    nonempty = exists & is.finite(size_bytes) & size_bytes > 0
  )

expected_modifier_curve_keys <- modifier_curve_grid %>%
  distinct(family, normalization_spec, modifier) %>%
  mutate(
    family = as.character(family),
    modifier = as.character(modifier)
  )

observed_modifier_curve_keys <- modifier_curve_status %>%
  transmute(
    family = as.character(family),
    normalization_spec,
    modifier = as.character(modifier)
  )

suppressed_modifier_curve_keys <- modifier_curve_status %>%
  filter(!estimable) %>%
  transmute(
    family = as.character(family),
    normalization_spec,
    modifier = as.character(modifier)
  )

modifier_annotation_keys <- modifier_curve_labels %>%
  transmute(
    family = as.character(family),
    normalization_spec,
    modifier = as.character(modifier)
  )

modifier_curve_checks <- list(
  expected_curve_panels_created = nrow(modifier_curve_status) == nrow(modifier_curve_grid),
  curve_status_keys_unique = !anyDuplicated(
    modifier_curve_status[c("family", "normalization_spec", "modifier")]
  ),
  curve_grid_rows_unique_within_panel = !anyDuplicated(
    modifier_curve_data[c("family", "normalization_spec", "modifier", "x")]
  ),
  one_hundred_rows_per_panel = modifier_curve_data %>%
    count(family, normalization_spec, modifier) %>%
    summarise(
      passed = n() == nrow(modifier_curve_grid) && all(n == 100L)
    ) %>%
    pull(passed),
  curve_data_has_expected_total_rows = nrow(modifier_curve_data) == 100L * nrow(modifier_curve_grid),
  all_expected_curve_panels_present = nrow(anti_join(
    expected_modifier_curve_keys,
    observed_modifier_curve_keys,
    by = c("family", "normalization_spec", "modifier")
  )) == 0L,
  no_unexpected_curve_panels_present = nrow(anti_join(
    observed_modifier_curve_keys,
    expected_modifier_curve_keys,
    by = c("family", "normalization_spec", "modifier")
  )) == 0L,
  curve_status_complete = !any(is.na(modifier_curve_status$estimable)),
  curve_alternative_labels_complete = !any(is.na(modifier_curve_data$normalization_label)),
  curve_grid_values_finite = all(is.finite(modifier_curve_data$x)),
  estimable_curves_have_finite_values = check_subset(
    modifier_curve_data,
    modifier_curve_data$estimable,
    function(data) {
      all(is.finite(data$ratio) & is.finite(data$conf_low) & is.finite(data$conf_high))
    }
  ),
  estimable_confidence_intervals_ordered = check_subset(
    modifier_curve_data,
    modifier_curve_data$estimable,
    function(data) {
      all(data$conf_low <= data$ratio & data$ratio <= data$conf_high)
    }
  ),
  suppressed_curves_have_missing_estimates = check_subset(
    modifier_curve_data,
    !modifier_curve_data$estimable,
    function(data) {
      all(is.na(data$ratio) & is.na(data$conf_low) & is.na(data$conf_high))
    }
  ),
  suppression_annotation_keys_unique = !anyDuplicated(modifier_annotation_keys),
  suppressed_curves_have_annotations = nrow(anti_join(
    suppressed_modifier_curve_keys,
    modifier_annotation_keys,
    by = c("family", "normalization_spec", "modifier")
  )) == 0L,
  modifier_curve_no_unexpected_suppression_annotations = nrow(anti_join(
    modifier_annotation_keys,
    suppressed_modifier_curve_keys,
    by = c("family", "normalization_spec", "modifier")
  )) == 0L,
  modifier_curve_outputs_exist = all(modifier_curve_output_qc$exists),
  modifier_curve_outputs_nonempty = all(modifier_curve_output_qc$nonempty)
)

if (!isTRUE(all(unlist(modifier_curve_checks, use.names = FALSE)))) {
  print(modifier_curve_status, n = Inf, width = Inf)
  print(modifier_curve_output_qc, n = Inf, width = Inf)
  print(modifier_curve_checks)
}

stopifnot(isTRUE(all(unlist(modifier_curve_checks, use.names = FALSE))))

normalization_analysis_diagnostics$mechanic_modification$single_modifier_curves <- list(
  display_criteria = list(
    minimum_events_per_fitted_parameter = MIN_MODIFIER_CURVE_EVENTS_PER_PARAMETER,
    requires_computationally_usable_interaction_model = TRUE,
    mechanic_contrast = "10% increase",
    grid_range = "5th to 95th percentile of the observed modifier",
    note = paste0(
      "The events-per-parameter threshold is a conservative visualization heuristic, ",
      "not an absolute model-validity boundary."
    )
  ),
  fit_grid = single_modifier_fit_grid,
  curve_grid = modifier_curve_grid,
  normalization_specifications = normalization_adjustment_specs,
  status = modifier_curve_status,
  curve_data = modifier_curve_data,
  suppression_annotations = modifier_curve_labels,
  layout = list(
    row_variable = "normalization_label",
    column_variable = "family",
    plot_blocks = c("Age modifier", "PFVC modifier"),
    width_inches = 16,
    height_inches = 32,
    visual_review_status = "Not established by runtime checks"
  ),
  output = modifier_curve_output_qc
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  modifier_curve_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)


# =============================================================================
# PART 2c -- Incremental information from PBW and PFVC
# =============================================================================
# Separate mechanic-plus-PBW and mechanic-plus-PFVC models are nonnested. The encompassing model includes both size
# references, allowing two nested one-parameter tests:
#
#   1. PFVC beyond PBW: mechanic + PBW versus mechanic + PBW + PFVC
#   2. PBW beyond PFVC: mechanic + PFVC versus mechanic + PBW + PFVC
#
# These tests evaluate incremental model information conditional on the specified covariates. A nonsignificant test does
# not establish equivalence, absence of information, or causal superiority. PBW and PFVC are strongly related predicted
# quantities, so collinearity diagnostics are retained.

encompassing_csv_path <- file.path(table_dir, paste0("norm_encompassing_", site_name, ".csv"))
encompassing_pdf_path <- file.path(figure_dir, paste0("norm_encompassing_", site_name, ".pdf"))

encompassing_common_vars <- unique(c(
  ".source_row", "hospitalization_id", "mortality_event_60", "pbw", "pfvc",
  normalization_adjustment_vars,
  "sofa_total", "sf10", "bmi", "age", "sex_category", "race_category"
))

prepare_encompassing_data <- function(family_name, mechanic_variable, data_name) {
  source_data <- data_map[[data_name]]
  required_vars <- unique(c(encompassing_common_vars, mechanic_variable))
  missing_vars <- setdiff(required_vars, names(source_data))
  
  if (length(missing_vars) > 0L) {
    stop("The ", family_name, " encompassing analysis is missing variable(s): ", paste(missing_vars, collapse = ", "))
  }
  
  analysis_data <- source_data %>%
    select(all_of(required_vars)) %>%
    drop_na(all_of(required_vars)) %>%
    filter(
      is.finite(.data[[mechanic_variable]]), .data[[mechanic_variable]] > 0,
      is.finite(pbw), pbw > 0,
      is.finite(pfvc), pfvc > 0
    )
  
  if (nrow(analysis_data) == 0L) {
    stop("No complete observations are available for the ", family_name, " encompassing analysis.")
  }
  
  analysis_data
}

encompassing_analysis_data <- pmap(
  exps,
  function(family, mx, data_name) {
    prepare_encompassing_data(
      family_name = family,
      mechanic_variable = mx,
      data_name = data_name
    )
  }
)

names(encompassing_analysis_data) <- exps$family

encompassing_sample_missingness_qc <- purrr::imap_dfr(
  encompassing_analysis_data,
  function(analysis_data, family_name) {
    specification_index <- match(family_name, exps$family)
    source_data <- data_map[[exps$data_name[[specification_index]]]]
    required_vars <- unique(c(encompassing_common_vars, exps$mx[[specification_index]]))
    
    purrr::map_dfr(required_vars, function(variable) {
      tibble(
        family = family_name,
        variable = variable,
        n_metric_frame = nrow(source_data),
        n_missing = sum(is.na(source_data[[variable]])),
        percent_missing = 100 * mean(is.na(source_data[[variable]]))
      )
    })
  }
)

encompassing_excluded_rows <- purrr::imap_dfr(
  encompassing_analysis_data,
  function(analysis_data, family_name) {
    specification_index <- match(family_name, exps$family)
    source_data <- data_map[[exps$data_name[[specification_index]]]]
    mechanic_variable <- exps$mx[[specification_index]]
    required_vars <- unique(c(encompassing_common_vars, mechanic_variable))
    
    excluded_indices <- which(!source_data$.source_row %in% analysis_data$.source_row)
    
    tibble(
      family = rep(family_name, length(excluded_indices)),
      source_row = source_data$.source_row[excluded_indices],
      hospitalization_id = source_data$hospitalization_id[excluded_indices],
      exclusion_reason = purrr::map_chr(excluded_indices, function(row_index) {
        missing_vars <- required_vars[purrr::map_lgl(
          required_vars,
          ~ is.na(source_data[[.x]][row_index])
        )]
        
        reasons <- character()
        
        if (length(missing_vars) > 0L) {
          reasons <- c(reasons, paste0("Missing: ", paste(missing_vars, collapse = ", ")))
        }
        
        for (variable in unique(c(mechanic_variable, "pbw", "pfvc"))) {
          value <- source_data[[variable]][[row_index]]
          
          if (!is.finite(value) || value <= 0) {
            reasons <- c(reasons, paste0(variable, " not finite and positive"))
          }
        }
        
        paste(reasons, collapse = "; ")
      })
    )
  }
)

encompassing_sample_qc <- purrr::imap_dfr(encompassing_analysis_data, function(data, family_name) {
  source_data_name <- exps$data_name[match(family_name, exps$family)]
  source_data <- data_map[[source_data_name]]
  
  tibble(
    family = family_name,
    n_metric_frame = nrow(source_data),
    n_analyzed = nrow(data),
    n_excluded = n_metric_frame - n_analyzed,
    percent_analyzed = 100 * n_analyzed / n_metric_frame,
    n_events = sum(data$mortality_event_60 == 1L),
    n_nonevents = sum(data$mortality_event_60 == 0L),
    event_percent = 100 * mean(data$mortality_event_60 == 1L)
  )
})

encompassing_sample_checks <- list(
  encompassing_sample_all_three_samples_created =
    identical(names(encompassing_analysis_data), exps$family),
  encompassing_sample_covers_all_alternatives = all(
    normalization_adjustment_vars %in% encompassing_common_vars
  ),
  encompassing_sample_exclusions_reconcile = all(purrr::map2_lgl(
    encompassing_sample_qc$family,
    encompassing_sample_qc$n_excluded,
    function(family_name, expected_n) {
      sum(encompassing_excluded_rows$family == family_name) == expected_n
    }
  )),
  encompassing_sample_exclusion_keys_unique = !anyDuplicated(
    encompassing_excluded_rows[c("family", "source_row")]
  ),
  encompassing_sample_exclusion_reasons_nonempty = all(
    nzchar(encompassing_excluded_rows$exclusion_reason)
  ),
  encompassing_sample_samples_nonempty =
    all(purrr::map_int(encompassing_analysis_data, nrow) > 0L),
  encompassing_sample_hospitalization_ids_unique = all(purrr::map_lgl(
    encompassing_analysis_data,
    ~ !anyDuplicated(.x$hospitalization_id)
  )),
  encompassing_sample_samples_complete = all(purrr::map_lgl(
    encompassing_analysis_data,
    ~ all(complete.cases(.x))
  )),
  encompassing_sample_samples_have_events = all(encompassing_sample_qc$n_events > 0L),
  encompassing_sample_samples_have_nonevents = all(encompassing_sample_qc$n_nonevents > 0L),
  encompassing_sample_outcomes_binary = all(purrr::map_lgl(
    encompassing_analysis_data,
    ~ all(.x$mortality_event_60 %in% c(0L, 1L))
  ))
)

if (!isTRUE(all(unlist(encompassing_sample_checks, use.names = FALSE)))) {
  print(encompassing_sample_qc, n = Inf, width = Inf)
  print(encompassing_sample_checks)
}

stopifnot(isTRUE(all(unlist(encompassing_sample_checks, use.names = FALSE))))

normalization_analysis_diagnostics$encompassing_sample_configuration <- list(
  policy = paste0(
    "One common complete-case sample per mechanic family across all seven alternative adjustment variables, ",
    "both adjustment levels, and all encompassing model forms. ",
    "Existing positive/finite mechanic, PBW, and PFVC requirements are preserved. ",
    "Exclusions affect these outcome models only, not the upstream or descriptive cohorts."
  ),
  required_common_variables = encompassing_common_vars,
  sample_summary = encompassing_sample_qc,
  variable_missingness = encompassing_sample_missingness_qc,
  excluded_rows = encompassing_excluded_rows,
  source_rows = purrr::map(
    encompassing_analysis_data,
    ~ .x %>% select(.source_row, hospitalization_id)
  ),
  checks = encompassing_sample_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

summarise_encompassing_model <- function(
    fitted_model, family_name, adjustment_id, model_name, normalization_spec
) {
  registry_row <- normalization_adjustment_specs %>%
    filter(.data$normalization_spec == .env$normalization_spec)
  
  stopifnot(nrow(registry_row) == 1L)
  
  outcome <- model.response(model.frame(fitted_model))
  model_matrix <- model.matrix(fitted_model)
  coefficient_table <- summary(fitted_model)$coefficients
  aliased_coefficients <- is.na(coef(fitted_model))
  n_events <- sum(outcome == 1L)
  n_parameters_excluding_intercept <- fitted_model$rank - 1L
  
  tibble(
    family = family_name,
    normalization_spec = normalization_spec,
    normalization_variable = registry_row$normalization_variable[[1]],
    normalization_label = registry_row$normalization_label[[1]],
    normalization_kind = registry_row$normalization_kind[[1]],
    adjustment = adjustment_id,
    adjusted = unname(adjustment_labels[[adjustment_id]]),
    model = model_name,
    formula = paste(deparse(formula(fitted_model)), collapse = " "),
    n_obs = nobs(fitted_model),
    n_events = n_events,
    n_nonevents = sum(outcome == 0L),
    n_parameters_excluding_intercept = n_parameters_excluding_intercept,
    events_per_parameter = n_events / n_parameters_excluding_intercept,
    converged = isTRUE(fitted_model$converged),
    boundary = isTRUE(fitted_model$boundary),
    iterations = fitted_model$iter,
    model_matrix_columns = ncol(model_matrix),
    model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    n_aliased_coefficients = sum(aliased_coefficients),
    coefficients_finite = all(is.finite(coef(fitted_model)[!aliased_coefficients])),
    standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
    minimum_fitted_probability = min(fitted(fitted_model)),
    maximum_fitted_probability = max(fitted(fitted_model)),
    n_extreme_fitted_probabilities = sum(
      fitted(fitted_model) < EXTREME_FITTED_PROBABILITY_THRESHOLD |
        fitted(fitted_model) > 1 - EXTREME_FITTED_PROBABILITY_THRESHOLD
    ),
    log_likelihood = as.numeric(logLik(fitted_model)),
    aic = AIC(fitted_model)
  ) %>%
    mutate(
      computationally_usable = converged & !boundary & full_rank & n_aliased_coefficients == 0L &
        coefficients_finite & standard_errors_finite & n_extreme_fitted_probabilities == 0L,
      sparse_events = events_per_parameter < MIN_PROGNOSTIC_EVENTS_PER_PARAMETER,
      highly_sparse_events = events_per_parameter < HIGHLY_SPARSE_PROGNOSTIC_EVENTS_PER_PARAMETER,
      model_status = case_when(
        !converged ~ "Suppressed: model did not converge",
        boundary ~ "Suppressed: boundary fit",
        !full_rank | n_aliased_coefficients > 0L ~ "Suppressed: rank-deficient model",
        !coefficients_finite | !standard_errors_finite ~ "Suppressed: nonfinite inference",
        n_extreme_fitted_probabilities > 0L ~ "Suppressed: possible separation",
        highly_sparse_events ~ "Exploratory: highly sparse events",
        sparse_events ~ "Exploratory: sparse events",
        TRUE ~ "Eligible for conventional interpretation"
      )
    )
}

fit_encompassing_models <- function(family_name, mechanic_variable, adjustment_id, normalization_spec) {
  analysis_data <- encompassing_analysis_data[[family_name]]
  
  registry_row <- normalization_adjustment_specs %>%
    filter(.data$normalization_spec == .env$normalization_spec)
  
  stopifnot(nrow(registry_row) == 1L)
  
  covariates <- get_normalization_adjustment(normalization_spec, adjustment_id)
  
  model_formulas <- list(
    pbw = as.formula(
      paste("mortality_event_60 ~ log(", mechanic_variable, ") + log(pbw) +", covariates)
    ),
    pfvc = as.formula(
      paste("mortality_event_60 ~ log(", mechanic_variable, ") + log(pfvc) +", covariates)
    ),
    both = as.formula(
      paste("mortality_event_60 ~ log(", mechanic_variable, ") + log(pbw) + log(pfvc) +", covariates)
    )
  )
  
  models <- purrr::map(model_formulas, function(model_formula) {
    glm(
      model_formula, data = analysis_data, family = binomial,
      x = TRUE, y = TRUE, model = TRUE
    )
  })
  
  model_diagnostics <- purrr::imap_dfr(models, function(fitted_model, model_name) {
    summarise_encompassing_model(
      fitted_model = fitted_model,
      family_name = family_name,
      adjustment_id = adjustment_id,
      model_name = model_name,
      normalization_spec = normalization_spec
    )
  })
  
  reference_rows <- rownames(model.frame(models$both))
  
  row_validation <- tibble(
    family = family_name,
    normalization_spec = normalization_spec,
    adjustment = adjustment_id,
    expected_n = nrow(analysis_data),
    n_distinct_model_sample_sizes = n_distinct(purrr::map_int(models, nobs)),
    all_models_use_expected_n = all(purrr::map_int(models, nobs) == nrow(analysis_data)),
    all_models_use_same_rows = all(purrr::map_lgl(
      models,
      ~ identical(rownames(model.frame(.x)), reference_rows)
    )),
    pfvc_beyond_pbw_df_difference = attr(logLik(models$both), "df") - attr(logLik(models$pbw), "df"),
    pbw_beyond_pfvc_df_difference = attr(logLik(models$both), "df") - attr(logLik(models$pfvc), "df")
  )
  
  comparison_specs <- tribble(
    ~direction, ~reduced_model,
    "PFVC beyond PBW", "pbw",
    "PBW beyond PFVC", "pfvc"
  )
  
  comparisons <- pmap_dfr(comparison_specs, function(direction, reduced_model) {
    reduced_fit <- models[[reduced_model]]
    expanded_fit <- models$both
    
    reduced_status <- model_diagnostics %>%
      filter(model == reduced_model)
    
    expanded_status <- model_diagnostics %>%
      filter(model == "both")
    
    stopifnot(nrow(reduced_status) == 1L, nrow(expanded_status) == 1L)
    
    same_rows <- identical(
      rownames(model.frame(reduced_fit)),
      rownames(model.frame(expanded_fit))
    )
    
    df_difference <- attr(logLik(expanded_fit), "df") - attr(logLik(reduced_fit), "df")
    
    computationally_usable <- reduced_status$computationally_usable[[1]] &&
      expanded_status$computationally_usable[[1]] &&
      same_rows &&
      df_difference == 1L
    
    comparison <- if (computationally_usable) {
      anova(reduced_fit, expanded_fit, test = "LRT")
    } else {
      NULL
    }
    
    comparison_status <- if (!computationally_usable) {
      "Suppressed: at least one model is computationally unstable"
    } else if (
      reduced_status$highly_sparse_events[[1]] ||
      expanded_status$highly_sparse_events[[1]]
    ) {
      "Exploratory: highly sparse events"
    } else if (
      reduced_status$sparse_events[[1]] ||
      expanded_status$sparse_events[[1]]
    ) {
      "Exploratory: sparse events"
    } else {
      "Eligible for conventional interpretation"
    }
    
    tibble(
      family = family_name,
      normalization_spec = normalization_spec,
      normalization_variable = registry_row$normalization_variable[[1]],
      normalization_label = registry_row$normalization_label[[1]],
      normalization_kind = registry_row$normalization_kind[[1]],
      adjustment = adjustment_id,
      adjusted = unname(adjustment_labels[[adjustment_id]]),
      direction = direction,
      reduced_model = reduced_model,
      expanded_model = "both",
      n_obs = nobs(reduced_fit),
      n_events = sum(model.response(model.frame(reduced_fit)) == 1L),
      same_rows = same_rows,
      df_difference = df_difference,
      computationally_usable = computationally_usable,
      comparison_status = comparison_status,
      lrt_statistic = if (computationally_usable) as.numeric(comparison$Deviance[2]) else NA_real_,
      lrt_df = if (computationally_usable) as.numeric(comparison$Df[2]) else NA_real_,
      lrt_p_value = if (computationally_usable) as.numeric(comparison$`Pr(>Chi)`[2]) else NA_real_,
      delta_aic_from_adding = if (computationally_usable) AIC(reduced_fit) - AIC(expanded_fit) else NA_real_
    )
  })
  
  both_model_usable <- model_diagnostics$computationally_usable[model_diagnostics$model == "both"]
  
  max_vif_both <- if (isTRUE(both_model_usable)) {
    vif_max(models$both)
  } else {
    NA_real_
  }
  
  list(
    family = family_name,
    mechanic_variable = mechanic_variable,
    normalization_spec = normalization_spec,
    normalization_variable = registry_row$normalization_variable[[1]],
    normalization_label = registry_row$normalization_label[[1]],
    normalization_kind = registry_row$normalization_kind[[1]],
    adjustment = adjustment_id,
    data = analysis_data,
    formulas = model_formulas,
    models = models,
    model_diagnostics = model_diagnostics,
    row_validation = row_validation,
    comparisons = comparisons,
    max_vif_both = max_vif_both
  )
}

encompassing_fit_grid <- crossing(
  family = exps$family,
  normalization_spec = normalization_adjustment_specs$normalization_spec,
  adjustment = names(adjustment_labels)
) %>%
  left_join(
    exps %>% select(family, mechanic_variable = mx),
    by = "family"
  )

message("Script 05 | fitting encompassing size-reference models")
encompassing_fits <- pmap(
  encompassing_fit_grid,
  function(family, normalization_spec, adjustment, mechanic_variable) {
    fit_encompassing_models(
      family_name = family,
      mechanic_variable = mechanic_variable,
      adjustment_id = adjustment,
      normalization_spec = normalization_spec
    )
  }
)

names(encompassing_fits) <- paste(
  encompassing_fit_grid$family,
  encompassing_fit_grid$normalization_spec,
  encompassing_fit_grid$adjustment,
  sep = " | "
)

encompassing_grid_checks <- list(
  encompassing_grid_has_expected_sets = nrow(encompassing_fit_grid) ==
    nrow(exps) * nrow(normalization_adjustment_specs) * length(adjustment_labels),
  encompassing_grid_keys_unique = !anyDuplicated(
    encompassing_fit_grid[c("family", "normalization_spec", "adjustment")]
  ),
  encompassing_grid_mechanic_variables_complete = !any(is.na(encompassing_fit_grid$mechanic_variable)),
  encompassing_fit_names_unique = !anyDuplicated(names(encompassing_fits)),
  encompassing_fit_count_matches_grid = length(encompassing_fits) == nrow(encompassing_fit_grid),
  encompassing_fit_metadata_matches_grid = all(purrr::map2_lgl(
    encompassing_fits,
    seq_len(nrow(encompassing_fit_grid)),
    function(fit, fit_index) {
      identical(fit$family, encompassing_fit_grid$family[[fit_index]]) &&
        identical(fit$normalization_spec, encompassing_fit_grid$normalization_spec[[fit_index]]) &&
        identical(fit$adjustment, encompassing_fit_grid$adjustment[[fit_index]]) &&
        identical(fit$mechanic_variable, encompassing_fit_grid$mechanic_variable[[fit_index]])
    }
  ))
)

stopifnot(isTRUE(all(unlist(encompassing_grid_checks, use.names = FALSE))))

normalization_analysis_diagnostics$encompassing_sample_configuration$fit_grid <- encompassing_fit_grid

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  encompassing_grid_checks
)

encompassing_model_qc <- purrr::map_dfr(encompassing_fits, "model_diagnostics")
encompassing_comparison_qc <- purrr::map_dfr(encompassing_fits, "comparisons")
encompassing_row_qc <- purrr::map_dfr(encompassing_fits, "row_validation")

encompassing_cross_specification_row_qc <- tibble(family = exps$family) %>%
  mutate(
    common_rows_preserved = purrr::map_lgl(family, function(family_name) {
      fit_indices <- which(encompassing_fit_grid$family == family_name)
      expected_rows <- rownames(encompassing_analysis_data[[family_name]])
      expected_sets <- nrow(normalization_adjustment_specs) * length(adjustment_labels)
      
      length(fit_indices) == expected_sets &&
        all(purrr::map_lgl(encompassing_fits[fit_indices], function(fit) {
          all(purrr::map_lgl(
            fit$models,
            ~ identical(rownames(model.frame(.x)), expected_rows)
          ))
        }))
    })
  )

encompassing_common_row_checks <- list(
  encompassing_common_row_validation_covers_all_families = setequal(
    encompassing_cross_specification_row_qc$family,
    exps$family
  ),
  encompassing_alternatives_and_adjustments_share_family_rows = all(
    encompassing_cross_specification_row_qc$common_rows_preserved
  )
)

stopifnot(isTRUE(all(unlist(encompassing_common_row_checks, use.names = FALSE))))

normalization_analysis_diagnostics$encompassing_sample_configuration$cross_specification_row_validation <-
  encompassing_cross_specification_row_qc

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  encompassing_common_row_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

encompassing_fit_checks <- list(
  encompassing_fit_expected_sets_fitted =
    length(encompassing_fits) == nrow(encompassing_fit_grid),
  encompassing_fit_model_keys_unique = !anyDuplicated(
    encompassing_model_qc[c("family", "normalization_spec", "adjustment", "model")]
  ),
  encompassing_fit_comparison_keys_unique = !anyDuplicated(
    encompassing_comparison_qc[c("family", "normalization_spec", "adjustment", "direction")]
  ),
  encompassing_fit_row_validation_matches_grid = nrow(encompassing_row_qc) ==
    nrow(encompassing_fit_grid),
  encompassing_fit_three_models_per_set = encompassing_model_qc %>%
    count(family, normalization_spec, adjustment) %>%
    summarise(passed = n() == nrow(encompassing_fit_grid) && all(n == 3L)) %>%
    pull(passed),
  encompassing_fit_expected_model_count = nrow(encompassing_model_qc) ==
    3L * nrow(encompassing_fit_grid),
  encompassing_fit_sets_use_one_sample_size = all(
    encompassing_row_qc$n_distinct_model_sample_sizes == 1L
  ),
  encompassing_fit_models_use_expected_n = all(encompassing_row_qc$all_models_use_expected_n),
  encompassing_fit_models_within_set_use_same_rows = all(
    encompassing_row_qc$all_models_use_same_rows
  ),
  encompassing_fit_nested_comparisons_add_one_parameter = all(
    encompassing_row_qc$pfvc_beyond_pbw_df_difference == 1L &
      encompassing_row_qc$pbw_beyond_pfvc_df_difference == 1L
  ),
  encompassing_fit_model_matrices_full_rank = all(encompassing_model_qc$full_rank),
  encompassing_fit_models_have_no_aliased_coefficients = all(
    encompassing_model_qc$n_aliased_coefficients == 0L
  ),
  encompassing_fit_coefficients_finite = all(encompassing_model_qc$coefficients_finite),
  encompassing_fit_standard_errors_finite = all(encompassing_model_qc$standard_errors_finite),
  encompassing_fit_probability_ranges_valid = all(
    is.finite(encompassing_model_qc$minimum_fitted_probability) &
      is.finite(encompassing_model_qc$maximum_fitted_probability) &
      encompassing_model_qc$minimum_fitted_probability >= 0 &
      encompassing_model_qc$maximum_fitted_probability <= 1
  ),
  encompassing_fit_two_comparisons_per_set = encompassing_comparison_qc %>%
    count(family, normalization_spec, adjustment) %>%
    summarise(passed = n() == nrow(encompassing_fit_grid) && all(n == 2L)) %>%
    pull(passed),
  encompassing_fit_expected_comparison_count = nrow(encompassing_comparison_qc) ==
    2L * nrow(encompassing_fit_grid),
  encompassing_fit_comparisons_use_same_rows = all(encompassing_comparison_qc$same_rows),
  encompassing_fit_usable_comparison_results_valid = check_subset(
    encompassing_comparison_qc,
    encompassing_comparison_qc$computationally_usable,
    function(data) {
      all(
        is.finite(data$lrt_statistic) &
          data$lrt_statistic >= 0 &
          is.finite(data$lrt_df) &
          data$lrt_df == 1 &
          is.finite(data$lrt_p_value) &
          data$lrt_p_value >= 0 &
          data$lrt_p_value <= 1 &
          is.finite(data$delta_aic_from_adding)
      )
    }
  ),
  encompassing_fit_unusable_comparison_results_suppressed = check_subset(
    encompassing_comparison_qc,
    !encompassing_comparison_qc$computationally_usable,
    function(data) {
      all(
        is.na(data$lrt_statistic) &
          is.na(data$lrt_df) &
          is.na(data$lrt_p_value) &
          is.na(data$delta_aic_from_adding)
      )
    }
  )
)

if (!isTRUE(all(unlist(encompassing_sample_checks, use.names = FALSE))) ||
    !isTRUE(all(unlist(encompassing_fit_checks, use.names = FALSE)))) {
  print(encompassing_sample_qc, n = Inf, width = Inf)
  print(encompassing_model_qc, n = Inf, width = Inf)
  print(encompassing_comparison_qc, n = Inf, width = Inf)
  print(encompassing_row_qc, n = Inf, width = Inf)
  print(encompassing_sample_checks)
  print(encompassing_fit_checks)
}

stopifnot(isTRUE(all(unlist(encompassing_sample_checks, use.names = FALSE))))
stopifnot(isTRUE(all(unlist(encompassing_fit_checks, use.names = FALSE))))

enc_tbl <- encompassing_fit_grid %>%
  pmap_dfr(function(family, normalization_spec, adjustment, mechanic_variable) {
    fit_name <- paste(family, normalization_spec, adjustment, sep = " | ")
    fit <- encompassing_fits[[fit_name]]
    
    if (is.null(fit)) {
      stop("The expected encompassing fit was not found: ", fit_name)
    }
    
    stopifnot(
      identical(fit$family, family),
      identical(fit$normalization_spec, normalization_spec),
      identical(fit$mechanic_variable, mechanic_variable)
    )
    
    pfvc_result <- fit$comparisons %>%
      filter(direction == "PFVC beyond PBW")
    
    pbw_result <- fit$comparisons %>%
      filter(direction == "PBW beyond PFVC")
    
    stopifnot(nrow(pfvc_result) == 1L, nrow(pbw_result) == 1L)
    
    overall_usable <- pfvc_result$computationally_usable[[1]] &&
      pbw_result$computationally_usable[[1]]
    
    comparison_status <- if (!overall_usable) {
      "Suppressed: computational instability or separation"
    } else if (
      str_detect(pfvc_result$comparison_status[[1]], fixed("highly sparse")) ||
      str_detect(pbw_result$comparison_status[[1]], fixed("highly sparse"))
    ) {
      "Exploratory: highly sparse events"
    } else if (
      str_detect(pfvc_result$comparison_status[[1]], fixed("sparse")) ||
      str_detect(pbw_result$comparison_status[[1]], fixed("sparse"))
    ) {
      "Exploratory: sparse events"
    } else {
      "Eligible for conventional interpretation"
    }
    
    tibble(
      site = site_name,
      family = family,
      normalization_spec = normalization_spec,
      normalization_variable = fit$normalization_variable,
      normalization_label = fit$normalization_label,
      normalization_kind = fit$normalization_kind,
      adjustment = adjustment,
      adjusted = unname(adjustment_labels[[adjustment]]),
      n = nrow(fit$data),
      n_events = sum(fit$data$mortality_event_60 == 1L),
      computationally_usable = overall_usable,
      comparison_status = comparison_status,
      pfvc_beyond_pbw_p = if (overall_usable) pfvc_result$lrt_p_value[[1]] else NA_real_,
      pbw_beyond_pfvc_p = if (overall_usable) pbw_result$lrt_p_value[[1]] else NA_real_,
      daic_pfvc_beyond_pbw = if (overall_usable) pfvc_result$delta_aic_from_adding[[1]] else NA_real_,
      daic_pbw_beyond_pfvc = if (overall_usable) pbw_result$delta_aic_from_adding[[1]] else NA_real_,
      max_vif_both = if (overall_usable) fit$max_vif_both else NA_real_
    )
  })

encompassing_summary_checks <- list(
  encompassing_summary_expected_rows_created = nrow(enc_tbl) == nrow(encompassing_fit_grid),
  encompassing_summary_keys_unique = !anyDuplicated(
    enc_tbl[c("site", "family", "normalization_spec", "adjustment")]
  ),
  encompassing_summary_all_alternatives_present = setequal(
    unique(enc_tbl$normalization_spec),
    normalization_adjustment_specs$normalization_spec
  ),
  encompassing_summary_sample_sizes_match_family_samples = all(purrr::map2_lgl(
    enc_tbl$family,
    enc_tbl$n,
    ~ .y == nrow(encompassing_analysis_data[[.x]])
  )),
  encompassing_summary_event_counts_match_family_samples = all(purrr::map2_lgl(
    enc_tbl$family,
    enc_tbl$n_events,
    ~ .y == sum(encompassing_analysis_data[[.x]]$mortality_event_60 == 1L)
  )),
  encompassing_summary_sample_sizes_positive = all(enc_tbl$n > 0L),
  encompassing_summary_event_counts_positive = all(enc_tbl$n_events > 0L),
  encompassing_summary_usable_p_values_valid = check_subset(
    enc_tbl,
    enc_tbl$computationally_usable,
    function(data) {
      all(
        data$pfvc_beyond_pbw_p >= 0 &
          data$pfvc_beyond_pbw_p <= 1 &
          data$pbw_beyond_pfvc_p >= 0 &
          data$pbw_beyond_pfvc_p <= 1
      )
    }
  ),
  encompassing_summary_usable_aic_and_vif_values_finite = check_subset(
    enc_tbl,
    enc_tbl$computationally_usable,
    function(data) {
      all(
        is.finite(data$daic_pfvc_beyond_pbw) &
          is.finite(data$daic_pbw_beyond_pfvc) &
          is.finite(data$max_vif_both)
      )
    }
  ),
  encompassing_summary_unusable_results_suppressed = check_subset(
    enc_tbl,
    !enc_tbl$computationally_usable,
    function(data) {
      all(
        is.na(data$pfvc_beyond_pbw_p) &
          is.na(data$pbw_beyond_pfvc_p) &
          is.na(data$daic_pfvc_beyond_pbw) &
          is.na(data$daic_pbw_beyond_pfvc) &
          is.na(data$max_vif_both)
      )
    }
  )
)

if (!isTRUE(all(unlist(encompassing_summary_checks, use.names = FALSE)))) {
  print(enc_tbl, n = Inf, width = Inf)
  print(encompassing_summary_checks)
}

stopifnot(isTRUE(all(unlist(encompassing_summary_checks, use.names = FALSE))))

write_csv(enc_tbl, encompassing_csv_path)

# =============================================================================
# Encompassing-analysis figure
# =============================================================================

encompassing_plot_data <- enc_tbl %>%
  filter(computationally_usable) %>%
  transmute(
    family = factor(family, levels = exps$family),
    normalization_spec,
    normalization_label = factor(
      normalization_label,
      levels = normalization_adjustment_specs$normalization_label
    ),
    adjustment,
    adjusted = factor(adjusted, levels = unname(adjustment_labels)),
    `PFVC beyond PBW` = daic_pfvc_beyond_pbw,
    `PBW beyond PFVC` = daic_pbw_beyond_pfvc
  ) %>%
  pivot_longer(
    cols = c(`PFVC beyond PBW`, `PBW beyond PFVC`),
    names_to = "direction",
    values_to = "delta_aic"
  )

encompassing_facet_scaffold <- encompassing_fit_grid %>%
  select(family, normalization_spec, adjustment) %>%
  left_join(
    normalization_adjustment_specs %>% select(normalization_spec, normalization_label),
    by = "normalization_spec"
  ) %>%
  mutate(
    family = factor(family, levels = exps$family),
    normalization_label = factor(
      normalization_label,
      levels = normalization_adjustment_specs$normalization_label
    ),
    adjusted = factor(
      unname(adjustment_labels[adjustment]),
      levels = unname(adjustment_labels)
    ),
    delta_aic = 0
  )

encompassing_suppressed_panels <- enc_tbl %>%
  filter(!computationally_usable) %>%
  transmute(
    family = factor(family, levels = exps$family),
    normalization_spec,
    normalization_label = factor(
      normalization_label,
      levels = normalization_adjustment_specs$normalization_label
    ),
    adjustment,
    adjusted = factor(adjusted, levels = unname(adjustment_labels)),
    delta_aic = 0,
    label = str_wrap(
      "Comparisons suppressed: existing computational-use criteria were not met.",
      width = 28
    )
  )
p_enc <- ggplot() +
  geom_blank(
    data = encompassing_facet_scaffold,
    aes(x = delta_aic, y = family)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_col(
    data = encompassing_plot_data,
    aes(x = delta_aic, y = family, fill = direction),
    position = position_dodge(width = 0.7),
    width = 0.6
  ) +
  geom_text(
    data = encompassing_suppressed_panels,
    aes(x = delta_aic, y = family, label = label),
    inherit.aes = FALSE,
    color = "#B2182B",
    size = 3.1,
    lineheight = 1.05,
    hjust = 0
  ) +
  facet_grid(normalization_label ~ adjusted, drop = FALSE) +
  scale_fill_manual(
    values = c(
      "PFVC beyond PBW" = "#009E73",
      "PBW beyond PFVC" = "#D55E00"
    ),
    name = NULL
  ) +
  labs(
    title = "Incremental prognostic information from PFVC and PBW",
    subtitle = str_wrap(
      paste0(
        site_name,
        " — improvement in AIC from adding the second size reference to a mechanic-plus-size model. ",
        "Positive values indicate lower AIC after adding the second term."
      ),
      width = 105
    ),
    x = "AIC improvement after adding the second size term",
    y = NULL,
    caption = str_wrap(
      paste0(
        "Comparisons are nested one-parameter tests on common samples. Nonsignificance does not establish equivalence ",
        "or absence of incremental information. Displayed results may remain exploratory because of sparse events."
      ),
      width = 115
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 8)),
    plot.caption = element_text(hjust = 0, margin = margin(t = 8)),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
  )

ggsave(
  filename = encompassing_pdf_path, plot = p_enc,
  width = 14, height = 24, units = "in", limitsize = TRUE
)

encompassing_output_qc <- tibble(
  output = c("Encompassing-analysis results", "Encompassing-analysis figure"),
  path = c(encompassing_csv_path, encompassing_pdf_path)
) %>%
  mutate(
    exists = file.exists(path),
    size_bytes = purrr::map_dbl(path, function(output_path) {
      if (file.exists(output_path)) file.info(output_path, extra_cols = FALSE)$size[[1]] else NA_real_
    }),
    nonempty = exists & is.finite(size_bytes) & size_bytes > 0
  )

expected_encompassing_keys <- encompassing_fit_grid %>%
  select(family, normalization_spec, adjustment)

observed_encompassing_scaffold_keys <- encompassing_facet_scaffold %>%
  transmute(family = as.character(family), normalization_spec, adjustment)

expected_encompassing_annotation_keys <- enc_tbl %>%
  filter(!computationally_usable) %>%
  select(family, normalization_spec, adjustment)

observed_encompassing_annotation_keys <- encompassing_suppressed_panels %>%
  transmute(family = as.character(family), normalization_spec, adjustment)

encompassing_output_checks <- list(
  encompassing_output_files_exist = all(encompassing_output_qc$exists),
  encompassing_output_files_nonempty = all(encompassing_output_qc$nonempty),
  encompassing_scaffold_keys_unique = !anyDuplicated(observed_encompassing_scaffold_keys),
  encompassing_all_expected_scaffold_rows_present = nrow(anti_join(
    expected_encompassing_keys,
    observed_encompassing_scaffold_keys,
    by = c("family", "normalization_spec", "adjustment")
  )) == 0L,
  encompassing_no_unexpected_scaffold_rows = nrow(anti_join(
    observed_encompassing_scaffold_keys,
    expected_encompassing_keys,
    by = c("family", "normalization_spec", "adjustment")
  )) == 0L,
  encompassing_scaffold_labels_complete = all(complete.cases(
    encompassing_facet_scaffold[c("family", "normalization_label", "adjusted")]
  )),
  encompassing_plotted_values_finite = all(is.finite(encompassing_plot_data$delta_aic)),
  encompassing_plot_labels_complete = all(complete.cases(
    encompassing_plot_data[c("family", "normalization_label", "adjusted")]
  )),
  encompassing_suppression_annotation_keys_unique = !anyDuplicated(observed_encompassing_annotation_keys),
  encompassing_all_suppressed_sets_annotated = nrow(anti_join(
    expected_encompassing_annotation_keys,
    observed_encompassing_annotation_keys,
    by = c("family", "normalization_spec", "adjustment")
  )) == 0L,
  encompassing_no_unexpected_suppression_annotations = nrow(anti_join(
    observed_encompassing_annotation_keys,
    expected_encompassing_annotation_keys,
    by = c("family", "normalization_spec", "adjustment")
  )) == 0L
)

if (!isTRUE(all(unlist(encompassing_output_checks, use.names = FALSE)))) {
  print(encompassing_output_qc, n = Inf, width = Inf)
  print(encompassing_output_checks)
}

stopifnot(isTRUE(all(unlist(encompassing_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$encompassing_analysis <- list(
  interpretation = paste0(
    "Nested one-parameter tests evaluate PFVC beyond PBW and PBW beyond PFVC conditional on the mechanic, selected ",
    "alternative adjustment variable, and specified covariates. ",
    "Nonsignificance does not establish equivalence or absence of information. ",
    "Summary results are suppressed when either direction fails the existing computational-use criteria; ",
    "direction-specific diagnostics remain retained."
  ),
  fit_grid = encompassing_fit_grid,
  normalization_specifications = normalization_adjustment_specs,
  common_sample_configuration = normalization_analysis_diagnostics$encompassing_sample_configuration,
  formulas = purrr::map(encompassing_fits, function(fit) {
    purrr::map_chr(fit$formulas, ~ paste(deparse(.x), collapse = " "))
  }),
  sample_summary = encompassing_sample_qc,
  source_rows = purrr::map(
    encompassing_analysis_data,
    ~ .x %>% select(.source_row, hospitalization_id)
  ),
  model_diagnostics = encompassing_model_qc,
  row_validation = encompassing_row_qc,
  cross_specification_row_validation = encompassing_cross_specification_row_qc,
  comparison_diagnostics = encompassing_comparison_qc,
  summary = enc_tbl,
  plotted_data = encompassing_plot_data,
  panel_scaffold = encompassing_facet_scaffold,
  suppressed_panels = encompassing_suppressed_panels,
  layout = list(
    row_variable = "normalization_label",
    column_variable = "adjusted",
    y_variable = "family",
    width_inches = 14,
    height_inches = 24,
    visual_review_status = "Not established by runtime checks"
  ),
  output = encompassing_output_qc
)

normalization_analysis_diagnostics$warnings$encompassing_analysis <- enc_tbl %>%
  filter(comparison_status != "Eligible for conventional interpretation")

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  encompassing_sample_checks,
  encompassing_fit_checks,
  encompassing_summary_checks,
  encompassing_output_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# Final Script 05 output and diagnostic validation
# =============================================================================

message("Script 05 | validating outputs and synchronizing diagnostics")
normalization_expansion_inventory <- tibble(
  component = c(
    "Prognostic models",
    "Prognostic mechanic-exposure coefficients",
    "Form-versus-physiology summary",
    "Demographic modification comparisons",
    "Demographic modification models",
    "Interaction ladder models",
    "Interaction ladder comparisons",
    "Interaction ladder summary",
    "Balanced single-modifier models",
    "Balanced single-modifier comparisons",
    "Single-modifier curve panels",
    "Single-modifier curve rows",
    "Encompassing models",
    "Encompassing comparisons",
    "Encompassing summary"
  ),
  expected_rows = c(
    nrow(prognostic_fit_grid),
    expected_exposure_coefficient_rows,
    nrow(decomposition_grid),
    nrow(demographic_modification_fit_grid),
    2L * nrow(demographic_modification_fit_grid),
    length(ladder_rhs) * nrow(ladder_fit_grid),
    nrow(ladder_comparison_specs) * nrow(ladder_fit_grid),
    nrow(ladder_fit_grid),
    length(single_modifier_formulas) * nrow(single_modifier_fit_grid),
    2L * nrow(single_modifier_fit_grid),
    nrow(modifier_curve_grid),
    100L * nrow(modifier_curve_grid),
    3L * nrow(encompassing_fit_grid),
    2L * nrow(encompassing_fit_grid),
    nrow(encompassing_fit_grid)
  ),
  observed_rows = c(
    nrow(prognostic_model_qc),
    nrow(prog_coefs),
    nrow(decomp_tbl),
    nrow(invar_tbl),
    nrow(demographic_modification_model_qc),
    nrow(ladder_model_qc),
    nrow(ladder_comparison_qc),
    nrow(ladder_tbl),
    nrow(single_modifier_model_qc),
    nrow(single_modifier_diagnostics),
    nrow(modifier_curve_status),
    nrow(modifier_curve_data),
    nrow(encompassing_model_qc),
    nrow(encompassing_comparison_qc),
    nrow(enc_tbl)
  )
) %>%
  mutate(count_matches = observed_rows == expected_rows)

normalization_expansion_tables <- list(
  prognostic_fit = prog_tbl,
  prognostic_coefficients = prog_coefs,
  decomposition = decomp_tbl,
  demographic_modification = invar_tbl,
  interaction_ladder = ladder_tbl,
  single_modifier = single_modifier_diagnostics,
  modifier_curves = modifier_curve_data,
  encompassing = enc_tbl
)

normalization_expansion_schema_qc <- purrr::imap_dfr(
  normalization_expansion_tables,
  function(data, component_name) {
    required_columns <- c(
      "site", "normalization_spec", "normalization_variable", "normalization_label", "normalization_kind"
    )
    missing_columns <- setdiff(required_columns, names(data))
    metadata_complete <- FALSE
    registry_matches <- FALSE
    
    if (length(missing_columns) == 0L) {
      metadata_complete <- all(complete.cases(data[required_columns]))
      
      observed_registry <- data %>%
        distinct(normalization_spec, normalization_variable, normalization_label, normalization_kind) %>%
        mutate(across(everything(), as.character))
      
      registry_matches <- nrow(anti_join(
        observed_registry,
        normalization_adjustment_specs,
        by = c("normalization_spec", "normalization_variable", "normalization_label", "normalization_kind")
      )) == 0L && setequal(
        observed_registry$normalization_spec,
        normalization_adjustment_specs$normalization_spec
      )
    }
    
    tibble(
      component = component_name,
      required_columns_present = length(missing_columns) == 0L,
      missing_columns = paste(missing_columns, collapse = "; "),
      metadata_complete = metadata_complete,
      registry_matches = registry_matches,
      site_matches = "site" %in% names(data) && all(data$site == site_name)
    )
  }
)

normalization_expansion_checks <- list(
  normalization_expansion_inventory_counts_match = all(normalization_expansion_inventory$count_matches),
  normalization_expansion_export_columns_present = all(
    normalization_expansion_schema_qc$required_columns_present
  ),
  normalization_expansion_export_metadata_complete = all(normalization_expansion_schema_qc$metadata_complete),
  normalization_expansion_export_registry_matches = all(normalization_expansion_schema_qc$registry_matches),
  normalization_expansion_export_sites_match = all(normalization_expansion_schema_qc$site_matches)
)

if (!isTRUE(all(unlist(normalization_expansion_checks, use.names = FALSE)))) {
  print(normalization_expansion_inventory, n = Inf, width = Inf)
  print(normalization_expansion_schema_qc, n = Inf, width = Inf)
  print(normalization_expansion_checks)
}

stopifnot(isTRUE(all(unlist(normalization_expansion_checks, use.names = FALSE))))

normalization_analysis_diagnostics$expansion_validation <- list(
  registry = normalization_adjustment_specs,
  inventory = normalization_expansion_inventory,
  export_schema_validation = normalization_expansion_schema_qc,
  checks = normalization_expansion_checks,
  interpretation = paste0(
    "These checks validate the configured expansion counts and alternative identifiers. ",
    "They do not establish statistical reliability, clinical validity, report compatibility, or pooling compatibility."
  )
)

normalization_analysis_diagnostics$metadata$n_normalization_specifications <-
  nrow(normalization_adjustment_specs)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  normalization_expansion_checks
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

ENCOMPASSING_VIF_WARNING_THRESHOLD <- 10

encompassing_collinearity_warnings <- enc_tbl %>%
  filter(computationally_usable, is.finite(max_vif_both), max_vif_both >= ENCOMPASSING_VIF_WARNING_THRESHOLD) %>%
  transmute(
    family, normalization_spec, normalization_variable, normalization_label, normalization_kind,
    adjustment, adjusted, n, n_events, max_vif_both,
    warning = "Both-size encompassing model has substantial collinearity.",
    interpretation = paste0(
      "Individual PBW and PFVC coefficients require cautious interpretation because the model predictors are related. ",
      "The warning does not change the fitted model or exclude observations. ",
      "Nested tests and AIC changes remain exploratory."
    )
  )

normalization_analysis_diagnostics$warnings$encompassing_collinearity <- encompassing_collinearity_warnings

if (nrow(encompassing_collinearity_warnings) > 0L) {
  message(
    "Script 05 | collinearity caution: ",
    nrow(encompassing_collinearity_warnings),
    " encompassing comparison sets have maximum VIF >= ",
    ENCOMPASSING_VIF_WARNING_THRESHOLD,
    "; individual size coefficients require caution. ",
    "Affected sets and VIF values retained in diagnostics."
  )
}

mechanical_power_status <- prognostic_model_qc %>%
  filter(family == "Mechanical power") %>%
  summarise(
    n_models = n(),
    n_nonconverged = sum(!converged),
    n_boundary = sum(boundary),
    n_with_extreme_fitted_probabilities = sum(n_extreme_fitted_probabilities > 0L),
    n_computationally_usable = sum(computationally_usable),
    minimum_events_per_parameter = min(events_per_parameter)
  )

usable_encompassing_vif <- enc_tbl$max_vif_both[
  enc_tbl$computationally_usable & is.finite(enc_tbl$max_vif_both)
]

maximum_usable_encompassing_vif <- if (length(usable_encompassing_vif) > 0L) {
  max(usable_encompassing_vif)
} else {
  NA_real_
}

final_interpretation <- c(
  "All retained data-integrity, model-integrity, suppression, and output checks passed.",
  paste0(
    dp_qc$n_pip_derived, " of ", dp_qc$n_dp_valid,
    " valid estimated driving-pressure values were PIP-derived fallbacks."
  ),
  paste0(
    "The minimum events per fitted parameter across candidate prognostic models was ",
    round(min(prognostic_model_qc$events_per_parameter), 2), "."
  ),
  paste0(
    "Among ", mechanical_power_status$n_models, " mechanical-power prognostic models, ",
    mechanical_power_status$n_computationally_usable, " met the prespecified computational criteria, ",
    mechanical_power_status$n_nonconverged, " did not converge, and ",
    mechanical_power_status$n_with_extreme_fitted_probabilities,
    " had extreme fitted probabilities."
  ),
  paste0(
    "The maximum VIF among computationally usable both-size encompassing models was ",
    ifelse(is.finite(maximum_usable_encompassing_vif), round(maximum_usable_encompassing_vif, 1), "not available"),
    "."
  ),
  "Events-per-parameter thresholds are interpretation and visualization heuristics, not absolute validity boundaries."
)

normalization_nonreport_outputs <- tibble(
  output = c(
    "Discordance summary",
    "Discordance demographic model",
    "Cohort-relative tercile-change summary",
    "Part 1 discordance figure",
    "Continuous-age discordance summary",
    "Age-pattern heterogeneity summary",
    "Cohort-relative tercile change versus age",
    "Part 1b age-pattern figure",
    "Prognostic fit comparison",
    "Prognostic coefficient estimates",
    "Form-versus-physiology decomposition",
    "Demographic effect-modification tests",
    "Prognostic fit figure",
    "Mechanic interaction ladder",
    "Single-modifier diagnostics",
    "Single-modifier curve data",
    "Single-modifier curve figure",
    "Encompassing-analysis results",
    "Encompassing-analysis figure"
  ),
  path = c(
    discordance_summary_path,
    discordance_demographics_path,
    reclassification_output_path,
    part1_figure_path,
    discordance_age_output_path,
    discordance_heterogeneity_output_path,
    reclassification_age_output_path,
    part1b_figure_path,
    prognostic_fit_output_path,
    prognostic_coefficient_output_path,
    form_physiology_output_path,
    demographic_modification_output_path,
    prognostic_figure_path,
    ladder_output_path,
    single_modifier_output_path,
    modifier_curve_csv_path,
    modifier_curve_pdf_path,
    encompassing_csv_path,
    encompassing_pdf_path
  )
) %>%
  mutate(
    exists = file.exists(path),
    size_bytes = purrr::map_dbl(path, function(output_path) {
      if (file.exists(output_path)) file.info(output_path, extra_cols = FALSE)$size[[1]] else NA_real_
    }),
    nonempty = exists & is.finite(size_bytes) & size_bytes > 0,
    extension = tools::file_ext(path)
  )

nonreport_output_checks <- list(
  expected_nonreport_output_count = nrow(normalization_nonreport_outputs) == 19L,
  all_nonreport_outputs_exist = all(normalization_nonreport_outputs$exists),
  all_nonreport_outputs_nonempty = all(normalization_nonreport_outputs$nonempty),
  nonreport_output_paths_unique = !anyDuplicated(normalization_nonreport_outputs$path)
)

if (!isTRUE(all(unlist(nonreport_output_checks, use.names = FALSE)))) {
  print(normalization_nonreport_outputs, n = Inf, width = Inf)
  print(nonreport_output_checks)
}

stopifnot(isTRUE(all(unlist(nonreport_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$nonreport_output_manifest <- normalization_nonreport_outputs
normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  nonreport_output_checks
)

normalization_analysis_diagnostics$metadata$report_input_stage <-
  "Substantive diagnostics complete; final report and manifest validation pending"
normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# Final Script 05 analysis-output validation
# =============================================================================
# Report rendering is separate. Validate analysis artifacts here without treating
# an older diagnostic PDF as proof of reporting completion for this run.

normalization_expected_outputs <- normalization_nonreport_outputs %>%
  mutate(
    exists = file.exists(path) & !dir.exists(path),
    size_bytes = purrr::map_dbl(path, function(output_path) {
      if (file.exists(output_path) && !dir.exists(output_path)) {
        file.info(output_path, extra_cols = FALSE)$size[[1]]
      } else {
        NA_real_
      }
    }),
    nonempty = exists & is.finite(size_bytes) & size_bytes > 0,
    extension = tools::file_ext(path)
  )

duplicate_output_paths <- normalization_expected_outputs %>%
  count(path, name = "n_uses") %>%
  filter(n_uses > 1L)

# Count files, not PDF pages: Part 1 is one multipage PDF artifact.
expected_output_type_counts <- tibble(
  extension = c("csv", "pdf"),
  expected_n = c(14L, 5L)
)

observed_output_type_counts <- normalization_expected_outputs %>%
  count(extension, name = "observed_n") %>%
  full_join(expected_output_type_counts, by = "extension") %>%
  mutate(
    observed_n = replace_na(observed_n, 0L),
    expected_n = replace_na(expected_n, 0L),
    count_matches = observed_n == expected_n
  )

final_output_checks <- list(
  expected_output_count = nrow(normalization_expected_outputs) == 19L,
  all_expected_outputs_exist = all(normalization_expected_outputs$exists),
  all_expected_outputs_nonempty = all(normalization_expected_outputs$nonempty),
  output_paths_unique = nrow(duplicate_output_paths) == 0L,
  expected_output_type_counts_match = all(observed_output_type_counts$count_matches)
)

if (!isTRUE(all(unlist(final_output_checks, use.names = FALSE)))) {
  print(normalization_expected_outputs, n = Inf, width = Inf)
  print(observed_output_type_counts, n = Inf)
  print(duplicate_output_paths, n = Inf, width = Inf)
  print(final_output_checks)
}

stopifnot(isTRUE(all(unlist(final_output_checks, use.names = FALSE))))

normalization_analysis_diagnostics$output_manifest <- normalization_expected_outputs
normalization_analysis_diagnostics$output_type_counts <- observed_output_type_counts

normalization_analysis_diagnostics$reporting <- list(
  status = "PENDING",
  source = file.path(config$utils_dir, "normalization_diagnostic_report.qmd"),
  intended_output = file.path(
    report_dir,
    paste0("normalization_diagnostic_report_", site_name, ".pdf")
  ),
  report_rendered_by_script05 = FALSE,
  report_inclusive_expected_file_count = 20L,
  report_inclusive_expected_type_counts = tibble(
    extension = c("csv", "pdf"),
    expected_n = c(14L, 6L)
  ),
  note = paste0(
    "Analysis artifacts are validated separately from reporting. The reporting workflow must render from the saved ",
    "diagnostics and validate the report-inclusive manifest. An existing PDF does not establish current-run completion."
  )
)

normalization_analysis_diagnostics$final_summary <- list(
  completed_at = Sys.time(),
  scope = "Script 05 analysis artifacts and diagnostics; reporting pending",
  n_cross_sectional = nrow(cross_sectional),
  n_deaths_cross_sectional = sum(cross_sectional$mortality_event_60 == 1L),
  n_expected_outputs = nrow(normalization_expected_outputs),
  n_warnings = sum(purrr::map_int(
    normalization_analysis_diagnostics$warnings,
    function(warning_object) {
      if (is.data.frame(warning_object)) nrow(warning_object) else length(warning_object)
    }
  )),
  overall_status = "PASS WITH INTERPRETATION WARNINGS",
  reporting_status = "PENDING",
  interpretation = final_interpretation
)

normalization_analysis_diagnostics$checks <- c(
  normalization_analysis_diagnostics$checks,
  final_output_checks
)

# Computational and interpretation status describe analysis, not report completion.
normalization_analysis_diagnostics$metadata$computational_status <- "PASS"
normalization_analysis_diagnostics$metadata$interpretation_status <- "WARNINGS PRESENT"
normalization_analysis_diagnostics$metadata$overall_status <- "PASS WITH INTERPRETATION WARNINGS"
normalization_analysis_diagnostics$metadata$reporting_status <- "PENDING"
normalization_analysis_diagnostics$metadata$report_input_stage <-
  "Analysis artifacts and substantive diagnostics validated; separate report rendering pending"
normalization_analysis_diagnostics$metadata$completed_at <- Sys.time()
normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()

save_normalization_diagnostics(normalization_analysis_diagnostics)

# =============================================================================
# Final diagnostic read-back
# =============================================================================

normalization_diagnostics_readback <- readRDS(normalization_analysis_diagnostics_path)

final_diagnostic_checks <- list(
  complete_diagnostic_readback_matches = identical(
    normalization_diagnostics_readback,
    normalization_analysis_diagnostics
  ),
  diagnostic_site_matches = identical(normalization_diagnostics_readback$metadata$site, site_name),
  diagnostic_computational_status_matches = identical(
    normalization_diagnostics_readback$metadata$computational_status,
    "PASS"
  ),
  diagnostic_status_matches = identical(
    normalization_diagnostics_readback$metadata$overall_status,
    "PASS WITH INTERPRETATION WARNINGS"
  ),
  diagnostic_reporting_status_pending = identical(
    normalization_diagnostics_readback$metadata$reporting_status,
    "PENDING"
  ),
  diagnostic_output_manifest_matches = identical(
    normalization_diagnostics_readback$output_manifest,
    normalization_expected_outputs
  ),
  all_saved_runtime_checks_pass = isTRUE(all(unlist(
    normalization_diagnostics_readback$checks,
    use.names = FALSE
  ))),
  final_summary_present = is.list(normalization_diagnostics_readback$final_summary)
)

final_diagnostic_check_summary <- tibble(
  check = names(final_diagnostic_checks),
  passed = vapply(final_diagnostic_checks, isTRUE, logical(1))
)

if (!all(final_diagnostic_check_summary$passed)) {
  print(final_diagnostic_check_summary, n = Inf)
}

stopifnot(all(final_diagnostic_check_summary$passed))

normalization_analysis_diagnostics$final_readback_validation <- list(
  checks = final_diagnostic_check_summary,
  note = "Checks validate the preceding analysis artifact; this summary is then retained and fully read back."
)

normalization_analysis_diagnostics$metadata$updated_at <- Sys.time()
save_normalization_diagnostics(normalization_analysis_diagnostics)

terminal_diagnostics_readback <- readRDS(normalization_analysis_diagnostics_path)

stopifnot(
  identical(terminal_diagnostics_readback, normalization_analysis_diagnostics),
  isTRUE(all(terminal_diagnostics_readback$final_readback_validation$checks$passed)),
  isTRUE(all(unlist(terminal_diagnostics_readback$checks, use.names = FALSE)))
)

# =============================================================================
# Synchronize Script 05 analysis diagnostics into the latest Script 04 model QC
# =============================================================================

script05_model_qc_section <- "normalization_analysis_script05"

if (!file.exists(script04_model_qc_path) || dir.exists(script04_model_qc_path)) {
  stop("The model QC RDS was not found for Script 05 synchronization: ", script04_model_qc_path)
}

latest_model_qc <- readRDS(script04_model_qc_path)

stopifnot(
  is.list(latest_model_qc),
  is.list(latest_model_qc$metadata),
  identical(latest_model_qc$metadata$site, site_name),
  is.list(latest_model_qc$checks),
  identical(terminal_diagnostics_readback$metadata$site, site_name),
  identical(terminal_diagnostics_readback$metadata$computational_status, "PASS")
)

existing_script05_section <- latest_model_qc[[script05_model_qc_section]]

if (!is.null(existing_script05_section)) {
  if (
    !is.list(existing_script05_section) ||
    !is.list(existing_script05_section$metadata) ||
    !identical(existing_script05_section$metadata$site, site_name) ||
    !identical(
      existing_script05_section$metadata$script,
      terminal_diagnostics_readback$metadata$script
    )
  ) {
    stop(
      "The existing model QC section `", script05_model_qc_section,
      "` does not match the expected Script 05 identity. It was not overwritten."
    )
  }
}

preserved_model_qc_sections <- setdiff(names(latest_model_qc), script05_model_qc_section)
model_qc_before_script05_sync <- latest_model_qc

latest_model_qc[[script05_model_qc_section]] <- terminal_diagnostics_readback

stopifnot(identical(
  latest_model_qc[preserved_model_qc_sections],
  model_qc_before_script05_sync[preserved_model_qc_sections]
))

saveRDS(latest_model_qc, script04_model_qc_path, compress = "xz")
model_qc_after_script05_sync <- readRDS(script04_model_qc_path)

script05_model_qc_sync_checks <- list(
  full_object_readback_matches = identical(model_qc_after_script05_sync, latest_model_qc),
  script05_section_matches_standalone_artifact = identical(
    model_qc_after_script05_sync[[script05_model_qc_section]],
    terminal_diagnostics_readback
  ),
  existing_sections_preserved = identical(
    model_qc_after_script05_sync[preserved_model_qc_sections],
    model_qc_before_script05_sync[preserved_model_qc_sections]
  ),
  site_matches = identical(model_qc_after_script05_sync$metadata$site, site_name)
)

script05_model_qc_sync_summary <- tibble(
  check = names(script05_model_qc_sync_checks),
  passed = vapply(script05_model_qc_sync_checks, isTRUE, logical(1))
)

if (!all(script05_model_qc_sync_summary$passed)) {
  print(script05_model_qc_sync_summary, n = Inf)
}

stopifnot(all(script05_model_qc_sync_summary$passed))


message("Script 05 analysis complete; reporting pending.")
