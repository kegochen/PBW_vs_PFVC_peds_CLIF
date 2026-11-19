# =============================================================================
# Script 04: Analysis
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================

rm(list = ls())

library(tidyverse)
library(arrow)
library(here)
library(gtsummary)
library(gt)
library(survival)
library(survminer)
library(splines)
library(patchwork)
library(parameters)
library(broom)
library(algorithmDiagnostics)
library(EValue)

save_script04_plot <- function(...) {
  withCallingHandlers(
    ggplot2::ggsave(...),
    message = function(condition) {
      text <- trimws(gsub(
        "[[:space:]]+", " ", conditionMessage(condition)
      ))

      routine_message <- text %in% c(
        "`geom_smooth()` using formula = 'y ~ x'",
        paste0(
          "Orientation inferred to be along y-axis; override with ",
          "`position_quasirandom(orientation = 'x')`"
        )
      )

      if (routine_message) {
        restart <- findRestart("muffleMessage")
        if (!is.null(restart)) invokeRestart(restart)
      }
    }
  )
}

# =============================================================================
# Load validated repository configuration
# =============================================================================

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

required_script04_input_directories <- c(output_dir, config$cohort_diagnostic_dir)
missing_script04_input_directories <- required_script04_input_directories[
  !dir.exists(required_script04_input_directories)
]

if (length(missing_script04_input_directories) > 0L) {
  stop(
    "Script 04 prerequisite directory(s) not found: ",
    paste(missing_script04_input_directories, collapse = "; "),
    ". Complete Scripts 01-03 before executing Script 04."
  )
}

consort_utility_path <- file.path(config$utils_dir, "consort_diagram.R")

if (!file.exists(consort_utility_path)) {
  stop("Repository CONSORT utility not found: ", consort_utility_path)
}

source(consort_utility_path)

# Establish explicit destinations for existing deliverables.
# Detailed diagnostic RDS files use their separate local diagnostic directories.
table_dir <- config$table_dir
figure_dir <- config$figure_dir
report_dir <- config$report_dir

script04_deliverable_directories <- c(table_dir, figure_dir, report_dir)

for (directory in script04_deliverable_directories) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

missing_script04_deliverable_directories <- script04_deliverable_directories[
  !dir.exists(script04_deliverable_directories)
]

if (length(missing_script04_deliverable_directories) > 0L) {
  stop(
    "Script 04 could not create deliverable directory(s): ",
    paste(missing_script04_deliverable_directories, collapse = "; ")
  )
}

rm(directory)

message("Script 04 starting | site: ", site_name)
message("Intermediate input directory: ", output_dir)
message("Final output directory: ", final_dir)

# =============================================================================
# Table 1 QC configuration
# =============================================================================
# Detailed diagnostic RDS files remain local under diagnostics/table1.

table1_diagnostic_dir <- config$table1_diagnostic_dir

dir.create(table1_diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(table1_diagnostic_dir)) {
  stop("Script 04 could not create the Table 1 diagnostic directory: ", table1_diagnostic_dir)
}

table1_qc_path <- file.path(
  table1_diagnostic_dir,
  paste0("table1_qc_", site_name, ".rds")
)

save_table1_qc <- function(qc, path) {
  if (!is.list(qc) || !is.list(qc$checks) || length(qc$checks) == 0L) {
    stop("Table 1 QC must contain a nonempty list of integrity checks.")
  }
  
  if (any(lengths(qc$checks) != 1L)) {
    stop("Every Table 1 integrity check must return exactly one value.")
  }
  
  check_values <- vapply(qc$checks, isTRUE, logical(1))
  
  if (!all(check_values)) {
    stop(
      "Table 1 QC cannot be saved as passed. Failed or indeterminate check(s): ",
      paste(names(qc$checks)[!check_values], collapse = ", ")
    )
  }
  
  saveRDS(qc, path, compress = "xz")
  
  if (!file.exists(path)) {
    stop("Table 1 QC RDS was not created: ", path)
  }
  
  qc_readback <- readRDS(path)
  
  if (!identical(qc_readback, qc)) {
    stop("Table 1 QC RDS failed its complete read-back comparison: ", path)
  }
  
  invisible(path)
}

# =============================================================================
# Load data
# =============================================================================

required_script04_input_files <- c(
  "analysis_cross_sectional.parquet",
  "analysis_all_timepoints.parquet",
  "analysis_broad_pfvc.parquet"
)

script04_input_paths <- file.path(output_dir, required_script04_input_files)
missing_script04_input_paths <- script04_input_paths[
  !file.exists(script04_input_paths)
]

if (length(missing_script04_input_paths) > 0L) {
  stop(
    "Script 04 required input file(s) not found: ",
    paste(missing_script04_input_paths, collapse = "; "),
    ". Complete Script 03 before executing Script 04."
  )
}

directory_script04_input_paths <- script04_input_paths[
  dir.exists(script04_input_paths)
]

if (length(directory_script04_input_paths) > 0L) {
  stop(
    "Script 04 input path(s) identify directories rather than files: ",
    paste(directory_script04_input_paths, collapse = "; ")
  )
}

cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))
analysis_all <- read_parquet(file.path(output_dir, "analysis_all_timepoints.parquet"))

message("Cross-sectional: ", nrow(cross_sectional), " observations")
message("All timepoints: ", nrow(analysis_all), " observations")

# =============================================================================
# Validate the Script 03 input contract
# =============================================================================

cohort_diagnostic_path <- file.path(config$cohort_diagnostic_dir,"cohort_diagnostics.rds")

if (!file.exists(cohort_diagnostic_path)) {
  stop("The combined Scripts 01-03 diagnostic RDS was not found: ", cohort_diagnostic_path)
}

cohort_diagnostics <- readRDS(cohort_diagnostic_path)

if (!is.list(cohort_diagnostics)) {
  stop("The combined Scripts 01-03 diagnostic RDS does not contain a list.")
}

if (!"quality_checks_script03" %in% names(cohort_diagnostics)) {
  stop("The combined diagnostic RDS does not contain quality_checks_script03.")
}

script03_qc <- cohort_diagnostics$quality_checks_script03

required_script03_sections <- c(
  "metadata", "cohort_funnel", "post_intubation_timepoints",
  "oxygenation_timepoint_matching", "psofa_timepoint_matching"
)

missing_script03_sections <- setdiff(required_script03_sections, names(script03_qc))

if (length(missing_script03_sections) > 0L) {
  stop(
    "The Script 03 diagnostic section is missing required component(s): ",
    paste(missing_script03_sections, collapse = ", ")
  )
}

required_script03_timepoint_sections <- c(
  "input_checks", "context_checks", "context_preservation_checks", "selected_tidal_volume_cleaning"
)

missing_script03_timepoint_sections <- setdiff(
  required_script03_timepoint_sections,
  names(script03_qc$post_intubation_timepoints)
)

if (length(missing_script03_timepoint_sections) > 0L) {
  stop(
    "The Script 03 timepoint diagnostics are missing: ",
    paste(missing_script03_timepoint_sections, collapse = ", ")
  )
}

script03_selected_tv_qc <- script03_qc$post_intubation_timepoints$selected_tidal_volume_cleaning
required_selected_tv_sections <- c("policy", "thresholds", "summary", "hospitalization_level", "retention", "checks")
missing_selected_tv_sections <- setdiff(required_selected_tv_sections, names(script03_selected_tv_qc))

if (length(missing_selected_tv_sections) > 0L) {
  stop("The Script 03 selected tidal-volume diagnostics are missing: ", paste(missing_selected_tv_sections, collapse = ", "))
}

allowed_intubation_identification_methods <- c(
  "Respiratory transition",
  "Initially invasive; first documented invasive support in ICU"
)

script04_input_required_cross_sectional_columns <- c(
  "hospitalization_id", "intubation_identification_method", "has_31500_in_hospitalization", "psofa_match_source",
  "tidal_volume_value", "tidal_volume_value_before_cleaning", "tidal_volume_cleaning_reason",
  "sf_ratio_match_offset_hours", "pf_ratio_match_offset_hours", "oi_match_offset_hours",
  "osi_match_offset_hours", "pbw", "pbw_nchs", "abw", "pfvc", "vtpbw", "vtnchs", "vtabw",
  "vtpfvc", "pbwpfvc", "nchspfvc", "abwpfvc"
)

script04_input_required_all_timepoint_columns <- c(
  "hospitalization_id", "intubation_identification_method", "has_31500_in_hospitalization", "psofa_match_source",
  "tidal_volume_value", "tidal_volume_value_before_cleaning", "tidal_volume_cleaning_reason",
  "sf_ratio_match_offset_hours", "pf_ratio_match_offset_hours", "oi_match_offset_hours",
  "osi_match_offset_hours", "pbw", "pbw_nchs", "abw", "pfvc", "vtpbw", "vtnchs", "vtabw",
  "vtpfvc", "pbwpfvc", "nchspfvc", "abwpfvc"
)

missing_cross_sectional_contract_columns <- setdiff(
  script04_input_required_cross_sectional_columns,
  names(cross_sectional)
)

missing_all_timepoint_contract_columns <- setdiff(
  script04_input_required_all_timepoint_columns,
  names(analysis_all)
)

script04_input_contract_summary <- tibble(
  source = c("Cross-sectional analysis", "All retained timepoints"),
  n_rows = c(nrow(cross_sectional), nrow(analysis_all)),
  n_hospitalizations = c(
    n_distinct(cross_sectional$hospitalization_id),
    n_distinct(analysis_all$hospitalization_id)
  )
)

script04_input_contract_checks <- list(
  script03_qc_status_passed = identical(script03_qc$metadata$overall_status, "PASS"),
  cross_sectional_contract_columns_present = length(missing_cross_sectional_contract_columns) == 0L,
  all_timepoint_contract_columns_present = length(missing_all_timepoint_contract_columns) == 0L,
  cross_sectional_hospitalization_ids_complete = !any(is.na(cross_sectional$hospitalization_id)),
  all_timepoint_hospitalization_ids_complete = !any(is.na(analysis_all$hospitalization_id)),
  cross_sectional_hospitalization_ids_unique =
    anyDuplicated(cross_sectional$hospitalization_id) == 0L,
  cross_sectional_count_matches_script03_funnel =
    nrow(cross_sectional) == last(script03_qc$cohort_funnel$n_hospitalizations),
  all_timepoint_hospitalizations_match_cross_sectional = setequal(
    unique(analysis_all$hospitalization_id),
    cross_sectional$hospitalization_id
  ),
  oxygenation_component_window_is_four_hours = isTRUE(
    script03_qc$oxygenation_timepoint_matching$component_match_window_hours == 4),
  oxygenation_timepoint_window_is_five_hours = isTRUE(
    script03_qc$oxygenation_timepoint_matching$timepoint_match_window_hours == 5),
  psofa_maximum_day_is_28 = isTRUE(
    script03_qc$psofa_timepoint_matching$maximum_day == 28),
  psofa_preadmission_grace_is_five_hours = isTRUE(
    script03_qc$psofa_timepoint_matching$preadmission_grace_hours == 5),
  oxygenation_matching_checks_passed = isTRUE(all(
    script03_qc$oxygenation_timepoint_matching$checks$passed
  )),
  psofa_matching_checks_passed = isTRUE(all(
    script03_qc$psofa_timepoint_matching$checks$passed
  )),
  intubation_input_checks_passed = isTRUE(all(
    script03_qc$post_intubation_timepoints$input_checks$passed
  )),
  intubation_context_checks_passed = isTRUE(all(
    script03_qc$post_intubation_timepoints$context_checks$passed
  )),
  cleaned_context_preservation_checks_passed = isTRUE(all(
    script03_qc$post_intubation_timepoints$context_preservation_checks$passed
  )),
  selected_tidal_volume_cleaning_checks_passed = isTRUE(all(
    script03_selected_tv_qc$checks$passed
  )),
  selected_tidal_volume_cleaning_policy_matches = identical(
    script03_selected_tv_qc$policy,
    "Clean the originally selected value; preserve rows and source; no source substitution"
  ),
  cross_sectional_identification_methods_valid = all(
    cross_sectional$intubation_identification_method %in% allowed_intubation_identification_methods
  ),
  all_timepoint_identification_methods_valid = all(
    analysis_all$intubation_identification_method %in% allowed_intubation_identification_methods
  ),
  cross_sectional_selected_volumes_retained = all(
    cross_sectional$tidal_volume_cleaning_reason == "Retained" &
      is.finite(cross_sectional$tidal_volume_value) & cross_sectional$tidal_volume_value > 0
  ),
  all_timepoint_cleaning_status_matches_missingness = identical(
    analysis_all$tidal_volume_cleaning_reason == "Retained",
    !is.na(analysis_all$tidal_volume_value)
  ),
  cross_sectional_procedure_evidence_complete =
    is.logical(cross_sectional$has_31500_in_hospitalization) &&
    !anyNA(cross_sectional$has_31500_in_hospitalization),
  all_timepoint_procedure_evidence_complete =
    is.logical(analysis_all$has_31500_in_hospitalization) &&
    !anyNA(analysis_all$has_31500_in_hospitalization)
)

non_scalar_input_contract_checks <- names(script04_input_contract_checks)[
  lengths(script04_input_contract_checks) != 1L
]

if (length(non_scalar_input_contract_checks) > 0L) {
  stop(
    "The following Script 04 input-contract checks did not return one value: ",
    paste(non_scalar_input_contract_checks, collapse = ", ")
  )
}

script04_input_contract_checks_table <- tibble(
  check = names(script04_input_contract_checks),
  passed = as.logical(unlist(script04_input_contract_checks, use.names = FALSE))
)

message("---- Script 04 upstream input contract ----")

if (length(missing_cross_sectional_contract_columns) > 0L) {
  message(
    "Missing cross-sectional contract columns: ",
    paste(missing_cross_sectional_contract_columns, collapse = ", ")
  )
}

if (length(missing_all_timepoint_contract_columns) > 0L) {
  message(
    "Missing all-timepoint contract columns: ",
    paste(missing_all_timepoint_contract_columns, collapse = ", ")
  )
}


stopifnot(isTRUE(all(script04_input_contract_checks_table$passed)))

# =============================================================================
# Table 1 cohort construction and QC
# =============================================================================

table1_required_vars <- c(
  "age_at_admission", "sex_category", "race_category", "height_cm",
  "pbw", "pbw_nchs", "abw", "pfvc", "sofa_total", "osi", "sf_ratio",
  "mortality_event_60", "vtpbw", "vtnchs", "vtabw", "vtpfvc",
  "pbwpfvc", "nchspfvc", "abwpfvc", "crs", "ers", "vfd_28"
)

missing_table1_vars <- setdiff(table1_required_vars, names(cross_sectional))

stopifnot(length(missing_table1_vars) == 0L)

# Check availability of PBW/PFVC, which defines the Table 1 strata.
pbwpfvc_qc <- cross_sectional %>%
  summarise(
    n_total = n(), n_missing = sum(is.na(pbwpfvc)), n_infinite = sum(is.infinite(pbwpfvc), na.rm = TRUE),
    n_finite = sum(is.finite(pbwpfvc)), pct_missing = 100 * mean(is.na(pbwpfvc))
  )


# Restrict Table 1 to patients who can be assigned to a PBW/PFVC tercile.
n_cross_sectional <- nrow(cross_sectional)

table1_source <- cross_sectional %>%
  filter(is.finite(pbwpfvc))

n_table1 <- nrow(table1_source)
n_excluded_table1 <- n_cross_sectional - n_table1

message(
  "Table 1 cohort: ", format(n_table1, big.mark = ","), " of ", format(n_cross_sectional, big.mark = ","),
  " patients included; ", format(n_excluded_table1, big.mark = ","),
  " excluded because PBW/PFVC was missing or nonfinite."
)

stopifnot(n_table1 > 0L)

# Create approximately equal-sized PBW/PFVC terciles.
table1_source <- table1_source %>%
  mutate(
    pbwpfvc_tercile_num = ntile(pbwpfvc, 3L),
    pbwpfvc_tercile = factor(pbwpfvc_tercile_num, levels = 1:3,
                             labels = c("T1 (Low)", "T2 (Mid)", "T3 (High)"))
  )

tercile_counts <- table1_source %>%
  count(pbwpfvc_tercile, name = "n") %>%
  mutate(percent = 100 * n / sum(n))

tercile_ranges <- table1_source %>%
  group_by(pbwpfvc_tercile) %>%
  summarise(
    n = n(), minimum = min(pbwpfvc), p25 = quantile(pbwpfvc, 0.25), median = median(pbwpfvc),
    p75 = quantile(pbwpfvc, 0.75), maximum = max(pbwpfvc), .groups = "drop"
  )

tercile_ties <- table1_source %>%
  distinct(pbwpfvc, pbwpfvc_tercile) %>%
  count(pbwpfvc, name = "n_terciles") %>%
  filter(n_terciles > 1L) %>%
  arrange(pbwpfvc)


stopifnot(
  sum(tercile_counts$n) == nrow(table1_source),
  !any(is.na(table1_source$pbwpfvc_tercile))
)

# Validate sex and race categories before reapplying their factor levels.
sex_counts <- table1_source %>%
  count(sex_category, sort = TRUE, .drop = FALSE)

race_counts <- table1_source %>%
  count(race_category, sort = TRUE, .drop = FALSE)

expected_sex <- c("Male", "Female")
expected_race <- c("WHITE", "BLACK", "OTHER")

unexpected_sex <- table1_source %>%
  filter(!is.na(sex_category), !sex_category %in% expected_sex) %>%
  distinct(sex_category) %>%
  pull(sex_category)

unexpected_race <- table1_source %>%
  filter(!is.na(race_category), !race_category %in% expected_race) %>%
  distinct(race_category) %>%
  pull(race_category)


stopifnot(length(unexpected_sex) == 0L, length(unexpected_race) == 0L)

table1_source <- table1_source %>%
  mutate(
    sex_category = factor(sex_category, levels = expected_sex),
    race_category = factor(race_category, levels = expected_race)
  )

# Validate 60-day mortality coding before assigning presentation labels.
mortality_counts <- table1_source %>%
  count(mortality_event_60, sort = TRUE, .drop = FALSE)

mortality_values <- sort(unique(na.omit(as.character(table1_source$mortality_event_60))))
mortality_coding_valid <- all(mortality_values %in% c("0", "1", "FALSE", "TRUE"))

mortality_validation <- table1_source %>%
  summarise(
    n_total = n(), n_alive_60 = sum(mortality_event_60 == 0, na.rm = TRUE),
    n_deceased_60 = sum(mortality_event_60 == 1, na.rm = TRUE), n_missing = sum(is.na(mortality_event_60)),
    mortality_percent = 100 * mean(mortality_event_60 == 1, na.rm = TRUE), coding_valid = mortality_coding_valid
  )


stopifnot(mortality_coding_valid)

# Check continuous variables for storage type, missingness, and infinite values.
table1_continuous_vars <- c(
  "age_at_admission", "height_cm", "pbw", "pbw_nchs", "abw", "pfvc", "sofa_total",
  "osi", "sf_ratio", "vtpbw", "vtnchs", "vtabw", "vtpfvc",
  "pbwpfvc", "nchspfvc", "abwpfvc", "crs", "ers", "vfd_28"
)

continuous_qc <- purrr::map_dfr(table1_continuous_vars, function(v) {
  x <- table1_source[[v]]
  
  tibble(
    variable = v, class = paste(class(x), collapse = "/"), n_total = length(x),
    n_nonmissing = sum(!is.na(x)), n_missing = sum(is.na(x)), pct_missing = 100 * mean(is.na(x)),
    n_infinite = if (is.numeric(x)) sum(is.infinite(x), na.rm = TRUE) else NA_integer_
  )
})


continuous_variables_numeric <- all(purrr::map_lgl(table1_source[table1_continuous_vars], is.numeric))
continuous_variables_finite <- !any(continuous_qc$n_infinite > 0, na.rm = TRUE)

stopifnot(continuous_variables_numeric, continuous_variables_finite)

# Summarize continuous-variable ranges and evaluate clear mathematical constraints.
continuous_range_qc <- purrr::map_dfr(table1_continuous_vars, function(v) {
  x <- table1_source[[v]]
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(tibble(
      variable = v, n_finite = 0L, minimum = NA_real_, p01 = NA_real_, p25 = NA_real_,
      median = NA_real_, p75 = NA_real_, p99 = NA_real_, maximum = NA_real_
    ))
  }
  
  tibble(
    variable = v, n_finite = length(x), minimum = min(x), p01 = quantile(x, 0.01, names = FALSE),
    p25 = quantile(x, 0.25, names = FALSE), median = median(x), p75 = quantile(x, 0.75, names = FALSE),
    p99 = quantile(x, 0.99, names = FALSE), maximum = max(x)
  )
})

mathematical_constraint_checks <- list(
  age_nonnegative = all(table1_source$age_at_admission >= 0, na.rm = TRUE),
  height_positive = all(table1_source$height_cm > 0, na.rm = TRUE),
  pbw_positive = all(table1_source$pbw > 0, na.rm = TRUE),
  pbw_nchs_positive = all(table1_source$pbw_nchs > 0, na.rm = TRUE),
  abw_positive = all(table1_source$abw > 0, na.rm = TRUE),
  pfvc_positive = all(table1_source$pfvc > 0, na.rm = TRUE),
  vtpbw_positive = all(table1_source$vtpbw > 0, na.rm = TRUE),
  vtnchs_positive = all(table1_source$vtnchs > 0, na.rm = TRUE),
  vtabw_positive = all(table1_source$vtabw > 0, na.rm = TRUE),
  vtpfvc_positive = all(table1_source$vtpfvc > 0, na.rm = TRUE),
  pbwpfvc_positive = all(table1_source$pbwpfvc > 0, na.rm = TRUE),
  nchspfvc_positive = all(table1_source$nchspfvc > 0, na.rm = TRUE),
  abwpfvc_positive = all(table1_source$abwpfvc > 0, na.rm = TRUE),
  vfd_28_within_0_to_28 = all(table1_source$vfd_28 >= 0 & table1_source$vfd_28 <= 28, na.rm = TRUE)
)


stopifnot(isTRUE(all(unlist(mathematical_constraint_checks))))

# Verify the internal consistency of mathematically related variables.
check_derived_variable <- function(observed, expected, derivation, relative_tolerance = 1e-6) {
  compared <- is.finite(observed) & is.finite(expected)
  n_compared <- sum(compared)
  
  if (n_compared == 0L) {
    return(tibble(
      derivation = derivation, n_compared = 0L, max_absolute_difference = NA_real_,
      max_relative_difference = NA_real_, relative_tolerance = relative_tolerance, passed = FALSE
    ))
  }
  
  absolute_difference <- abs(observed[compared] - expected[compared])
  relative_difference <- absolute_difference / pmax(abs(expected[compared]), .Machine$double.eps)
  
  tibble(
    derivation = derivation, n_compared = n_compared, max_absolute_difference = max(absolute_difference),
    max_relative_difference = max(relative_difference), relative_tolerance = relative_tolerance,
    passed = all(relative_difference <= relative_tolerance)
  )
}

derivation_qc <- bind_rows(
  check_derived_variable(
    observed = table1_source$pbwpfvc,
    expected = table1_source$pbw / table1_source$pfvc,
    derivation = "pbwpfvc = pbw / pfvc"
  ),
  check_derived_variable(
    observed = table1_source$nchspfvc,
    expected = table1_source$pbw_nchs / table1_source$pfvc,
    derivation = "nchspfvc = pbw_nchs / pfvc"
  ),
  check_derived_variable(
    observed = table1_source$abwpfvc,
    expected = table1_source$abw / table1_source$pfvc,
    derivation = "abwpfvc = abw / pfvc"
  ),
  check_derived_variable(
    observed = table1_source$vtpfvc,
    expected = table1_source$vtpbw * table1_source$pbw / (10 * table1_source$pfvc),
    derivation = "vtpfvc = vtpbw * pbw / (10 * pfvc)"
  ),
  check_derived_variable(
    observed = table1_source$vtnchs,
    expected = table1_source$vtpbw * table1_source$pbw / table1_source$pbw_nchs,
    derivation = "vtnchs = vtpbw * pbw / pbw_nchs"
  ),
  check_derived_variable(
    observed = table1_source$vtabw,
    expected = table1_source$vtpbw * table1_source$pbw / table1_source$abw,
    derivation = "vtabw = vtpbw * pbw / abw"
  ),
  check_derived_variable(
    observed = table1_source$ers,
    expected = 1000 / table1_source$crs,
    derivation = "ers = 1000 / crs"
  )
)


derived_variables_consistent <- nrow(derivation_qc) == 7L && all(derivation_qc$passed)

stopifnot(derived_variables_consistent)

# Summarize missingness within each PBW/PFVC tercile.
table1_missingness_by_tercile <- purrr::map_dfr(table1_continuous_vars, function(v) {
  table1_source %>%
    group_by(pbwpfvc_tercile) %>%
    summarise(
      n_total = n(), n_nonmissing = sum(!is.na(.data[[v]])), n_missing = sum(is.na(.data[[v]])),
      pct_missing = 100 * mean(is.na(.data[[v]])), .groups = "drop"
    ) %>%
    mutate(variable = v, .before = 1)
})

missingness_reconciliation <- table1_missingness_by_tercile %>%
  group_by(variable) %>%
  summarise(
    n_total_from_terciles = sum(n_total), n_missing_from_terciles = sum(n_missing), .groups = "drop"
  ) %>%
  left_join(
    continuous_qc %>%
      select(variable, n_total_overall = n_total, n_missing_overall = n_missing),
    by = "variable"
  ) %>%
  mutate(
    total_n_matches = n_total_from_terciles == n_total_overall,
    missing_n_matches = n_missing_from_terciles == n_missing_overall
  )

missingness_by_tercile_reconciles <- all(
  missingness_reconciliation$total_n_matches & missingness_reconciliation$missing_n_matches
)


stopifnot(missingness_by_tercile_reconciles)

table1_checks <- list(
  required_variables_present = length(missing_table1_vars) == 0L,
  table1_cohort_nonempty = n_table1 > 0L,
  pbwpfvc_all_finite = all(is.finite(table1_source$pbwpfvc)),
  tercile_assignment_complete = !any(is.na(table1_source$pbwpfvc_tercile)),
  tercile_counts_sum_to_table1_n = sum(tercile_counts$n) == nrow(table1_source),
  sex_categories_valid = length(unexpected_sex) == 0L,
  race_categories_valid = length(unexpected_race) == 0L,
  mortality_coding_valid = mortality_coding_valid,
  continuous_variables_numeric = continuous_variables_numeric,
  continuous_variables_have_no_infinite_values = continuous_variables_finite,
  age_nonnegative = mathematical_constraint_checks$age_nonnegative,
  height_positive = mathematical_constraint_checks$height_positive,
  pbw_positive = mathematical_constraint_checks$pbw_positive,
  pbw_nchs_positive = mathematical_constraint_checks$pbw_nchs_positive,
  abw_positive = mathematical_constraint_checks$abw_positive,
  pfvc_positive = mathematical_constraint_checks$pfvc_positive,
  vtpbw_positive = mathematical_constraint_checks$vtpbw_positive,
  vtnchs_positive = mathematical_constraint_checks$vtnchs_positive,
  vtabw_positive = mathematical_constraint_checks$vtabw_positive,
  vtpfvc_positive = mathematical_constraint_checks$vtpfvc_positive,
  pbwpfvc_positive = mathematical_constraint_checks$pbwpfvc_positive,
  nchspfvc_positive = mathematical_constraint_checks$nchspfvc_positive,
  abwpfvc_positive = mathematical_constraint_checks$abwpfvc_positive,
  vfd_28_within_0_to_28 = mathematical_constraint_checks$vfd_28_within_0_to_28,
  derived_variables_consistent = derived_variables_consistent,
  missingness_by_tercile_reconciles = missingness_by_tercile_reconciles
)

table1_anchor_provenance_qc <- table1_source %>%
  count(intubation_identification_method, has_31500_in_hospitalization, name = "n_hospitalizations")

table1_qc <- list(
  metadata = list(
    site = site_name, generated_at = Sys.time(),
    source_file = file.path(output_dir, "analysis_cross_sectional.parquet"),
    qc_scope = "Table 1 analysis and presentation", cohort_construction_qc_included = FALSE
  ),
  upstream_input_contract = list(
    summary = script04_input_contract_summary,
    missing_cross_sectional_columns = missing_cross_sectional_contract_columns,
    missing_all_timepoint_columns = missing_all_timepoint_contract_columns,
    checks = script04_input_contract_checks_table
  ),
  respiratory_anchor_provenance = table1_anchor_provenance_qc,
  selected_tidal_volume_cleaning = list(
    policy = script03_selected_tv_qc$policy,
    thresholds = script03_selected_tv_qc$thresholds,
    upstream_retention = script03_selected_tv_qc$retention,
    checks = script03_selected_tv_qc$checks
  ),
  required_variables = list(
    required = table1_required_vars,
    missing = missing_table1_vars
  ),
  pbwpfvc_eligibility = pbwpfvc_qc,
  table1_cohort = tibble(
    n_cross_sectional = n_cross_sectional,
    n_table1 = n_table1,
    n_excluded_missing_or_nonfinite_pbwpfvc = n_excluded_table1
  ),
  terciles = list(
    counts = tercile_counts,
    ranges = tercile_ranges,
    tied_values_split_across_terciles = tercile_ties,
    n_tied_values_split_across_terciles = nrow(tercile_ties)
  ),
  categorical_variables = list(
    sex_counts = sex_counts,
    race_counts = race_counts,
    unexpected_sex_categories = as.character(unexpected_sex),
    unexpected_race_categories = as.character(unexpected_race)
  ),
  mortality = list(
    counts = mortality_counts,
    validation = mortality_validation,
    observed_values = mortality_values
  ),
  continuous_variables = list(
    missingness_and_finiteness = continuous_qc,
    missingness_by_tercile = table1_missingness_by_tercile,
    missingness_reconciliation = missingness_reconciliation,
    descriptive_ranges = continuous_range_qc,
    mathematical_constraints = mathematical_constraint_checks,
    derived_variable_consistency = derivation_qc,
    all_missing_variables = continuous_qc %>%
      filter(n_nonmissing == 0L) %>%
      pull(variable),
    variables_with_missingness = continuous_qc %>%
      filter(n_missing > 0L) %>%
      arrange(desc(pct_missing))
  ),
  checks = c(
    set_names(
      as.list(script04_input_contract_checks_table$passed),
      script04_input_contract_checks_table$check
    ),
    table1_checks
  )
)

duplicate_table1_qc_check_names <- unique(
  names(table1_qc$checks)[duplicated(names(table1_qc$checks))]
)

if (length(duplicate_table1_qc_check_names) > 0L) {
  stop(
    "table1_qc contains duplicate check name(s): ",
    paste(duplicate_table1_qc_check_names, collapse = ", ")
  )
}

save_table1_qc(table1_qc, table1_qc_path)

message("Table 1 QC saved to: ", table1_qc_path)

# Shared pediatric age groups for downstream descriptive diagnostics.
# Intervals include their lower boundary and exclude their upper boundary.
# These groups do not replace continuous age in the regression adjustment sets.
cross_sectional <- cross_sectional %>%
  mutate(
    age_group = cut(
      age_at_admission,
      breaks = c(3, 6, 9, 13, Inf),
      labels = c("3-5", "6-8", "9-12", "13-17"),
      right = FALSE
    )
  )

# Reference categories are male and white. These were set in Script 03 but are reapplied here so the references remain
# explicit and robust to the Parquet round trip. SF ratio is scaled per 10 units for downstream regression models.
cross_sectional <- cross_sectional %>%
  mutate(
    sex_category = factor(sex_category, levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    sf10 = sf_ratio / 10
  )

# =============================================================================
# 4a. Table 1
# =============================================================================

table1_data <- table1_source %>%
  select(
    pbwpfvc_tercile, age_at_admission, sex_category, race_category, height_cm,
    pbw, pbw_nchs, abw, pfvc, sofa_total, sf_ratio, mortality_event_60,
    vtpbw, vtnchs, vtabw, vtpfvc, pbwpfvc, nchspfvc, abwpfvc,
    crs, ers, vfd_28
  ) %>%
  mutate(
    mortality_event_60 = factor(
      mortality_event_60,
      levels = c(0, 1),
      labels = c("Alive at 60 days", "Deceased by 60 days")
    )
  )

table1 <- table1_data %>%
  tbl_summary(
    by = pbwpfvc_tercile, missing = "ifany", missing_text = "Missing",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
    label = list(
      age_at_admission ~ "Age (years)",
      sex_category ~ "Sex",
      race_category ~ "Race",
      height_cm ~ "Height (cm)",
      pbw ~ "McLaren predicted body weight (kg)",
      pbw_nchs ~ "NCHS Standards weight (kg)",
      abw ~ "Actual body weight (kg)",
      pfvc ~ "Predicted FVC (L)",
      sofa_total ~ "pSOFA total",
      sf_ratio ~ "SF ratio",
      mortality_event_60 ~ "60-day mortality status",
      vtpbw ~ "VT/McLaren PBW (mL/kg)",
      vtnchs ~ "VT/NCHS weight (mL/kg)",
      vtabw ~ "VT/ABW (mL/kg)",
      vtpfvc ~ "VT/PFVC (%)",
      pbwpfvc ~ "McLaren PBW/PFVC (kg/L)",
      nchspfvc ~ "NCHS weight/PFVC (kg/L)",
      abwpfvc ~ "ABW/PFVC (kg/L)",
      crs ~ "Compliance (mL/cmH2O)",
      ers ~ "Elastance (cmH2O/L)",
      vfd_28 ~ "28-day ventilator-free days"
    )
  ) %>%
  add_p(
    test = list(all_continuous() ~ "kruskal.test"),
    pvalue_fun = label_style_pvalue(digits = 3)
  ) %>%
  add_overall(last = FALSE) %>%
  bold_labels() %>%
  modify_spanning_header(
    all_stat_cols(stat_0 = FALSE) ~ "**Predicted body weight / predicted FVC (PBW/PFVC), tercile**"
  ) %>%
  modify_header(
    label ~ "**Characteristic**",
    stat_0 ~ "**Overall**, N = {N}",
    p.value ~ "**P value**"
  ) %>%
  modify_caption(paste0("**Table 1. Patient characteristics by PBW/PFVC tercile — ", site_name, "**"))

table1_gt <- table1 %>%
  as_gt() %>%
  tab_source_note(
    source_note = paste0(
      "Continuous variables are reported as median (25th, 75th percentile); categorical variables are reported as n (%). ",
      "Continuous-variable P values use the Kruskal-Wallis test. Missing values are displayed when present. ",
      "PBW/PFVC terciles were created using dplyr::ntile()."
    )
  )

table1_output_path <- file.path(table_dir, paste0("table1_", site_name, ".html"))

gt::gtsave(table1_gt, table1_output_path)

table1_output_exists <- file.exists(table1_output_path)
table1_output_size_bytes <- if (table1_output_exists) {
  file.info(table1_output_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

stopifnot(
  table1_output_exists,
  is.finite(table1_output_size_bytes),
  table1_output_size_bytes > 0
)

message("Table 1 generated and saved to: ", table1_output_path)

# Save the final Table 1 structure, statistical tests, rendered values, and output diagnostics in the Table 1 QC RDS.
table1_test_qc <- table1$table_body %>%
  filter(!is.na(p.value)) %>%
  select(any_of(c("variable", "label", "test_name", "test", "p.value")))

table1_rendered_qc <- table1$table_body %>%
  select(any_of(c(
    "variable", "var_type", "row_type", "label", "stat_0", "stat_1", "stat_2", "stat_3", "p.value"
  )))

table1_qc$final_table <- list(
  output_path = table1_output_path, output_exists = table1_output_exists,
  output_size_bytes = table1_output_size_bytes, n_source_patients = nrow(table1_data),
  tests = table1_test_qc, rendered_body = table1_rendered_qc
)

table1_qc$checks$table1_output_exists <- table1_output_exists
table1_qc$checks$table1_output_nonempty <- table1_output_exists &&
  is.finite(table1_output_size_bytes) &&
  table1_output_size_bytes > 0
table1_qc$checks$table1_overall_n_matches_source <- nrow(table1_data) == n_table1
table1_qc$metadata$updated_at <- Sys.time()

save_table1_qc(table1_qc, table1_qc_path)

message("Final Table 1 diagnostics added to: ", table1_qc_path)

# =============================================================================
# 4b. Model specifications
# =============================================================================

# The primary adjustment set includes age in years, sex, race, pSOFA, and SF ratio per 10 units. The minimally adjusted
# set excludes demographics but retains illness severity. These names describe statistical adjustment and do not, by
# themselves, imply that either set is a sufficient causal adjustment set.
covariates <- "race_category + age_at_admission + sex_category + sofa_total + sf10"
covariates_unadj <- "sofa_total + sf10"

# Readable labels for covariates in rendered regression tables. Age is modeled per year; SF ratio is modeled per 10 units.
covar_labels <- list(age_at_admission ~ "Age at admission (per year)", sf10 ~ "SF ratio (per 10)")

# =============================================================================
# Model QC: adjustment covariates
# =============================================================================
# This QC file is separate from the cohort-construction QC and Table 1 QC files.

# Detailed model diagnostics remain local under diagnostics/models.
model_diagnostic_dir <- config$model_diagnostic_dir

dir.create(model_diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(model_diagnostic_dir)) {
  stop("Script 04 could not create the model diagnostic directory: ", model_diagnostic_dir)
}

model_qc_path <- file.path(
  model_diagnostic_dir,
  paste0("model_qc_", site_name, ".rds")
)

save_model_qc <- function(qc, path) {
  if (!is.list(qc) || !is.list(qc$checks) || length(qc$checks) == 0L) {
    stop("Model QC must contain a nonempty list of integrity checks.")
  }
  
  if (any(lengths(qc$checks) != 1L)) {
    stop("Every model integrity check must return exactly one value.")
  }
  
  check_names <- names(qc$checks)
  
  if (
    is.null(check_names) ||
    anyNA(check_names) ||
    any(!nzchar(trimws(check_names))) ||
    anyDuplicated(check_names)
  ) {
    stop("Model integrity checks must have unique, nonempty names.")
  }
  
  check_values <- vapply(qc$checks, isTRUE, logical(1))
  
  if (!all(check_values)) {
    stop(
      "Model QC cannot be saved as passed. Failed or indeterminate check(s): ",
      paste(check_names[!check_values], collapse = ", ")
    )
  }
  
  saveRDS(qc, path, compress = "xz")
  
  if (!file.exists(path)) {
    stop("Model QC RDS was not created: ", path)
  }
  
  qc_readback <- readRDS(path)
  
  if (!identical(qc_readback, qc)) {
    section_names <- union(names(qc), names(qc_readback))
    differing_sections <- section_names[
      !vapply(
        section_names,
        function(section) identical(qc[[section]], qc_readback[[section]]),
        logical(1)
      )
    ]
    
    comparison <- all.equal(qc_readback, qc)
    
    message(
      "Model-QC differing top-level sections: ",
      paste(differing_sections, collapse = ", ")
    )
    print(comparison)
    
    # Inspect the live objects before evidence serialization can normalize them.
    describe_live_mismatch <- function(expected, observed, label) {
      if (identical(expected, observed)) {
        return(invisible(NULL))
      }
      
      message("Live mismatch at: ", label)
      message(
        "Types: expected=", typeof(expected),
        "; readback=", typeof(observed)
      )
      message(
        "Classes: expected=", paste(class(expected), collapse = "/"),
        "; readback=", paste(class(observed), collapse = "/")
      )
      message(
        "Attributes identical: ",
        identical(attributes(expected), attributes(observed))
      )
      message(
        "Identical ignoring environment identity: ",
        identical(expected, observed, ignore.environment = TRUE)
      )
      message(
        "Identical ignoring bytecode: ",
        identical(expected, observed, ignore.bytecode = TRUE)
      )
      
      print(all.equal(observed, expected, tolerance = 0))
      
      if (
        is.list(expected) && is.list(observed) &&
        length(expected) == length(observed) &&
        identical(names(expected), names(observed))
      ) {
        for (index in seq_along(expected)) {
          component <- if (is.null(names(expected))) {
            as.character(index)
          } else {
            names(expected)[[index]]
          }
          
          describe_live_mismatch(
            expected[[index]],
            observed[[index]],
            paste0(label, "$", component)
          )
        }
      }
    }
    
    for (section in differing_sections) {
      describe_live_mismatch(
        qc[[section]],
        qc_readback[[section]],
        section
      )
    }
    
    evidence <- list(
      generated_at = Sys.time(),
      target_path = path,
      differing_top_level_sections = differing_sections,
      comparison = comparison,
      expected_object = qc,
      readback_object = qc_readback
    )
    
    evidence_path <- tempfile(
      pattern = "model_qc_readback_failure_",
      tmpdir = dirname(path),
      fileext = ".rds"
    )
    
    evidence_saved <- tryCatch(
      {
        saveRDS(evidence, evidence_path, compress = "xz")
        TRUE
      },
      error = function(error) {
        message(
          "Could not save model-QC read-back failure evidence: ",
          conditionMessage(error)
        )
        FALSE
      }
    )
    
    if (evidence_saved) {
      message(
        "Model-QC failure evidence written, not independently read-back validated: ",
        evidence_path
      )
    }
    
    stop("Model QC RDS failed its complete read-back comparison: ", path)
  }
  
  invisible(path)
}

model_covariate_vars <- c(
  "race_category", "age_at_admission", "sex_category", "sofa_total", "sf_ratio", "sf10", "bmi"
)

missing_model_covariates <- setdiff(model_covariate_vars, names(cross_sectional))

stopifnot(length(missing_model_covariates) == 0L)

model_covariate_qc <- purrr::map_dfr(model_covariate_vars, function(v) {
  x <- cross_sectional[[v]]
  x_numeric <- if (is.numeric(x)) x else NULL
  
  tibble(
    variable = v, class = paste(class(x), collapse = "/"), n_total = length(x),
    n_nonmissing = sum(!is.na(x)), n_missing = sum(is.na(x)), pct_missing = 100 * mean(is.na(x)),
    n_unique_nonmissing = n_distinct(x, na.rm = TRUE),
    n_infinite = if (is.numeric(x)) sum(is.infinite(x), na.rm = TRUE) else NA_integer_,
    minimum = if (is.numeric(x) && any(is.finite(x_numeric))) {
      min(x_numeric[is.finite(x_numeric)])
    } else {
      NA_real_
    },
    maximum = if (is.numeric(x) && any(is.finite(x_numeric))) {
      max(x_numeric[is.finite(x_numeric)])
    } else {
      NA_real_
    }
  )
})

sf10_compared <- is.finite(cross_sectional$sf10) & is.finite(cross_sectional$sf_ratio)

sf10_max_absolute_difference <- if (any(sf10_compared)) {
  max(abs(cross_sectional$sf10[sf10_compared] - cross_sectional$sf_ratio[sf10_compared] / 10))
} else {
  NA_real_
}

sf10_derivation_valid <- any(sf10_compared) &&
  is.finite(sf10_max_absolute_difference) &&
  sf10_max_absolute_difference <= 1e-12

model_covariate_complete_cases <- tibble(
  adjustment_set = c("Primary adjusted", "Minimally adjusted", "Primary adjusted plus BMI"),
  variables = c(
    "race_category + age_at_admission + sex_category + sofa_total + sf10",
    "sofa_total + sf10",
    "race_category + age_at_admission + sex_category + sofa_total + sf10 + bmi"
  ),
  n_complete = c(
    sum(complete.cases(cross_sectional[c(
      "race_category", "age_at_admission", "sex_category", "sofa_total", "sf10"
    )])),
    sum(complete.cases(cross_sectional[c("sofa_total", "sf10")])),
    sum(complete.cases(cross_sectional[c(
      "race_category", "age_at_admission", "sex_category", "sofa_total", "sf10", "bmi"
    )]))
  )
) %>%
  mutate(
    n_total = nrow(cross_sectional), n_excluded = n_total - n_complete,
    pct_complete = 100 * n_complete / n_total
  )

model_covariate_checks <- list(
  required_covariates_present = length(missing_model_covariates) == 0L,
  numeric_covariates_numeric = all(purrr::map_lgl(
    cross_sectional[c("age_at_admission", "sofa_total", "sf_ratio", "sf10", "bmi")],
    is.numeric
  )),
  numeric_covariates_have_no_infinite_values = !any(model_covariate_qc$n_infinite > 0, na.rm = TRUE),
  covariates_have_variation = all(model_covariate_qc$n_unique_nonmissing >= 2L),
  sex_levels_valid = identical(levels(cross_sectional$sex_category), c("Male", "Female")),
  race_levels_valid = identical(levels(cross_sectional$race_category), c("WHITE", "BLACK", "OTHER")),
  sf10_derivation_valid = sf10_derivation_valid,
  primary_adjustment_set_has_complete_cases = model_covariate_complete_cases$n_complete[
    model_covariate_complete_cases$adjustment_set == "Primary adjusted"
  ] > 0L,
  bmi_adjustment_set_has_complete_cases = model_covariate_complete_cases$n_complete[
    model_covariate_complete_cases$adjustment_set == "Primary adjusted plus BMI"
  ] > 0L
)



stopifnot(isTRUE(all(unlist(model_covariate_checks))))

model_qc <- list(
  metadata = list(
    site = site_name, generated_at = Sys.time(),
    source_file = file.path(output_dir, "analysis_cross_sectional.parquet"),
    qc_scope = "Regression model specification and fitting",
    cohort_construction_qc_included = FALSE, table1_qc_included = FALSE
  ),
  upstream_input_contract = list(
    summary = script04_input_contract_summary,
    missing_cross_sectional_columns = missing_cross_sectional_contract_columns,
    missing_all_timepoint_columns = missing_all_timepoint_contract_columns,
    checks = script04_input_contract_checks_table
  ),
  adjustment_covariates = list(
    primary_formula = covariates, minimally_adjusted_formula = covariates_unadj,
    required_variables = model_covariate_vars, missing_variables = missing_model_covariates,
    variable_summary = model_covariate_qc, complete_case_summary = model_covariate_complete_cases,
    sf10_derivation = tibble(
      n_compared = sum(sf10_compared), max_absolute_difference = sf10_max_absolute_difference,
      passed = sf10_derivation_valid
    )
  ),
  checks = c(
    set_names(
      as.list(script04_input_contract_checks_table$passed),
      script04_input_contract_checks_table$check
    ),
    model_covariate_checks
  )
)

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Initial model QC saved to: ", model_qc_path)

# Per the original paper, any model whose outcome or exposure is derived from driving pressure is additionally adjusted
# for BMI. Models without a driving-pressure component use the standard covariate set.
DP_DERIVED <- c(
  "dp", "ers", "crs", "ers_pbw", "ers_pfvc", "mechanical_power", "mp_pbw", "mp_pfvc"
)

uses_dp <- function(...) {
  vars <- trimws(unlist(strsplit(paste(c(...), collapse = " + "), "\\+")))
  any(vars %in% DP_DERIVED)
}

model_covariates <- function(..., adjusted = TRUE) {
  base <- if (adjusted) covariates else covariates_unadj
  if (uses_dp(...)) paste(base, "+ bmi") else base
}

# Define all exposure specifications.
exposure_specs <- list(
  vtpbw = "vtpbw",
  vtpfvc_vtpbw = "vtpfvc + vtpbw",
  vtpbw_pfvc = "vtpbw + pfvc",
  vtpbw_pbwpfvc = "vtpbw + pbwpfvc",
  
  vtnchs = "vtnchs",
  vtpfvc_vtnchs = "vtpfvc + vtnchs",
  vtnchs_pfvc = "vtnchs + pfvc",
  vtnchs_nchspfvc = "vtnchs + nchspfvc",
  
  vtabw = "vtabw",
  vtpfvc_vtabw = "vtpfvc + vtabw",
  vtabw_pfvc = "vtabw + pfvc",
  vtabw_abwpfvc = "vtabw + abwpfvc",
  
  vtpfvc = "vtpfvc"
)

exposure_labels <- c(
  vtpbw = "VT/McLaren PBW",
  vtpfvc_vtpbw = "VT/PFVC + VT/McLaren PBW",
  vtpbw_pfvc = "VT/McLaren PBW + PFVC",
  vtpbw_pbwpfvc = "VT/McLaren PBW + McLaren PBW/PFVC",
  
  vtnchs = "VT/NCHS weight",
  vtpfvc_vtnchs = "VT/PFVC + VT/NCHS weight",
  vtnchs_pfvc = "VT/NCHS weight + PFVC",
  vtnchs_nchspfvc = "VT/NCHS weight + NCHS weight/PFVC",
  
  vtabw = "VT/ABW",
  vtpfvc_vtabw = "VT/PFVC + VT/ABW",
  vtabw_pfvc = "VT/ABW + PFVC",
  vtabw_abwpfvc = "VT/ABW + ABW/PFVC",
  
  vtpfvc = "VT/PFVC"
)

# =============================================================================
# Model QC: exposure specifications
# =============================================================================

exposure_vars <- unique(unlist(purrr::map(exposure_specs, function(spec) {
  all.vars(as.formula(paste("~", spec)))
})))

missing_exposure_vars <- setdiff(exposure_vars, names(cross_sectional))

stopifnot(length(missing_exposure_vars) == 0L)

exposure_variable_qc <- purrr::map_dfr(exposure_vars, function(v) {
  x <- cross_sectional[[v]]
  finite_values <- x[is.finite(x)]
  
  tibble(
    variable = v, class = paste(class(x), collapse = "/"), n_total = length(x),
    n_nonmissing = sum(!is.na(x)), n_missing = sum(is.na(x)), pct_missing = 100 * mean(is.na(x)),
    n_unique_nonmissing = n_distinct(x, na.rm = TRUE), n_infinite = sum(is.infinite(x), na.rm = TRUE),
    minimum = if (length(finite_values) > 0L) min(finite_values) else NA_real_,
    median = if (length(finite_values) > 0L) median(finite_values) else NA_real_,
    maximum = if (length(finite_values) > 0L) max(finite_values) else NA_real_
  )
})

exposure_specification_qc <- purrr::imap_dfr(exposure_specs, function(spec, spec_name) {
  terms <- all.vars(as.formula(paste("~", spec)))
  adjustment_vars <- all.vars(as.formula(paste("~", model_covariates(spec))))
  matrix_vars <- unique(c(terms, adjustment_vars))
  
  analysis_data <- cross_sectional %>%
    select(all_of(matrix_vars)) %>%
    drop_na()
  
  stopifnot(nrow(analysis_data) > 0L)
  
  model_matrix <- model.matrix(
    as.formula(paste("~", spec, "+", model_covariates(spec))),
    data = analysis_data
  )
  
  tibble(
    model_spec = spec_name, label = unname(exposure_labels[[spec_name]]), formula = spec,
    exposure_terms = paste(terms, collapse = " + "), n_exposure_terms = length(terms),
    n_complete = nrow(analysis_data), n_model_matrix_columns = ncol(model_matrix),
    model_matrix_rank = qr(model_matrix)$rank, full_rank = qr(model_matrix)$rank == ncol(model_matrix)
  )
})

exposure_correlation_matrix <- cross_sectional %>%
  select(all_of(exposure_vars)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson")

exposure_correlation <- exposure_correlation_matrix %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  as_tibble()

exposure_correlation_long <- exposure_correlation_matrix %>%
  as.data.frame() %>%
  rownames_to_column("variable_1") %>%
  pivot_longer(-variable_1, names_to = "variable_2", values_to = "correlation") %>%
  filter(variable_1 < variable_2) %>%
  arrange(desc(abs(correlation))) %>%
  as_tibble()

maximum_absolute_exposure_correlation <- if (nrow(exposure_correlation_long) > 0L) {
  max(abs(exposure_correlation_long$correlation), na.rm = TRUE)
} else {
  NA_real_
}

exposure_checks <- list(
  required_exposure_variables_present = length(missing_exposure_vars) == 0L,
  exposure_variables_numeric = all(purrr::map_lgl(cross_sectional[exposure_vars], is.numeric)),
  exposure_variables_have_no_infinite_values = !any(exposure_variable_qc$n_infinite > 0L),
  exposure_variables_have_variation = all(exposure_variable_qc$n_unique_nonmissing >= 2L),
  exposure_specifications_have_complete_cases = all(exposure_specification_qc$n_complete > 0L),
  exposure_model_matrices_full_rank = all(exposure_specification_qc$full_rank),
  exposure_labels_match_specifications = identical(names(exposure_specs), names(exposure_labels))
)


stopifnot(isTRUE(all(unlist(exposure_checks))))

model_qc$exposures <- list(
  specifications = exposure_specs, labels = exposure_labels, required_variables = exposure_vars,
  missing_variables = missing_exposure_vars, variable_summary = exposure_variable_qc,
  specification_summary = exposure_specification_qc, correlation_matrix = exposure_correlation,
  pairwise_correlations = exposure_correlation_long,
  maximum_absolute_correlation = maximum_absolute_exposure_correlation
)

model_qc$checks <- c(model_qc$checks, exposure_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Exposure-specification QC added to: ", model_qc_path)

# =============================================================================
# 4c. Logistic regression — 60-day mortality
# =============================================================================

mortality_outcome <- "mortality_event_60"

stopifnot(mortality_outcome %in% names(cross_sectional))

mortality_values <- sort(unique(na.omit(cross_sectional[[mortality_outcome]])))
has_mortality_variation <- length(mortality_values) > 1L

mortality_outcome_qc <- tibble(
  outcome = mortality_outcome, n_total = nrow(cross_sectional),
  n_nonmissing = sum(!is.na(cross_sectional[[mortality_outcome]])),
  n_missing = sum(is.na(cross_sectional[[mortality_outcome]])),
  n_events = sum(cross_sectional[[mortality_outcome]] == 1, na.rm = TRUE),
  n_nonevents = sum(cross_sectional[[mortality_outcome]] == 0, na.rm = TRUE),
  observed_values = paste(mortality_values, collapse = ", "),
  coding_valid = all(mortality_values %in% c(0, 1)),
  has_variation = has_mortality_variation
)


mortality_outcome_checks <- list(
  mortality_outcome_present = mortality_outcome %in% names(cross_sectional),
  mortality_outcome_coding_valid = mortality_outcome_qc$coding_valid,
  mortality_outcome_has_variation = mortality_outcome_qc$has_variation,
  mortality_outcome_has_events = mortality_outcome_qc$n_events > 0L,
  mortality_outcome_has_nonevents = mortality_outcome_qc$n_nonevents > 0L
)

stopifnot(isTRUE(all(unlist(mortality_outcome_checks))))

model_qc$mortality <- list(
  outcome = mortality_outcome,
  outcome_summary = mortality_outcome_qc
)

model_qc$checks <- c(model_qc$checks, mortality_outcome_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Mortality-outcome QC added to: ", model_qc_path)

# Use one common complete-case sample for all adjusted mortality models so likelihoods and AIC values are comparable.
mortality_analysis_vars <- unique(c(
  mortality_outcome, exposure_vars, all.vars(as.formula(paste("~", covariates)))
))

mortality_analysis_data <- cross_sectional %>%
  mutate(.source_row = row_number()) %>%
  select(.source_row, all_of(mortality_analysis_vars)) %>%
  drop_na(all_of(mortality_analysis_vars))

mortality_sample_qc <- mortality_analysis_data %>%
  summarise(
    n_total_source = nrow(cross_sectional), n_analyzed = n(),
    n_excluded = n_total_source - n_analyzed, pct_analyzed = 100 * n_analyzed / n_total_source,
    n_events = sum(.data[[mortality_outcome]] == 1L),
    n_nonevents = sum(.data[[mortality_outcome]] == 0L),
    event_percent = 100 * mean(.data[[mortality_outcome]] == 1L)
  )

mortality_exclusion_qc <- cross_sectional %>%
  mutate(
    included_in_mortality_analysis = complete.cases(across(all_of(mortality_analysis_vars))),
    mortality_status = case_when(
      .data[[mortality_outcome]] == 0L ~ "Alive at 60 days",
      .data[[mortality_outcome]] == 1L ~ "Deceased by 60 days",
      TRUE ~ "Missing"
    )
  ) %>%
  count(mortality_status, included_in_mortality_analysis, name = "n") %>%
  group_by(mortality_status) %>%
  mutate(percent_within_status = 100 * n / sum(n)) %>%
  ungroup()

mortality_sample_checks <- list(
  common_sample_nonempty = nrow(mortality_analysis_data) > 0L,
  common_sample_has_events = mortality_sample_qc$n_events > 0L,
  common_sample_has_nonevents = mortality_sample_qc$n_nonevents > 0L,
  common_sample_has_complete_data = all(complete.cases(mortality_analysis_data)),
  common_sample_outcome_is_binary = all(mortality_analysis_data[[mortality_outcome]] %in% c(0L, 1L))
)


stopifnot(isTRUE(all(unlist(mortality_sample_checks))))

model_qc$mortality$analysis_sample <- list(
  required_variables = mortality_analysis_vars,
  sample_summary = mortality_sample_qc,
  inclusion_by_outcome = mortality_exclusion_qc
)

model_qc$checks <- c(model_qc$checks, mortality_sample_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Mortality common-sample QC added to: ", model_qc_path)

if (has_mortality_variation) {
  mortality_models <- map(exposure_specs, function(exposure_spec) {
    formula_str <- paste(mortality_outcome, "~", exposure_spec, "+", model_covariates(exposure_spec))
    glm(as.formula(formula_str), data = mortality_analysis_data, family = binomial)
  })
  
  mortality_model_qc <- purrr::imap_dfr(mortality_models, function(model, model_name) {
    outcome_values <- model.response(model.frame(model))
    model_matrix <- model.matrix(model)
    coefficient_table <- summary(model)$coefficients
    n_events <- sum(outcome_values == 1L)
    n_nonevents <- sum(outcome_values == 0L)
    n_parameters_excluding_intercept <- model$rank - 1L
    
    tibble(
      model_spec = model_name, label = unname(exposure_labels[[model_name]]),
      formula = paste(deparse(formula(model)), collapse = " "), n_obs = nobs(model),
      n_events = n_events, n_nonevents = n_nonevents,
      n_parameters_excluding_intercept = n_parameters_excluding_intercept,
      events_per_parameter = n_events / n_parameters_excluding_intercept,
      converged = isTRUE(model$converged), boundary = isTRUE(model$boundary),
      iterations = model$iter, model_matrix_columns = ncol(model_matrix),
      model_matrix_rank = qr(model_matrix)$rank,
      full_rank = qr(model_matrix)$rank == ncol(model_matrix),
      coefficients_finite = all(is.finite(coef(model))),
      standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
      maximum_absolute_coefficient = max(abs(coef(model))),
      minimum_fitted_probability = min(fitted(model)),
      maximum_fitted_probability = max(fitted(model)),
      n_extreme_fitted_probabilities = sum(fitted(model) < 1e-6 | fitted(model) > 1 - 1e-6),
      log_likelihood = as.numeric(logLik(model)), AIC = AIC(model)
    )
  })
  
  mortality_coefficient_qc <- purrr::imap_dfr(mortality_models, function(model, model_name) {
    coefficient_table <- summary(model)$coefficients
    
    as_tibble(coefficient_table, rownames = "term") %>%
      transmute(
        model_spec = model_name, label = unname(exposure_labels[[model_name]]), term,
        coefficient = Estimate, standard_error = `Std. Error`, z_value = `z value`,
        p_value = `Pr(>|z|)`, odds_ratio = exp(coefficient),
        conf_low = exp(coefficient - 1.96 * standard_error),
        conf_high = exp(coefficient + 1.96 * standard_error)
      )
  })
  
  mortality_same_analysis_rows <- all(purrr::map_lgl(
    mortality_models,
    ~ identical(rownames(model.frame(.x)), rownames(model.frame(mortality_models[[1]])))
  ))
  
  mortality_model_checks <- list(
    all_models_fitted = length(mortality_models) == length(exposure_specs),
    all_models_use_common_n = all(mortality_model_qc$n_obs == nrow(mortality_analysis_data)),
    all_models_use_same_rows = mortality_same_analysis_rows,
    all_models_converged = all(mortality_model_qc$converged),
    no_models_at_boundary = !any(mortality_model_qc$boundary),
    all_model_matrices_full_rank = all(mortality_model_qc$full_rank),
    all_coefficients_finite = all(mortality_model_qc$coefficients_finite),
    all_standard_errors_finite = all(mortality_model_qc$standard_errors_finite),
    no_extreme_fitted_probabilities = all(mortality_model_qc$n_extreme_fitted_probabilities == 0L)
  )
  
  # Evaluate influence without automatically excluding observations that exceed heuristic thresholds.
  mortality_influence_qc <- purrr::imap_dfr(mortality_models, function(model, model_name) {
    n <- nobs(model)
    p <- model$rank
    cooks_values <- cooks.distance(model)
    leverage_values <- hatvalues(model)
    standardized_residuals <- rstandard(model)
    dfbeta_values <- dfbetas(model)
    maximum_absolute_dfbeta <- apply(abs(dfbeta_values), 1L, max, na.rm = TRUE)
    
    cooks_threshold <- 4 / n
    leverage_threshold <- 2 * p / n
    dfbeta_threshold <- 2 / sqrt(n)
    
    tibble(
      model_spec = model_name, label = unname(exposure_labels[[model_name]]),
      analysis_row = seq_len(n), source_row = mortality_analysis_data$.source_row,
      outcome = model.response(model.frame(model)),
      cooks_distance = unname(cooks_values), cooks_threshold = cooks_threshold,
      high_cooks_distance = cooks_distance > cooks_threshold,
      leverage = unname(leverage_values), leverage_threshold = leverage_threshold,
      high_leverage = leverage > leverage_threshold,
      standardized_residual = unname(standardized_residuals),
      large_standardized_residual = abs(standardized_residual) > 3,
      maximum_absolute_dfbeta = maximum_absolute_dfbeta, dfbeta_threshold = dfbeta_threshold,
      high_dfbeta = maximum_absolute_dfbeta > dfbeta_threshold,
      flagged_by_any_measure = high_cooks_distance | high_leverage | large_standardized_residual | high_dfbeta
    )
  })
  
  mortality_influence_summary <- mortality_influence_qc %>%
    group_by(model_spec, label) %>%
    summarise(
      n_obs = n(), n_high_cooks_distance = sum(high_cooks_distance),
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
  
  mortality_flagged_observations <- mortality_influence_qc %>%
    filter(flagged_by_any_measure) %>%
    arrange(model_spec, desc(cooks_distance))
  
  mortality_influence_checks <- list(
    influence_diagnostics_have_expected_rows = all(
      mortality_influence_summary$n_obs == nrow(mortality_analysis_data)
    ),
    source_rows_are_complete = !any(is.na(mortality_influence_qc$source_row)),
    source_rows_are_valid = all(mortality_influence_qc$source_row %in% seq_len(nrow(cross_sectional))),
    source_rows_unique_within_model = mortality_influence_qc %>%
      group_by(model_spec) %>%
      summarise(rows_unique = n_distinct(source_row) == n(), .groups = "drop") %>%
      pull(rows_unique) %>%
      all(),
    cooks_distances_finite = all(is.finite(mortality_influence_qc$cooks_distance)),
    leverage_values_finite = all(is.finite(mortality_influence_qc$leverage)),
    standardized_residuals_finite = all(is.finite(mortality_influence_qc$standardized_residual)),
    dfbeta_values_finite = all(is.finite(mortality_influence_qc$maximum_absolute_dfbeta))
  )
  
  
  stopifnot(
    isTRUE(all(unlist(mortality_model_checks))),
    isTRUE(all(unlist(mortality_influence_checks)))
  )
  
  model_qc$mortality$model_diagnostics <- mortality_model_qc
  model_qc$mortality$coefficient_diagnostics <- mortality_coefficient_qc
  model_qc$mortality$events_per_parameter_below_10 <- mortality_model_qc %>%
    select(model_spec, label, events_per_parameter) %>%
    mutate(below_10 = events_per_parameter < 10)
  
  model_qc$mortality$influence_diagnostics <- list(
    summary = mortality_influence_summary,
    flagged_observations = mortality_flagged_observations,
    thresholds = tibble(
      diagnostic = c("Cook's distance", "Leverage", "Absolute standardized residual", "Absolute DFBETA"),
      threshold = c("4 / n", "2 × model rank / n", "3", "2 / sqrt(n)")
    )
  )
  
  model_qc$checks <- c(model_qc$checks, mortality_model_checks, mortality_influence_checks)
  model_qc$metadata$updated_at <- Sys.time()
  
  save_model_qc(model_qc, model_qc_path)
  
  stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))
  
  message("Mortality-model diagnostics added to: ", model_qc_path)
  
  mortality_regression_labels <- c(
    vtpbw = "VT/McLaren PBW (mL/kg)",
    vtnchs = "VT/NCHS weight (mL/kg)",
    vtabw = "VT/ABW (mL/kg)",
    vtpfvc = "VT/PFVC (%)",
    pfvc = "Predicted FVC (L)",
    pbwpfvc = "McLaren PBW/PFVC (kg/L)",
    nchspfvc = "NCHS weight/PFVC (kg/L)",
    abwpfvc = "ABW/PFVC (kg/L)",
    race_category = "Race",
    age_at_admission = "Age at admission (per year)",
    sex_category = "Sex",
    sofa_total = "pSOFA total",
    sf10 = "SF ratio (per 10)"
  )
  
  mortality_tables <- map(mortality_models, function(model) {
    tbl_regression(model, exponentiate = TRUE) %>%
      modify_table_body(function(table_body) {
        table_body %>%
          mutate(
            replacement_label = unname(mortality_regression_labels[variable]),
            label = if_else(row_type == "label" & !is.na(replacement_label), replacement_label, label)
          ) %>%
          select(-replacement_label)
      }) %>%
      bold_p()
  })
  
  mortality_aic_qc <- purrr::imap_dfr(mortality_models, function(model, model_name) {
    tibble(
      model_spec = model_name, label = unname(exposure_labels[[model_name]]),
      n_obs = nobs(model), n_parameters = attr(logLik(model), "df"),
      log_likelihood = as.numeric(logLik(model)), AIC = AIC(model)
    )
  })
  
  
  mortality_table <- tbl_merge(
    mortality_tables, tab_spanner = exposure_labels, quiet = TRUE
  ) %>%
    modify_caption(paste0("**Adjusted logistic regression models for 60-day mortality — ", site_name, "**"))
  
  mortality_table_gt <- mortality_table %>%
    as_gt() %>%
    tab_source_note(
      source_note = paste0(
        "All models use the same complete-case sample: N = ", nrow(mortality_analysis_data), ", including ",
        mortality_sample_qc$n_events, " deaths by 60 days. Estimates are odds ratios with 95% confidence intervals. ",
        "Models adjust for age, sex, race, pSOFA total, and SF ratio per 10 units."
      )
    ) %>%
    tab_source_note(
      source_note = paste0(
        "Only terms included in a given exposure specification are displayed in that model column. ",
        "Blank cells indicate that the term was not included."
      )
    )
  
  mortality_output_path <- file.path(table_dir, paste0("regression_mortality_", site_name, ".html"))
  
  gt::gtsave(mortality_table_gt, mortality_output_path)
  
  mortality_output_exists <- file.exists(mortality_output_path)
  mortality_output_info <- file.info(mortality_output_path, extra_cols = FALSE)
  mortality_output_size_bytes <- if (mortality_output_exists) mortality_output_info$size[[1]] else NA_real_
  
  mortality_output_checks <- list(
    mortality_output_exists = mortality_output_exists,
    mortality_output_nonempty = mortality_output_exists &&
      is.finite(mortality_output_size_bytes) &&
      mortality_output_size_bytes > 0,
    mortality_aic_uses_common_n = all(mortality_aic_qc$n_obs == nrow(mortality_analysis_data))
  )
  
  model_qc$mortality$output <- list(
    path = mortality_output_path, exists = mortality_output_exists,
    size_bytes = mortality_output_size_bytes, AIC = mortality_aic_qc,
    common_sample_n = nrow(mortality_analysis_data), event_n = mortality_sample_qc$n_events
  )
  
  model_qc$checks <- c(model_qc$checks, mortality_output_checks)
  model_qc$metadata$updated_at <- Sys.time()
  
  save_model_qc(model_qc, model_qc_path)
  
  stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))
  
  message("Mortality regression table saved to: ", mortality_output_path)
} else {
  mortality_models <- NULL
  
  message(
    "Skipping 60-day mortality regression: `", mortality_outcome,
    "` does not contain both event and nonevent observations."
  )
}

# =============================================================================
# 4d. Linear regression — continuous outcomes
# =============================================================================

# NOTE: 28-day VFDs are not modeled here as a continuous outcome.
continuous_outcomes <- list(
  ers = list(var = "ers", label = "Elastance"),
  crs = list(var = "crs", label = "Compliance"),
  dp = list(var = "dp", label = "Driving pressure"),
  ers_pbw = list(var = "ers_pbw", label = "Ers x PBW"),
  ers_pfvc = list(var = "ers_pfvc", label = "Ers x PFVC"),
  mp = list(var = "mechanical_power", label = "Mechanical power"),
  mp_pbw = list(var = "mp_pbw", label = "MP / PBW"),
  mp_pfvc = list(var = "mp_pfvc", label = "MP / PFVC")
)

# =============================================================================
# Model QC: continuous outcomes
# =============================================================================

continuous_outcome_vars <- purrr::map_chr(continuous_outcomes, "var")
missing_continuous_outcomes <- setdiff(continuous_outcome_vars, names(cross_sectional))

stopifnot(length(missing_continuous_outcomes) == 0L)

continuous_outcome_qc <- purrr::imap_dfr(continuous_outcomes, function(outcome, outcome_name) {
  x <- cross_sectional[[outcome$var]]
  finite_values <- x[is.finite(x)]
  
  tibble(
    outcome_name = outcome_name, variable = outcome$var, label = outcome$label,
    class = paste(class(x), collapse = "/"), n_total = length(x),
    n_nonmissing = sum(!is.na(x)), n_missing = sum(is.na(x)), pct_missing = 100 * mean(is.na(x)),
    n_infinite = sum(is.infinite(x), na.rm = TRUE), n_unique_nonmissing = n_distinct(x, na.rm = TRUE),
    minimum = if (length(finite_values) > 0L) min(finite_values) else NA_real_,
    median = if (length(finite_values) > 0L) median(finite_values) else NA_real_,
    maximum = if (length(finite_values) > 0L) max(finite_values) else NA_real_,
    bmi_adjusted = uses_dp(outcome$var)
  )
})

continuous_sample_qc <- purrr::imap_dfr(continuous_outcomes, function(outcome, outcome_name) {
  covariate_formula <- model_covariates(outcome$var)
  covariate_vars <- all.vars(as.formula(paste("~", covariate_formula)))
  required_vars <- unique(c(outcome$var, exposure_vars, covariate_vars))
  
  analysis_data <- cross_sectional %>%
    select(all_of(required_vars)) %>%
    drop_na()
  
  tibble(
    outcome_name = outcome_name, variable = outcome$var, label = outcome$label,
    covariates = covariate_formula, bmi_adjusted = "bmi" %in% covariate_vars,
    n_total = nrow(cross_sectional), n_analyzed = nrow(analysis_data),
    n_excluded = n_total - n_analyzed, pct_analyzed = 100 * n_analyzed / n_total,
    n_unique_outcome = n_distinct(analysis_data[[outcome$var]])
  )
})

# Validate mathematical constraints without redefining any upstream-derived outcome.
continuous_outcome_constraint_qc <- tibble(
  outcome = c("ers", "crs", "dp", "ers_pbw", "ers_pfvc", "mechanical_power", "mp_pbw", "mp_pfvc"),
  constraint = c(
    "Greater than 0", "Greater than 0", "Greater than 0", "Greater than 0",
    "Greater than 0", "Greater than or equal to 0", "Greater than or equal to 0", "Greater than or equal to 0"
  ),
  n_nonmissing = c(
    sum(!is.na(cross_sectional$ers)),
    sum(!is.na(cross_sectional$crs)),
    sum(!is.na(cross_sectional$dp)),
    sum(!is.na(cross_sectional$ers_pbw)),
    sum(!is.na(cross_sectional$ers_pfvc)),
    sum(!is.na(cross_sectional$mechanical_power)),
    sum(!is.na(cross_sectional$mp_pbw)),
    sum(!is.na(cross_sectional$mp_pfvc))
  ),
  n_violations = c(
    sum(cross_sectional$ers <= 0, na.rm = TRUE),
    sum(cross_sectional$crs <= 0, na.rm = TRUE),
    sum(cross_sectional$dp <= 0, na.rm = TRUE),
    sum(cross_sectional$ers_pbw <= 0, na.rm = TRUE),
    sum(cross_sectional$ers_pfvc <= 0, na.rm = TRUE),
    sum(cross_sectional$mechanical_power < 0, na.rm = TRUE),
    sum(cross_sectional$mp_pbw < 0, na.rm = TRUE),
    sum(cross_sectional$mp_pfvc < 0, na.rm = TRUE)
  )
)

# Confirm that the regenerated Script 03 output retains the driving-pressure provenance variables.
dp_provenance_vars <- c("dp_raw", "dp_source_raw", "dp_source", "dp_type", "dp_invalid_reason")
missing_dp_provenance_vars <- setdiff(dp_provenance_vars, names(cross_sectional))

dp_source_qc <- if (length(missing_dp_provenance_vars) == 0L) {
  cross_sectional %>%
    count(dp_source, dp_type, name = "n", .drop = FALSE) %>%
    mutate(percent = 100 * n / sum(n))
} else {
  tibble(dp_source = character(), dp_type = character(), n = integer(), percent = numeric())
}

continuous_outcome_checks <- list(
  required_outcomes_present = length(missing_continuous_outcomes) == 0L,
  outcomes_numeric = all(purrr::map_lgl(cross_sectional[continuous_outcome_vars], is.numeric)),
  outcomes_have_no_infinite_values = !any(continuous_outcome_qc$n_infinite > 0L),
  outcomes_have_variation = all(continuous_outcome_qc$n_unique_nonmissing >= 2L),
  outcome_samples_nonempty = all(continuous_sample_qc$n_analyzed > 0L),
  analyzed_outcomes_have_variation = all(continuous_sample_qc$n_unique_outcome >= 2L),
  continuous_outcome_constraints_satisfied = all(continuous_outcome_constraint_qc$n_violations == 0L),
  dp_provenance_variables_present = length(missing_dp_provenance_vars) == 0L
)


stopifnot(isTRUE(all(unlist(continuous_outcome_checks))))

model_qc$continuous_outcomes <- list(
  definitions = continuous_outcomes,
  variable_summary = continuous_outcome_qc,
  common_sample_summary = continuous_sample_qc,
  mathematical_constraints = continuous_outcome_constraint_qc,
  driving_pressure_provenance = list(
    required_variables = dp_provenance_vars,
    missing_variables = missing_dp_provenance_vars,
    source_summary = dp_source_qc,
    derivation_note = paste0(
      "Derived upstream in Script 03 using observed plateau pressure minus set PEEP when available, ",
      "with observed PIP and set PIP minus set PEEP as prespecified pediatric fallback estimates. ",
      "Nonpositive candidate values are set to missing upstream."
    )
  )
)

model_qc$checks <- c(model_qc$checks, continuous_outcome_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Script 04 | continuous-outcome inputs validated")

if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop("The sandwich package is required for HC3 robust standard errors.")
}

tidy_lm_hc3 <- function(x, conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE, ...) {
  coefficient_estimates <- coef(x)
  robust_vcov <- sandwich::vcovHC(x, type = "HC3")
  robust_standard_errors <- sqrt(diag(robust_vcov))
  residual_df <- df.residual(x)
  critical_value <- qt(1 - (1 - conf.level) / 2, df = residual_df)
  test_statistics <- coefficient_estimates / robust_standard_errors
  p_values <- 2 * pt(abs(test_statistics), df = residual_df, lower.tail = FALSE)
  
  results <- tibble(
    term = names(coefficient_estimates), estimate = unname(coefficient_estimates),
    std.error = unname(robust_standard_errors), statistic = unname(test_statistics),
    p.value = unname(p_values),
    conf.low = unname(coefficient_estimates - critical_value * robust_standard_errors),
    conf.high = unname(coefficient_estimates + critical_value * robust_standard_errors)
  )
  
  if (isTRUE(exponentiate)) {
    results <- results %>%
      mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high))
  }
  
  results
}

continuous_regression_labels <- c(
  vtpbw = "VT/McLaren PBW (mL/kg)",
  vtnchs = "VT/NCHS weight (mL/kg)",
  vtabw = "VT/ABW (mL/kg)",
  vtpfvc = "VT/PFVC (%)",
  pfvc = "Predicted FVC (L)",
  pbwpfvc = "McLaren PBW/PFVC (kg/L)",
  nchspfvc = "NCHS weight/PFVC (kg/L)",
  abwpfvc = "ABW/PFVC (kg/L)",
  race_category = "Race",
  age_at_admission = "Age at admission (per year)",
  sex_category = "Sex",
  sofa_total = "pSOFA total",
  sf10 = "SF ratio (per 10)",
  bmi = "BMI (kg/m²)"
)

continuous_models <- list()
continuous_tables <- list()
continuous_analysis_data <- list()
continuous_model_qc <- list()
continuous_aic_qc <- list()
continuous_hc3_qc <- list()

for (outcome_name in names(continuous_outcomes)) {
  outcome_var <- continuous_outcomes[[outcome_name]]$var
  outcome_label <- continuous_outcomes[[outcome_name]]$label
  covars_used <- model_covariates(outcome_var)
  covariate_vars <- all.vars(as.formula(paste("~", covars_used)))
  required_vars <- unique(c(outcome_var, exposure_vars, covariate_vars))
  
  # Use one complete-case sample for all exposure specifications within this outcome.
  analysis_data_outcome <- cross_sectional %>%
    mutate(.source_row = row_number()) %>%
    select(.source_row, all_of(required_vars)) %>%
    drop_na(all_of(required_vars))
  
  stopifnot(
    nrow(analysis_data_outcome) > 0L,
    n_distinct(analysis_data_outcome[[outcome_var]]) >= 2L
  )
  
  message("Script 04 | continuous models and tables: ", outcome_name)
  models_for_outcome <- purrr::map(exposure_specs, function(exposure_spec) {
    formula_str <- paste(outcome_var, "~", exposure_spec, "+", covars_used)
    lm(as.formula(formula_str), data = analysis_data_outcome)
  })
  
  model_qc_for_outcome <- purrr::imap_dfr(models_for_outcome, function(model, model_name) {
    model_matrix <- model.matrix(model)
    coefficient_table <- summary(model)$coefficients
    aliased_coefficients <- is.na(coef(model))
    
    tibble(
      outcome_name = outcome_name, outcome_variable = outcome_var, outcome_label = outcome_label,
      model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
      formula = paste(deparse(formula(model)), collapse = " "), n_obs = nobs(model),
      n_parameters = model$rank, residual_df = df.residual(model),
      model_matrix_columns = ncol(model_matrix), model_matrix_rank = qr(model_matrix)$rank,
      full_rank = qr(model_matrix)$rank == ncol(model_matrix),
      n_aliased_coefficients = sum(aliased_coefficients),
      coefficients_finite = all(is.finite(coef(model)[!aliased_coefficients])),
      standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
      residuals_finite = all(is.finite(residuals(model))),
      fitted_values_finite = all(is.finite(fitted(model))),
      residual_standard_error = sigma(model), r_squared = summary(model)$r.squared,
      adjusted_r_squared = summary(model)$adj.r.squared,
      log_likelihood = as.numeric(logLik(model)), AIC = AIC(model)
    )
  })
  
  same_rows_for_outcome <- all(purrr::map_lgl(
    models_for_outcome,
    ~ identical(rownames(model.frame(.x)), rownames(model.frame(models_for_outcome[[1]])))
  ))
  
  expected_sample_n <- continuous_sample_qc$n_analyzed[
    continuous_sample_qc$outcome_name == outcome_name
  ]
  
  model_checks_for_outcome <- list(
    all_models_fitted = length(models_for_outcome) == length(exposure_specs),
    all_models_use_expected_n = all(model_qc_for_outcome$n_obs == expected_sample_n),
    all_models_use_common_n = all(model_qc_for_outcome$n_obs == nrow(analysis_data_outcome)),
    all_models_use_same_rows = same_rows_for_outcome,
    all_model_matrices_full_rank = all(model_qc_for_outcome$full_rank),
    no_aliased_coefficients = all(model_qc_for_outcome$n_aliased_coefficients == 0L),
    all_coefficients_finite = all(model_qc_for_outcome$coefficients_finite),
    all_standard_errors_finite = all(model_qc_for_outcome$standard_errors_finite),
    all_residuals_finite = all(model_qc_for_outcome$residuals_finite),
    all_fitted_values_finite = all(model_qc_for_outcome$fitted_values_finite),
    all_models_have_positive_residual_df = all(model_qc_for_outcome$residual_df > 0L),
    all_aic_values_finite = all(is.finite(model_qc_for_outcome$AIC))
  )
  
  
  stopifnot(isTRUE(all(unlist(model_checks_for_outcome))))
  
  hc3_for_outcome <- purrr::imap_dfr(models_for_outcome, function(model, model_name) {
    tidy_lm_hc3(model) %>%
      mutate(
        outcome_name = outcome_name, outcome_variable = outcome_var, outcome_label = outcome_label,
        model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
        n_obs = nobs(model), residual_df = df.residual(model), .before = 1
      )
  })
  
  tables_for_outcome <- purrr::imap(models_for_outcome, function(model, model_name) {
    regression_table <- tbl_regression(model, tidy_fun = tidy_lm_hc3)

    regression_table <- modify_table_body(regression_table, function(table_body) {
      table_body %>%
        mutate(
          replacement_label = unname(continuous_regression_labels[variable]),
          label = if_else(row_type == "label" & !is.na(replacement_label), replacement_label, label)
        ) %>%
        select(-replacement_label)
    })

    regression_table <- bold_p(regression_table)
    regression_table
  })
  
  aic_for_outcome <- model_qc_for_outcome %>%
    select(
      outcome_name, outcome_variable, outcome_label, model_spec, exposure_label, n_obs,
      n_parameters, log_likelihood, AIC
    )
  
  continuous_analysis_data[[outcome_name]] <- analysis_data_outcome
  continuous_models[[outcome_name]] <- models_for_outcome
  continuous_tables[[outcome_name]] <- tables_for_outcome
  continuous_model_qc[[outcome_name]] <- model_qc_for_outcome
  continuous_aic_qc[[outcome_name]] <- aic_for_outcome
  continuous_hc3_qc[[outcome_name]] <- hc3_for_outcome
  
  
  output_path <- file.path(table_dir, paste0("regression_", outcome_name, "_", site_name, ".html"))
  
  continuous_table_gt <- tbl_merge(tables_for_outcome, tab_spanner = exposure_labels, quiet = TRUE) %>%
    modify_caption(paste0("**Adjusted linear regression models for ", outcome_label, " — ", site_name, "**")) %>%
    as_gt() %>%
    tab_source_note(
      source_note = paste0(
        "Coefficients are unstandardized linear-regression estimates. Confidence intervals and P values use HC3 ",
        "heteroskedasticity-robust standard errors. All exposure specifications use the same complete-case sample: N = ",
        nrow(analysis_data_outcome), "."
      )
    )
  
  gt::gtsave(continuous_table_gt, output_path)

  stopifnot(file.exists(output_path))
  message("Script 04 | continuous output saved: ", outcome_name)
}

# =============================================================================
# Model QC: combined continuous-model results
# =============================================================================

message("Script 04 | continuous-model loop returned; validating combined results")
continuous_model_qc_combined <- bind_rows(continuous_model_qc)
continuous_aic_qc_combined <- bind_rows(continuous_aic_qc)
continuous_hc3_qc_combined <- bind_rows(continuous_hc3_qc)

continuous_output_qc <- purrr::imap_dfr(continuous_outcomes, function(outcome, outcome_name) {
  output_path <- file.path(table_dir, paste0("regression_", outcome_name, "_", site_name, ".html"))
  output_exists <- file.exists(output_path)
  output_size_bytes <- if (output_exists) file.info(output_path, extra_cols = FALSE)$size[[1]] else NA_real_
  
  tibble(
    outcome_name = outcome_name, outcome_variable = outcome$var, outcome_label = outcome$label,
    output_path = output_path, output_exists = output_exists, output_size_bytes = output_size_bytes,
    output_nonempty = output_exists && is.finite(output_size_bytes) && output_size_bytes > 0
  )
})

continuous_basic_model_checks <- list(
  continuous_all_expected_outcomes_fitted = identical(names(continuous_models), names(continuous_outcomes)),
  continuous_expected_models_fitted_per_outcome = all(
    purrr::map_int(continuous_models, length) == length(exposure_specs)
  ),
  continuous_models_use_common_outcome_specific_samples = continuous_model_qc_combined %>%
    group_by(outcome_name) %>%
    summarise(one_sample_size = n_distinct(n_obs) == 1L, .groups = "drop") %>%
    pull(one_sample_size) %>%
    all(),
  continuous_model_matrices_full_rank = all(continuous_model_qc_combined$full_rank),
  continuous_models_have_no_aliased_coefficients = all(
    continuous_model_qc_combined$n_aliased_coefficients == 0L
  ),
  continuous_model_coefficients_finite = all(continuous_model_qc_combined$coefficients_finite),
  continuous_model_standard_errors_finite = all(continuous_model_qc_combined$standard_errors_finite),
  continuous_model_residuals_finite = all(continuous_model_qc_combined$residuals_finite),
  continuous_model_fitted_values_finite = all(continuous_model_qc_combined$fitted_values_finite),
  continuous_model_aic_values_finite = all(is.finite(continuous_model_qc_combined$AIC)),
  continuous_hc3_results_present = nrow(continuous_hc3_qc_combined) > 0L,
  continuous_hc3_estimates_finite = all(is.finite(continuous_hc3_qc_combined$estimate)),
  continuous_hc3_standard_errors_finite = all(is.finite(continuous_hc3_qc_combined$std.error)),
  continuous_hc3_confidence_intervals_finite = all(
    is.finite(continuous_hc3_qc_combined$conf.low) & is.finite(continuous_hc3_qc_combined$conf.high)
  ),
  continuous_hc3_p_values_valid = all(
    continuous_hc3_qc_combined$p.value >= 0 & continuous_hc3_qc_combined$p.value <= 1
  ),
  continuous_model_outputs_exist = all(continuous_output_qc$output_exists),
  continuous_model_outputs_nonempty = all(continuous_output_qc$output_nonempty)
)



stopifnot(isTRUE(all(unlist(continuous_basic_model_checks))))

model_qc$continuous_outcomes$basic_model_diagnostics <- continuous_model_qc_combined
model_qc$continuous_outcomes$HC3_coefficient_results <- continuous_hc3_qc_combined
model_qc$continuous_outcomes$AIC <- continuous_aic_qc_combined
model_qc$continuous_outcomes$output_files <- continuous_output_qc
model_qc$continuous_outcomes$analysis_sample_source_rows <- purrr::imap(
  continuous_analysis_data,
  ~ tibble(outcome_name = .y, source_row = .x$.source_row)
)

model_qc$checks <- c(model_qc$checks, continuous_basic_model_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Basic continuous-model diagnostics added to: ", model_qc_path)

# =============================================================================
# Model QC: continuous-model residuals, specification, and influence
# =============================================================================

if (!requireNamespace("lmtest", quietly = TRUE)) {
  stop("The lmtest package is required for continuous-model diagnostics.")
}

safe_shapiro_test <- function(x, max_n = 5000L) {
  x <- x[is.finite(x)]
  
  if (length(x) < 3L) {
    stop("Shapiro-Wilk testing requires at least three finite observations.")
  }
  
  if (length(x) > max_n) {
    x <- x[seq_len(max_n)]
  }
  
  shapiro.test(x)
}

continuous_residual_qc <- purrr::imap_dfr(continuous_models, function(models_for_outcome, outcome_name) {
  purrr::imap_dfr(models_for_outcome, function(model, model_name) {
    residual_values <- residuals(model)
    standardized_residuals <- rstandard(model)
    studentized_residuals <- rstudent(model)
    cooks_values <- cooks.distance(model)
    leverage_values <- hatvalues(model)
    n <- nobs(model)
    p <- model$rank
    
    residual_z <- as.numeric(scale(residual_values))
    residual_skewness <- mean(residual_z^3)
    residual_excess_kurtosis <- mean(residual_z^4) - 3
    
    shapiro_result <- safe_shapiro_test(residual_values)
    bp_result <- lmtest::bptest(model)
    reset_result <- lmtest::resettest(model, power = 2:3, type = "fitted")
    
    tibble(
      outcome_name = outcome_name, outcome_label = continuous_outcomes[[outcome_name]]$label,
      model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
      n_obs = n, n_parameters = p, residual_skewness = residual_skewness,
      residual_excess_kurtosis = residual_excess_kurtosis,
      shapiro_w = unname(shapiro_result$statistic), shapiro_p_value = shapiro_result$p.value,
      breusch_pagan_statistic = unname(bp_result$statistic),
      breusch_pagan_df = unname(bp_result$parameter), breusch_pagan_p_value = bp_result$p.value,
      reset_statistic = unname(reset_result$statistic),
      reset_df1 = unname(reset_result$parameter[["df1"]]),
      reset_df2 = unname(reset_result$parameter[["df2"]]), reset_p_value = reset_result$p.value,
      maximum_absolute_standardized_residual = max(abs(standardized_residuals)),
      maximum_absolute_studentized_residual = max(abs(studentized_residuals)),
      n_absolute_studentized_residual_gt_3 = sum(abs(studentized_residuals) > 3),
      cooks_threshold = 4 / n, maximum_cooks_distance = max(cooks_values),
      n_cooks_distance_above_threshold = sum(cooks_values > 4 / n),
      leverage_threshold = 2 * p / n, maximum_leverage = max(leverage_values),
      n_leverage_above_threshold = sum(leverage_values > 2 * p / n),
      diagnostics_finite = all(is.finite(c(
        residual_skewness, residual_excess_kurtosis, shapiro_result$statistic, shapiro_result$p.value,
        bp_result$statistic, bp_result$p.value, reset_result$statistic, reset_result$p.value,
        standardized_residuals, studentized_residuals, cooks_values, leverage_values
      )))
    )
  })
})

continuous_influence_qc <- purrr::imap_dfr(continuous_models, function(models_for_outcome, outcome_name) {
  source_rows <- continuous_analysis_data[[outcome_name]]$.source_row
  
  purrr::imap_dfr(models_for_outcome, function(model, model_name) {
    n <- nobs(model)
    p <- model$rank
    cooks_values <- cooks.distance(model)
    leverage_values <- hatvalues(model)
    studentized_residuals <- rstudent(model)
    
    tibble(
      outcome_name = outcome_name, outcome_label = continuous_outcomes[[outcome_name]]$label,
      model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
      analysis_row = seq_len(n), source_row = source_rows,
      cooks_distance = unname(cooks_values), cooks_threshold = 4 / n,
      high_cooks_distance = cooks_distance > cooks_threshold,
      leverage = unname(leverage_values), leverage_threshold = 2 * p / n,
      high_leverage = leverage > leverage_threshold,
      studentized_residual = unname(studentized_residuals),
      large_studentized_residual = abs(studentized_residual) > 3,
      flagged_by_any_measure = high_cooks_distance | high_leverage | large_studentized_residual
    )
  })
})

continuous_influence_summary <- continuous_influence_qc %>%
  group_by(outcome_name, outcome_label, model_spec, exposure_label) %>%
  summarise(
    n_obs = n(), n_high_cooks_distance = sum(high_cooks_distance),
    n_high_leverage = sum(high_leverage),
    n_large_studentized_residual = sum(large_studentized_residual),
    n_flagged_by_any_measure = sum(flagged_by_any_measure),
    maximum_cooks_distance = max(cooks_distance), maximum_leverage = max(leverage),
    maximum_absolute_studentized_residual = max(abs(studentized_residual)), .groups = "drop"
  )

continuous_flagged_observations <- continuous_influence_qc %>%
  filter(flagged_by_any_measure) %>%
  arrange(outcome_name, model_spec, desc(cooks_distance))

continuous_diagnostic_flags <- continuous_residual_qc %>%
  transmute(
    outcome_name, outcome_label, model_spec, exposure_label,
    nonnormal_residual_test = shapiro_p_value < 0.05,
    heteroskedasticity_test = breusch_pagan_p_value < 0.05,
    specification_test = reset_p_value < 0.05,
    large_studentized_residuals_present = n_absolute_studentized_residual_gt_3 > 0L,
    influential_observations_present = n_cooks_distance_above_threshold > 0L,
    high_leverage_observations_present = n_leverage_above_threshold > 0L
  )

continuous_model_sample_reference <- continuous_model_qc_combined %>%
  select(outcome_name, model_spec, expected_n_obs = n_obs)

continuous_diagnostic_sample_check <- continuous_residual_qc %>%
  left_join(continuous_model_sample_reference, by = c("outcome_name", "model_spec")) %>%
  mutate(sample_size_matches = n_obs == expected_n_obs)

continuous_residual_checks <- list(
  expected_number_of_models_diagnosed = nrow(continuous_residual_qc) ==
    length(continuous_outcomes) * length(exposure_specs),
  diagnostic_sample_sizes_match_models = all(continuous_diagnostic_sample_check$sample_size_matches),
  all_diagnostics_finite = all(continuous_residual_qc$diagnostics_finite),
  influence_source_rows_complete = !any(is.na(continuous_influence_qc$source_row)),
  influence_source_rows_valid = all(continuous_influence_qc$source_row %in% seq_len(nrow(cross_sectional))),
  influence_rows_unique_within_model = continuous_influence_qc %>%
    group_by(outcome_name, model_spec) %>%
    summarise(rows_unique = n_distinct(source_row) == n(), .groups = "drop") %>%
    pull(rows_unique) %>%
    all()
)

continuous_diagnostic_summary <- continuous_diagnostic_flags %>%
  group_by(outcome_name, outcome_label) %>%
  summarise(
    n_models = n(), n_nonnormal = sum(nonnormal_residual_test),
    n_heteroskedastic = sum(heteroskedasticity_test),
    n_reset_flagged = sum(specification_test),
    n_with_large_residuals = sum(large_studentized_residuals_present),
    n_with_influential_observations = sum(influential_observations_present),
    n_with_high_leverage = sum(high_leverage_observations_present), .groups = "drop"
  )



stopifnot(isTRUE(all(unlist(continuous_residual_checks))))

model_qc$continuous_outcomes$residual_diagnostics <- continuous_residual_qc
model_qc$continuous_outcomes$diagnostic_flags <- continuous_diagnostic_flags
model_qc$continuous_outcomes$diagnostic_summary <- continuous_diagnostic_summary
model_qc$continuous_outcomes$influence_summary <- continuous_influence_summary
model_qc$continuous_outcomes$flagged_observations <- continuous_flagged_observations
model_qc$continuous_outcomes$diagnostic_thresholds <- tibble(
  diagnostic = c("Cook's distance", "Leverage", "Absolute externally studentized residual"),
  threshold = c("4 / n", "2 × model rank / n", "3")
)

model_qc$checks <- c(model_qc$checks, continuous_residual_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Continuous-model residual and influence diagnostics added to: ", model_qc_path)

# =============================================================================
# Sensitivity analysis: log-transformed elastance and compliance
# =============================================================================

log_sensitivity_outcomes <- list(
  ers = list(variable = "ers", label = "Elastance"),
  crs = list(variable = "crs", label = "Compliance")
)

log_sensitivity_models <- list()
log_sensitivity_results <- list()
log_sensitivity_diagnostics <- list()

for (outcome_name in names(log_sensitivity_outcomes)) {
  outcome_var <- log_sensitivity_outcomes[[outcome_name]]$variable
  outcome_label <- log_sensitivity_outcomes[[outcome_name]]$label
  analysis_data_outcome <- continuous_analysis_data[[outcome_name]]
  covars_used <- model_covariates(outcome_var)
  
  expected_sample_n <- continuous_sample_qc$n_analyzed[
    continuous_sample_qc$outcome_name == outcome_name
  ]
  
  stopifnot(
    all(analysis_data_outcome[[outcome_var]] > 0),
    nrow(analysis_data_outcome) == expected_sample_n
  )
  
  models_for_outcome <- purrr::map(exposure_specs, function(exposure_spec) {
    formula_str <- paste0("log(", outcome_var, ") ~ ", exposure_spec, " + ", covars_used)
    lm(as.formula(formula_str), data = analysis_data_outcome)
  })
  
  results_for_outcome <- purrr::imap_dfr(models_for_outcome, function(model, model_name) {
    tidy_lm_hc3(model) %>%
      mutate(
        outcome_name = outcome_name, outcome_label = outcome_label,
        model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
        n_obs = nobs(model), .before = 1
      ) %>%
      mutate(
        percent_change = 100 * (exp(estimate) - 1),
        percent_change_conf_low = 100 * (exp(conf.low) - 1),
        percent_change_conf_high = 100 * (exp(conf.high) - 1)
      )
  })
  
  diagnostics_for_outcome <- purrr::imap_dfr(models_for_outcome, function(model, model_name) {
    residual_values <- residuals(model)
    residual_z <- as.numeric(scale(residual_values))
    shapiro_result <- safe_shapiro_test(residual_values)
    bp_result <- lmtest::bptest(model)
    reset_result <- lmtest::resettest(model, power = 2:3, type = "fitted")
    
    tibble(
      outcome_name = outcome_name, outcome_label = outcome_label,
      model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
      n_obs = nobs(model), residual_skewness = mean(residual_z^3),
      residual_excess_kurtosis = mean(residual_z^4) - 3,
      shapiro_p_value = shapiro_result$p.value,
      breusch_pagan_p_value = bp_result$p.value,
      reset_p_value = reset_result$p.value,
      maximum_absolute_studentized_residual = max(abs(rstudent(model))),
      maximum_cooks_distance = max(cooks.distance(model)),
      maximum_leverage = max(hatvalues(model))
    )
  })
  
  log_sensitivity_models[[outcome_name]] <- models_for_outcome
  log_sensitivity_results[[outcome_name]] <- results_for_outcome
  log_sensitivity_diagnostics[[outcome_name]] <- diagnostics_for_outcome
}

log_sensitivity_results_combined <- bind_rows(log_sensitivity_results)
log_sensitivity_diagnostics_combined <- bind_rows(log_sensitivity_diagnostics)

expected_log_model_n <- continuous_sample_qc %>%
  filter(outcome_name %in% names(log_sensitivity_outcomes)) %>%
  select(outcome_name, expected_n = n_analyzed)

log_sensitivity_sample_qc <- log_sensitivity_diagnostics_combined %>%
  left_join(expected_log_model_n, by = "outcome_name") %>%
  mutate(sample_size_matches = n_obs == expected_n)

# Because Ers = 1000 / Crs, the log-scale exposure coefficients should be exact sign reversals.
log_reciprocal_check <- log_sensitivity_results_combined %>%
  filter(term != "(Intercept)") %>%
  select(outcome_name, model_spec, term, estimate) %>%
  pivot_wider(names_from = outcome_name, values_from = estimate) %>%
  filter(!is.na(ers), !is.na(crs)) %>%
  mutate(
    coefficient_sum = ers + crs,
    reciprocal_relationship_valid = abs(coefficient_sum) <= 1e-10
  )

log_sensitivity_checks <- list(
  expected_log_models_fitted = sum(purrr::map_int(log_sensitivity_models, length)) ==
    length(log_sensitivity_outcomes) * length(exposure_specs),
  log_models_use_expected_n = all(log_sensitivity_sample_qc$sample_size_matches),
  log_hc3_estimates_finite = all(is.finite(log_sensitivity_results_combined$estimate)),
  log_hc3_standard_errors_finite = all(is.finite(log_sensitivity_results_combined$std.error)),
  log_hc3_confidence_intervals_finite = all(
    is.finite(log_sensitivity_results_combined$conf.low) &
      is.finite(log_sensitivity_results_combined$conf.high)
  ),
  log_hc3_p_values_valid = all(
    log_sensitivity_results_combined$p.value >= 0 &
      log_sensitivity_results_combined$p.value <= 1
  ),
  log_ers_crs_reciprocal_relationship_valid = all(log_reciprocal_check$reciprocal_relationship_valid)
)

log_diagnostic_comparison <- continuous_residual_qc %>%
  filter(outcome_name %in% names(log_sensitivity_outcomes)) %>%
  select(
    outcome_name, model_spec, raw_skewness = residual_skewness,
    raw_excess_kurtosis = residual_excess_kurtosis,
    raw_shapiro_p_value = shapiro_p_value,
    raw_breusch_pagan_p_value = breusch_pagan_p_value,
    raw_reset_p_value = reset_p_value
  ) %>%
  left_join(
    log_sensitivity_diagnostics_combined %>%
      select(
        outcome_name, model_spec, log_skewness = residual_skewness,
        log_excess_kurtosis = residual_excess_kurtosis,
        log_shapiro_p_value = shapiro_p_value,
        log_breusch_pagan_p_value = breusch_pagan_p_value,
        log_reset_p_value = reset_p_value
      ),
    by = c("outcome_name", "model_spec")
  )



stopifnot(isTRUE(all(unlist(log_sensitivity_checks))))

model_qc$continuous_outcomes$log_sensitivity <- list(
  coefficient_results = log_sensitivity_results_combined,
  residual_diagnostics = log_sensitivity_diagnostics_combined,
  sample_size_validation = log_sensitivity_sample_qc,
  raw_vs_log_diagnostic_comparison = log_diagnostic_comparison,
  reciprocal_relationship_check = log_reciprocal_check,
  interpretation = paste0(
    "For log-outcome models, 100 * (exp(beta) - 1) is the estimated percentage change in the outcome per one-unit ",
    "increase in the predictor. Log-elastance and log-compliance coefficients are exact sign reversals because ",
    "ers = 1000 / crs."
  )
)

model_qc$checks <- c(model_qc$checks, log_sensitivity_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Log-scale elastance and compliance sensitivity analyses added to: ", model_qc_path)

# =============================================================================
# 4d2. 28-day VFDs — competing-risks (Fine-Gray) regression
# =============================================================================

# Event of interest = successful final liberation (vfd_status 1); competing event = death within 28 days
# (vfd_status 2); vfd_status 0 represents no documented successful liberation before censoring.

vfd_required_outcome_vars <- c("vfd_time", "vfd_status", "vfd_28")
missing_vfd_outcome_vars <- setdiff(vfd_required_outcome_vars, names(cross_sectional))

stopifnot(length(missing_vfd_outcome_vars) == 0L)

vfd_status_values <- sort(unique(na.omit(cross_sectional$vfd_status)))

vfd_outcome_qc <- cross_sectional %>%
  summarise(
    n_total = n(), n_complete_outcome = sum(complete.cases(vfd_time, vfd_status, vfd_28)),
    n_missing_time = sum(is.na(vfd_time)), n_missing_status = sum(is.na(vfd_status)),
    n_missing_vfd_28 = sum(is.na(vfd_28)), n_status_0 = sum(vfd_status == 0L, na.rm = TRUE),
    n_liberation = sum(vfd_status == 1L, na.rm = TRUE), n_death = sum(vfd_status == 2L, na.rm = TRUE),
    minimum_time = min(vfd_time, na.rm = TRUE), maximum_time = max(vfd_time, na.rm = TRUE),
    minimum_vfd_28 = min(vfd_28, na.rm = TRUE), maximum_vfd_28 = max(vfd_28, na.rm = TRUE)
  )

vfd_outcome_checks <- list(
  required_vfd_outcomes_present = length(missing_vfd_outcome_vars) == 0L,
  vfd_status_codes_valid = all(vfd_status_values %in% c(0L, 1L, 2L)),
  vfd_outcomes_complete = all(complete.cases(cross_sectional[vfd_required_outcome_vars])),
  vfd_times_positive = all(cross_sectional$vfd_time > 0),
  vfd_times_within_28_days = all(cross_sectional$vfd_time <= 28),
  vfd_28_within_range = all(cross_sectional$vfd_28 >= 0 & cross_sectional$vfd_28 <= 28),
  deaths_have_zero_vfd = all(cross_sectional$vfd_28[cross_sectional$vfd_status == 2L] == 0),
  status_0_has_zero_vfd = all(cross_sectional$vfd_28[cross_sectional$vfd_status == 0L] == 0),
  liberation_vfd_matches_time = all(
    abs(
      cross_sectional$vfd_28[cross_sectional$vfd_status == 1L] -
        (28 - cross_sectional$vfd_time[cross_sectional$vfd_status == 1L])
    ) <= 1e-8
  ),
  liberation_events_present = vfd_outcome_qc$n_liberation > 0L,
  competing_deaths_present = vfd_outcome_qc$n_death > 0L
)


stopifnot(isTRUE(all(unlist(vfd_outcome_checks))))

vfd_cr_covariates <- covariates # VFD exposures are not DP-derived, so BMI is not added.

vfd_covariate_vars <- all.vars(as.formula(paste("~", vfd_cr_covariates)))
vfd_analysis_vars <- unique(c(vfd_required_outcome_vars, exposure_vars, vfd_covariate_vars))

vfd_analysis_data <- cross_sectional %>%
  mutate(.source_row = row_number(), .patient_id = row_number()) %>%
  select(.source_row, .patient_id, all_of(vfd_analysis_vars)) %>%
  drop_na(all_of(vfd_analysis_vars))

vfd_sample_qc <- vfd_analysis_data %>%
  summarise(
    n_total_source = nrow(cross_sectional), n_analyzed = n(),
    n_excluded = n_total_source - n_analyzed, pct_analyzed = 100 * n_analyzed / n_total_source,
    n_status_0 = sum(vfd_status == 0L), n_liberation = sum(vfd_status == 1L),
    n_death = sum(vfd_status == 2L)
  )

vfd_sample_by_status_qc <- cross_sectional %>%
  mutate(
    included_in_vfd_analysis = complete.cases(across(all_of(vfd_analysis_vars))),
    vfd_status_label = case_when(
      vfd_status == 0L ~ "No documented successful liberation",
      vfd_status == 1L ~ "Successful final liberation",
      vfd_status == 2L ~ "Death within 28 days",
      TRUE ~ "Missing or unexpected"
    )
  ) %>%
  count(vfd_status, vfd_status_label, included_in_vfd_analysis, name = "n") %>%
  group_by(vfd_status, vfd_status_label) %>%
  mutate(percent_within_status = 100 * n / sum(n)) %>%
  ungroup()

vfd_sample_checks <- list(
  vfd_common_sample_nonempty = nrow(vfd_analysis_data) > 0L,
  vfd_common_sample_complete = all(complete.cases(vfd_analysis_data[vfd_analysis_vars])),
  vfd_common_sample_has_liberations = vfd_sample_qc$n_liberation > 0L,
  vfd_common_sample_has_deaths = vfd_sample_qc$n_death > 0L,
  vfd_patient_ids_unique = !anyDuplicated(vfd_analysis_data$.patient_id),
  vfd_source_rows_unique = !anyDuplicated(vfd_analysis_data$.source_row)
)


# Audit censoring times among patients without documented successful liberation or death.
vfd_censoring_qc <- cross_sectional %>%
  filter(vfd_status == 0L) %>%
  summarise(
    n_censored = n(), n_censored_before_day_28 = sum(vfd_time < 28),
    n_censored_at_day_28 = sum(vfd_time == 28),
    pct_censored_before_day_28 = 100 * mean(vfd_time < 28),
    minimum_censoring_time = min(vfd_time), median_censoring_time = median(vfd_time),
    maximum_censoring_time = max(vfd_time)
  )

vfd_censoring_time_distribution <- cross_sectional %>%
  filter(vfd_status == 0L) %>%
  mutate(
    censoring_category = case_when(
      vfd_time < 7 ~ "Before day 7",
      vfd_time < 14 ~ "Day 7 to <14",
      vfd_time < 21 ~ "Day 14 to <21",
      vfd_time < 28 ~ "Day 21 to <28",
      vfd_time == 28 ~ "Day 28",
      TRUE ~ "Unexpected"
    )
  ) %>%
  count(censoring_category, name = "n") %>%
  mutate(percent = 100 * n / sum(n))

vfd_analysis_censoring_qc <- vfd_analysis_data %>%
  filter(vfd_status == 0L) %>%
  summarise(
    n_censored = n(), n_censored_before_day_28 = sum(vfd_time < 28),
    n_censored_at_day_28 = sum(vfd_time == 28),
    pct_censored_before_day_28 = 100 * mean(vfd_time < 28),
    minimum_censoring_time = min(vfd_time), median_censoring_time = median(vfd_time),
    maximum_censoring_time = max(vfd_time)
  )

vfd_censoring_checks <- list(
  full_cohort_censoring_count_matches = vfd_censoring_qc$n_censored == vfd_outcome_qc$n_status_0,
  analysis_censoring_count_matches = vfd_analysis_censoring_qc$n_censored == vfd_sample_qc$n_status_0,
  full_cohort_censoring_times_valid = all(
    cross_sectional$vfd_time[cross_sectional$vfd_status == 0L] > 0 &
      cross_sectional$vfd_time[cross_sectional$vfd_status == 0L] <= 28
  ),
  analysis_censoring_times_valid = all(
    vfd_analysis_data$vfd_time[vfd_analysis_data$vfd_status == 0L] > 0 &
      vfd_analysis_data$vfd_time[vfd_analysis_data$vfd_status == 0L] <= 28
  )
)


stopifnot(
  isTRUE(all(unlist(vfd_outcome_checks))),
  isTRUE(all(unlist(vfd_sample_checks))),
  isTRUE(all(unlist(vfd_censoring_checks)))
)

model_qc$vfd_competing_risks <- list(
  outcome_summary = vfd_outcome_qc,
  observed_status_values = vfd_status_values,
  required_variables = vfd_analysis_vars,
  analysis_sample = vfd_sample_qc,
  inclusion_by_status = vfd_sample_by_status_qc,
  censoring = list(
    full_cohort_summary = vfd_censoring_qc,
    full_cohort_time_distribution = vfd_censoring_time_distribution,
    analysis_sample_summary = vfd_analysis_censoring_qc
  )
)

model_qc$checks <- c(model_qc$checks, vfd_outcome_checks, vfd_sample_checks, vfd_censoring_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("VFD competing-risk outcome and sample QC added to: ", model_qc_path)

fit_vfd_finegray <- function(exposure_spec) {
  model_rhs <- paste(exposure_spec, "+", vfd_cr_covariates)
  
  finegray_data <- vfd_analysis_data %>%
    mutate(vfd_status_f = factor(
      vfd_status, levels = c(0L, 1L, 2L),
      labels = c("censored", "extubation", "death")
    ))
  
  finegray_formula <- as.formula(
    paste("survival::Surv(vfd_time, vfd_status_f) ~ .patient_id +", model_rhs)
  )
  
  finegray_expanded <- survival::finegray(
    finegray_formula, data = finegray_data, etype = "extubation", id = .patient_id
  )
  
  model_formula <- as.formula(
    paste0("survival::Surv(fgstart, fgstop, fgstatus) ~ ", model_rhs, " + cluster(.patient_id)")
  )
  
  model <- survival::coxph(
    model_formula, data = finegray_expanded, weights = fgwt,
    robust = TRUE, x = TRUE, model = TRUE
  )
  
  attr(model, "n_patients") <- nrow(finegray_data)
  attr(model, "n_finegray_rows") <- nrow(finegray_expanded)
  attr(model, "n_liberation") <- sum(finegray_data$vfd_status == 1L)
  attr(model, "n_death") <- sum(finegray_data$vfd_status == 2L)
  attr(model, "n_censored") <- sum(finegray_data$vfd_status == 0L)
  
  model
}

vfd_cr_models <- purrr::map(exposure_specs, fit_vfd_finegray)

vfd_model_qc <- purrr::imap_dfr(vfd_cr_models, function(model, model_name) {
  coefficient_table <- summary(model)$coefficients
  model_matrix <- model$x
  
  tibble(
    model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
    formula = paste(deparse(formula(model)), collapse = " "),
    n_patients = attr(model, "n_patients"), n_finegray_rows = attr(model, "n_finegray_rows"),
    n_liberation = attr(model, "n_liberation"), n_death = attr(model, "n_death"),
    n_censored = attr(model, "n_censored"), n_coefficients = length(coef(model)),
    model_matrix_columns = ncol(model_matrix), model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    coefficients_finite = all(is.finite(coef(model))),
    robust_standard_errors_finite = all(is.finite(coefficient_table[, "robust se"])),
    log_likelihood = as.numeric(logLik(model)), AIC = AIC(model)
  )
})

vfd_same_patient_counts <- all(vfd_model_qc$n_patients == nrow(vfd_analysis_data))
vfd_same_expanded_row_counts <- n_distinct(vfd_model_qc$n_finegray_rows) == 1L

vfd_model_checks <- list(
  vfd_all_models_fitted = length(vfd_cr_models) == length(exposure_specs),
  vfd_all_models_use_common_patient_n = vfd_same_patient_counts,
  vfd_all_models_use_common_expanded_n = vfd_same_expanded_row_counts,
  vfd_all_models_have_expected_liberations = all(vfd_model_qc$n_liberation == vfd_sample_qc$n_liberation),
  vfd_all_models_have_expected_deaths = all(vfd_model_qc$n_death == vfd_sample_qc$n_death),
  vfd_all_models_have_expected_censoring = all(vfd_model_qc$n_censored == vfd_sample_qc$n_status_0),
  vfd_model_matrices_full_rank = all(vfd_model_qc$full_rank),
  vfd_coefficients_finite = all(vfd_model_qc$coefficients_finite),
  vfd_robust_standard_errors_finite = all(vfd_model_qc$robust_standard_errors_finite),
  vfd_aic_values_finite = all(is.finite(vfd_model_qc$AIC))
)


stopifnot(isTRUE(all(unlist(vfd_model_checks))))

vfd_regression_labels <- c(
  vtpbw = "VT/McLaren PBW (mL/kg)",
  vtnchs = "VT/NCHS weight (mL/kg)",
  vtabw = "VT/ABW (mL/kg)",
  vtpfvc = "VT/PFVC (%)",
  pfvc = "Predicted FVC (L)",
  pbwpfvc = "McLaren PBW/PFVC (kg/L)",
  nchspfvc = "NCHS weight/PFVC (kg/L)",
  abwpfvc = "ABW/PFVC (kg/L)",
  race_category = "Race",
  age_at_admission = "Age at admission (per year)",
  sex_category = "Sex",
  sofa_total = "pSOFA total",
  sf10 = "SF ratio (per 10)"
)

vfd_cr_tables <- purrr::map(vfd_cr_models, function(model) {
  tbl_regression(model, exponentiate = TRUE) %>%
    modify_table_body(function(table_body) {
      table_body %>%
        mutate(
          replacement_label = unname(vfd_regression_labels[variable]),
          label = if_else(row_type == "label" & !is.na(replacement_label), replacement_label, label)
        ) %>%
        select(-replacement_label)
    }) %>%
    modify_header(estimate ~ "**SHR**") %>%
    bold_p()
})


vfd_table <- tbl_merge(vfd_cr_tables, tab_spanner = exposure_labels, quiet = TRUE) %>%
  modify_caption(paste0("**Fine–Gray regression models for successful liberation by day 28 — ", site_name, "**"))

vfd_table_gt <- vfd_table %>%
  as_gt() %>%
  tab_source_note(
    source_note = paste0(
      "Subdistribution hazard ratios greater than 1 indicate a higher cumulative incidence of successful liberation. ",
      "Death within 28 days is the competing event. All models use the same complete-case sample: N = ",
      nrow(vfd_analysis_data), ", including ", vfd_sample_qc$n_liberation, " liberation events, ",
      vfd_sample_qc$n_death, " deaths, and ", vfd_sample_qc$n_status_0,
      " patients without documented successful liberation."
    )
  ) %>%
  tab_source_note(
    source_note = paste0(
      "Robust standard errors account for repeated rows created by the Fine–Gray expansion. ",
      vfd_analysis_censoring_qc$n_censored_before_day_28,
      " patients without documented liberation were censored before day 28."
    )
  )

vfd_output_path <- file.path(table_dir, paste0("regression_vfd28_", site_name, ".html"))

gt::gtsave(vfd_table_gt, vfd_output_path)

vfd_output_exists <- file.exists(vfd_output_path)
vfd_output_size_bytes <- if (vfd_output_exists) {
  file.info(vfd_output_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

vfd_output_checks <- list(
  vfd_output_exists = vfd_output_exists,
  vfd_output_nonempty = vfd_output_exists &&
    is.finite(vfd_output_size_bytes) &&
    vfd_output_size_bytes > 0
)

model_qc$vfd_competing_risks$model_diagnostics <- vfd_model_qc
model_qc$vfd_competing_risks$output <- list(
  path = vfd_output_path, exists = vfd_output_exists, size_bytes = vfd_output_size_bytes,
  common_sample_n = nrow(vfd_analysis_data), liberation_n = vfd_sample_qc$n_liberation,
  death_n = vfd_sample_qc$n_death, censored_n = vfd_sample_qc$n_status_0
)

model_qc$checks <- c(model_qc$checks, vfd_model_checks, vfd_output_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Fine-Gray regression table saved to: ", vfd_output_path)

# =============================================================================
# Model QC: proportional subdistribution-hazards diagnostics
# =============================================================================

vfd_zph_objects <- purrr::map(vfd_cr_models, function(model) {
  survival::cox.zph(model, transform = "km", terms = TRUE, singledf = FALSE)
})

vfd_zph_qc <- purrr::imap_dfr(vfd_zph_objects, function(zph_object, model_name) {
  zph_table <- as.data.frame(zph_object$table) %>%
    rownames_to_column("term") %>%
    as_tibble()
  
  names(zph_table) <- c("term", "chisq", "df", "p_value")
  
  zph_table %>%
    mutate(
      model_spec = model_name,
      exposure_label = unname(exposure_labels[[model_name]]),
      .before = 1
    )
})

vfd_zph_summary <- vfd_zph_qc %>%
  mutate(
    is_global_test = term == "GLOBAL",
    proportional_hazards_flag = p_value < 0.05
  ) %>%
  group_by(model_spec, exposure_label) %>%
  summarise(
    n_term_tests = sum(!is_global_test),
    n_term_flags = sum(proportional_hazards_flag & !is_global_test),
    global_p_value = p_value[is_global_test][1],
    global_test_flag = global_p_value < 0.05,
    .groups = "drop"
  )

vfd_zph_checks <- list(
  vfd_zph_all_models_evaluated = length(vfd_zph_objects) == length(vfd_cr_models),
  vfd_zph_all_models_have_global_test = vfd_zph_qc %>%
    group_by(model_spec) %>%
    summarise(has_global = any(term == "GLOBAL"), .groups = "drop") %>%
    pull(has_global) %>%
    all(),
  vfd_zph_statistics_finite = all(is.finite(vfd_zph_qc$chisq)),
  vfd_zph_degrees_of_freedom_valid = all(is.finite(vfd_zph_qc$df) & vfd_zph_qc$df > 0),
  vfd_zph_p_values_valid = all(vfd_zph_qc$p_value >= 0 & vfd_zph_qc$p_value <= 1)
)


stopifnot(isTRUE(all(unlist(vfd_zph_checks))))

model_qc$vfd_competing_risks$proportional_subdistribution_hazards <- list(
  tests = vfd_zph_qc,
  summary = vfd_zph_summary,
  interpretation = paste0(
    "These are proportional-hazards diagnostics applied to the weighted Cox representation of the Fine-Gray model. ",
    "A P value below 0.05 flags possible time variation in a subdistribution hazard ratio; it is not treated as a ",
    "computational model failure."
  )
)

model_qc$checks <- c(model_qc$checks, vfd_zph_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Fine-Gray proportional-hazards diagnostics added to: ", model_qc_path)

# =============================================================================
# 4e. AIC comparison across all models and outcomes
# =============================================================================

aic_results <- list()

# Mortality
if (!is.null(mortality_models)) {
  aic_results[["Mortality"]] <- tibble(
    exposure = exposure_labels,
    AIC = map_dbl(mortality_models, AIC),
    is_reference = exposure == "VT/McLaren PBW"
  )
}

# Continuous outcomes. Exclude normalized-mechanics outcomes from the AIC comparison and heatmap because the outcome
# shares a predicted-size variable with size-containing exposures. Their standalone regression tables and long-format
# results are still produced.
AIC_EXCLUDE <- c("ers_pbw", "ers_pfvc", "mp_pbw", "mp_pfvc")

for (outcome_name in setdiff(names(continuous_outcomes), AIC_EXCLUDE)) {
  aic_results[[continuous_outcomes[[outcome_name]]$label]] <- tibble(
    exposure = exposure_labels,
    AIC = map_dbl(continuous_models[[outcome_name]], AIC),
    is_reference = exposure == "VT/McLaren PBW"
  )
}

# AICs for 28-day VFD Fine-Gray models are compared within the outcome and referenced to VT/McLaren PBW.
aic_results[["28-day VFDs"]] <- tibble(
  exposure = exposure_labels,
  AIC = map_dbl(vfd_cr_models, AIC),
  is_reference = exposure == "VT/McLaren PBW"
)

# Validate the model sets and analysis samples used for within-outcome AIC comparisons.
aic_input_qc <- list()

if (!is.null(mortality_models)) {
  aic_input_qc[["Mortality"]] <- tibble(
    outcome = "Mortality",
    model_spec = names(mortality_models),
    exposure = unname(exposure_labels[names(mortality_models)]),
    n_patients = purrr::map_int(mortality_models, nobs),
    n_likelihood_rows = purrr::map_int(mortality_models, nobs),
    AIC = purrr::map_dbl(mortality_models, AIC),
    likelihood_type = "Binomial likelihood"
  )
}

for (outcome_name in setdiff(names(continuous_outcomes), AIC_EXCLUDE)) {
  models_for_outcome <- continuous_models[[outcome_name]]
  outcome_label <- continuous_outcomes[[outcome_name]]$label
  
  aic_input_qc[[outcome_label]] <- tibble(
    outcome = outcome_label,
    model_spec = names(models_for_outcome),
    exposure = unname(exposure_labels[names(models_for_outcome)]),
    n_patients = purrr::map_int(models_for_outcome, nobs),
    n_likelihood_rows = purrr::map_int(models_for_outcome, nobs),
    AIC = purrr::map_dbl(models_for_outcome, AIC),
    likelihood_type = "Gaussian likelihood"
  )
}

aic_input_qc[["28-day VFDs"]] <- tibble(
  outcome = "28-day VFDs",
  model_spec = names(vfd_cr_models),
  exposure = unname(exposure_labels[names(vfd_cr_models)]),
  n_patients = purrr::map_int(vfd_cr_models, ~ attr(.x, "n_patients")),
  n_likelihood_rows = purrr::map_int(vfd_cr_models, ~ attr(.x, "n_finegray_rows")),
  AIC = purrr::map_dbl(vfd_cr_models, AIC),
  likelihood_type = "Fine-Gray weighted Cox partial likelihood"
)

aic_input_qc <- bind_rows(aic_input_qc)

aic_outcome_qc <- aic_input_qc %>%
  group_by(outcome) %>%
  summarise(
    n_models = n(),
    n_distinct_patient_counts = n_distinct(n_patients),
    n_distinct_likelihood_row_counts = n_distinct(n_likelihood_rows),
    all_aic_finite = all(is.finite(AIC)),
    n_reference_models = sum(exposure == "VT/McLaren PBW"),
    minimum_AIC = min(AIC),
    maximum_AIC = max(AIC),
    AIC_range = maximum_AIC - minimum_AIC,
    .groups = "drop"
  )

aic_input_checks <- list(
  expected_outcomes_present = identical(sort(names(aic_results)), sort(unique(aic_input_qc$outcome))),
  expected_candidate_models_per_outcome = all(aic_outcome_qc$n_models == length(exposure_specs)),
  common_patient_count_within_outcome = all(aic_outcome_qc$n_distinct_patient_counts == 1L),
  common_likelihood_rows_within_outcome = all(aic_outcome_qc$n_distinct_likelihood_row_counts == 1L),
  all_aic_values_finite = all(aic_outcome_qc$all_aic_finite),
  exactly_one_reference_model_per_outcome = all(aic_outcome_qc$n_reference_models == 1L)
)


stopifnot(isTRUE(all(unlist(aic_input_checks))))

model_qc$aic_comparison <- list(
  input_models = aic_input_qc,
  outcome_summary = aic_outcome_qc,
  excluded_outcomes = AIC_EXCLUDE,
  reference_model = "VT/McLaren PBW",
  interpretation = paste0(
    "AIC comparisons are made only within outcomes. Mortality uses a binomial likelihood, continuous outcomes use ",
    "Gaussian likelihoods, and 28-day VFD models use the weighted Cox partial-likelihood representation of Fine-Gray ",
    "regression. Fine-Gray AIC comparisons should be interpreted as exploratory."
  )
)

model_qc$checks <- c(model_qc$checks, aic_input_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("AIC-comparison input QC added to: ", model_qc_path)

# Evidence ratios are referenced to the VT/McLaren PBW alone model within each outcome:
# ER = exp(-0.5 * (AIC_model - AIC_VT/McLaren PBW)). ER > 1 indicates more support than
# VT/McLaren PBW alone, and ER = 1 for the VT/McLaren PBW model.
# Values displayed in the heatmap are truncated to [0.001, 1000].
ER_FLOOR <- 0.001
ER_CEIL <- 1000

aic_all <- bind_rows(aic_results, .id = "outcome") %>%
  group_by(outcome) %>%
  mutate(
    aic_ref = AIC[is_reference][1],
    delta_AIC = AIC - aic_ref,
    log_evidence_ratio = -0.5 * delta_AIC,
    evidence_ratio = exp(log_evidence_ratio),
    evidence_ratio_trunc = exp(pmin(pmax(log_evidence_ratio, log(ER_FLOOR)), log(ER_CEIL))),
    er_label = case_when(
      log_evidence_ratio > log(ER_CEIL) ~ ">1000",
      log_evidence_ratio < log(ER_FLOOR) ~ "<0.001",
      TRUE ~ trimws(formatC(evidence_ratio, format = "g", digits = 2))
    )
  ) %>%
  ungroup() %>%
  mutate(
    outcome = factor(
      outcome,
      levels = intersect(
        c("Mortality", "28-day VFDs", "Compliance", "Elastance", "Driving pressure", "Mechanical power"),
        unique(outcome)
      )
    )
  )

aic_evidence_ratio_qc <- aic_all %>%
  group_by(outcome) %>%
  summarise(
    n_models = n(),
    n_reference_models = sum(is_reference),
    reference_delta_aic = delta_AIC[is_reference][1],
    reference_evidence_ratio = evidence_ratio[is_reference][1],
    all_delta_aic_finite = all(is.finite(delta_AIC)),
    all_log_evidence_ratios_finite = all(is.finite(log_evidence_ratio)),
    all_evidence_ratios_finite = all(is.finite(evidence_ratio)),
    all_truncated_ratios_within_limits = all(
      evidence_ratio_trunc >= ER_FLOOR & evidence_ratio_trunc <= ER_CEIL
    ),
    .groups = "drop"
  )

aic_evidence_ratio_checks <- list(
  all_outcome_labels_retained = !any(is.na(aic_all$outcome)),
  expected_number_of_rows = nrow(aic_all) == length(aic_results) * length(exposure_specs),
  expected_models_per_outcome = all(aic_evidence_ratio_qc$n_models == length(exposure_specs)),
  exactly_one_reference_per_outcome = all(aic_evidence_ratio_qc$n_reference_models == 1L),
  reference_delta_aic_is_zero = all(abs(aic_evidence_ratio_qc$reference_delta_aic) <= 1e-12),
  reference_evidence_ratio_is_one = all(abs(aic_evidence_ratio_qc$reference_evidence_ratio - 1) <= 1e-12),
  delta_aic_values_finite = all(aic_evidence_ratio_qc$all_delta_aic_finite),
  log_evidence_ratios_finite = all(aic_evidence_ratio_qc$all_log_evidence_ratios_finite),
  evidence_ratios_finite = all(aic_evidence_ratio_qc$all_evidence_ratios_finite),
  truncated_evidence_ratios_within_limits = all(aic_evidence_ratio_qc$all_truncated_ratios_within_limits)
)


stopifnot(isTRUE(all(unlist(aic_evidence_ratio_checks))))

model_qc$aic_comparison$evidence_ratios <- aic_all
model_qc$aic_comparison$evidence_ratio_qc <- aic_evidence_ratio_qc
model_qc$aic_comparison$cautions <- c(
  "Evidence ratios describe relative AIC support within an outcome and are not causal effect measures.",
  "Fine-Gray AIC uses a weighted Cox partial-likelihood representation and is exploratory.",
  "Gaussian AIC results for compliance and elastance should be interpreted with their residual diagnostics."
)

model_qc$checks <- c(model_qc$checks, aic_evidence_ratio_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("AIC evidence-ratio QC added to: ", model_qc_path)

# Use a prespecified exposure order so heatmap rows remain consistent across sites. ggplot displays the first factor
# level at the bottom of a discrete y-axis, so the intended top-to-bottom order is reversed in the factor levels.
exposure_order_top_to_bottom <- unname(exposure_labels)

stopifnot(setequal(unique(as.character(aic_all$exposure)), exposure_order_top_to_bottom))

aic_all <- aic_all %>%
  mutate(exposure = factor(as.character(exposure), levels = rev(exposure_order_top_to_bottom)))

message("AIC comparison (evidence ratios vs VT/McLaren PBW alone within each outcome):")

# Divergent log10 color scale: white = 1, blue = less support, and red = more support than VT/PBW alone.
er_heatmap <- ggplot(aic_all, aes(x = outcome, y = exposure, fill = log10(evidence_ratio_trunc))) +
  geom_tile(color = "grey80", linewidth = 0.5) +
  geom_text(aes(label = er_label), size = 3.5) +
  scale_fill_gradient2(
    name = "Evidence ratio\n(vs VT/McLaren PBW)",
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(log10(ER_FLOOR), log10(ER_CEIL)),
    breaks = -3:3, labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")
  ) +
  labs(
    title = "Relative AIC support across models and outcomes",
    subtitle = stringr::str_wrap(
      paste0(
        "Each cell compares a model with VT/McLaren PBW alone within the same outcome; ",
        "values are truncated to [0.001, 1000]. ",
        "Fine–Gray and diagnostically flagged Gaussian comparisons are exploratory."
      ),
      width = 95
    ),
    x = "Outcome", y = "Exposure specification",
    caption = stringr::str_wrap(
      paste0(
        "Values above 1 indicate stronger relative model support than VT/McLaren PBW alone; ",
        "values below 1 indicate weaker support. A value of 1 indicates equal support. ",
        "Values compare models within each outcome and do not represent patient risk or causal effects."
      ),
      width = 110
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9, margin = margin(t = 12))
  )

aic_output_path <- file.path(table_dir, paste0("aic_comparison_all_", site_name, ".csv"))
er_heatmap_output_path <- file.path(figure_dir,paste0("evidence_ratio_heatmap_all_", site_name, ".pdf"))

write_csv(aic_all, aic_output_path)
save_script04_plot(er_heatmap_output_path, er_heatmap, width = 10, height = 13)

aic_output_checks <- list(
  aic_csv_exists = file.exists(aic_output_path),
  aic_heatmap_exists = file.exists(er_heatmap_output_path),
  aic_csv_nonempty = file.exists(aic_output_path) &&
    file.info(aic_output_path, extra_cols = FALSE)$size[[1]] > 0,
  aic_heatmap_nonempty = file.exists(er_heatmap_output_path) &&
    file.info(er_heatmap_output_path, extra_cols = FALSE)$size[[1]] > 0
)

model_qc$aic_comparison$output <- list(
  csv_path = aic_output_path,
  heatmap_path = er_heatmap_output_path,
  exposure_order_top_to_bottom = exposure_order_top_to_bottom,
  outcome_order_left_to_right = levels(aic_all$outcome)
)

model_qc$checks <- c(model_qc$checks, aic_output_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("AIC comparison and evidence-ratio heatmap saved.")

# =============================================================================
# 4f. Survival analysis
# =============================================================================

# Initialize survival objects so an interactive rerun cannot reuse models or plots created from an earlier dataset.
cox_model <- NULL
cox_model_nchspfvc <- NULL
cox_model_abwpfvc <- NULL
cox_model_pfvc <- NULL

cox_model_race_stratified <- NULL
cox_model_nchspfvc_race_stratified <- NULL
cox_model_abwpfvc_race_stratified <- NULL
cox_model_pfvc_race_stratified <- NULL

cox_model_spline <- NULL
cox_model_nchspfvc_spline <- NULL
cox_model_abwpfvc_spline <- NULL
cox_model_pfvc_spline <- NULL

cox_unadj <- NULL
cox_nchspfvc_unadj <- NULL
cox_abwpfvc_unadj <- NULL
cox_pfvc_unadj <- NULL

km_fit <- NULL
km_plot <- NULL

# Event = all-cause death within 60 days, in- or out-of-hospital, derived in Script 03 from the patient-level death
# timestamp. Survivors are censored administratively at day 60 rather than at hospital discharge.
survival_required_vars <- c("surv_time", "mortality_event_60")
missing_survival_vars <- setdiff(survival_required_vars, names(cross_sectional))

stopifnot(length(missing_survival_vars) == 0L)

survival_event_values <- sort(unique(na.omit(cross_sectional$mortality_event_60)))

survival_endpoint_qc <- cross_sectional %>%
  summarise(
    n_total = n(), n_complete = sum(complete.cases(surv_time, mortality_event_60)),
    n_missing_time = sum(is.na(surv_time)), n_missing_event = sum(is.na(mortality_event_60)),
    n_events = sum(mortality_event_60 == 1L, na.rm = TRUE),
    n_censored = sum(mortality_event_60 == 0L, na.rm = TRUE),
    n_nonfinite_times = sum(!is.na(surv_time) & !is.finite(surv_time)),
    n_nonpositive_times = sum(surv_time <= 0, na.rm = TRUE),
    n_times_above_60 = sum(surv_time > 60, na.rm = TRUE),
    n_events_after_day_60 = sum(mortality_event_60 == 1L & surv_time > 60, na.rm = TRUE),
    n_survivors_not_censored_at_60 = sum(mortality_event_60 == 0L & surv_time != 60, na.rm = TRUE),
    minimum_time = min(surv_time, na.rm = TRUE), median_time = median(surv_time, na.rm = TRUE),
    maximum_time = max(surv_time, na.rm = TRUE)
  )

survival_time_by_status_qc <- cross_sectional %>%
  mutate(
    survival_status = case_when(
      mortality_event_60 == 0L ~ "Censored alive at day 60",
      mortality_event_60 == 1L ~ "Death by day 60",
      TRUE ~ "Missing or unexpected"
    )
  ) %>%
  group_by(mortality_event_60, survival_status) %>%
  summarise(
    n = n(), minimum_time = min(surv_time, na.rm = TRUE),
    median_time = median(surv_time, na.rm = TRUE), maximum_time = max(surv_time, na.rm = TRUE),
    .groups = "drop"
  )

survival_endpoint_checks <- list(
  required_survival_variables_present = length(missing_survival_vars) == 0L,
  survival_outcomes_complete = survival_endpoint_qc$n_complete == nrow(cross_sectional),
  survival_event_codes_valid = identical(survival_event_values, c(0L, 1L)),
  survival_times_finite = survival_endpoint_qc$n_nonfinite_times == 0L,
  survival_times_positive = survival_endpoint_qc$n_nonpositive_times == 0L,
  survival_times_within_60_days = survival_endpoint_qc$n_times_above_60 == 0L,
  all_events_occur_by_day_60 = survival_endpoint_qc$n_events_after_day_60 == 0L,
  all_survivors_censored_at_day_60 = survival_endpoint_qc$n_survivors_not_censored_at_60 == 0L,
  survival_event_count_matches_mortality = survival_endpoint_qc$n_events ==
    sum(cross_sectional$mortality_event_60 == 1L, na.rm = TRUE)
)


stopifnot(isTRUE(all(unlist(survival_endpoint_checks))))

# Retain the full outcome-valid cohort for Kaplan-Meier analysis. A separate common complete-case sample is created for
# the adjusted Cox models.
surv_data <- cross_sectional %>%
  filter(complete.cases(surv_time, mortality_event_60)) %>%
  mutate(event = as.integer(mortality_event_60))

n_deaths <- sum(surv_data$event)
n_censored <- sum(surv_data$event == 0L)

stopifnot(
  nrow(surv_data) == survival_endpoint_qc$n_complete,
  n_deaths == survival_endpoint_qc$n_events,
  n_censored == survival_endpoint_qc$n_censored
)

model_qc$survival <- list(
  endpoint_summary = survival_endpoint_qc,
  time_by_status = survival_time_by_status_qc,
  observed_event_values = survival_event_values
)

model_qc$checks <- c(model_qc$checks, survival_endpoint_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message(
  "Survival endpoint validated: ", nrow(surv_data), " patients, ", n_deaths,
  " deaths by day 60, and ", n_censored, " censored alive at day 60."
)

message("Survival-endpoint QC added to: ", model_qc_path)

# Use one common complete-case sample for both adjusted Cox models so coefficients, likelihoods, and AIC values are based
# on the same patients.
cox_adjustment_vars <- c(
  "age_at_admission", "sex_category", "race_category", "sf10", "sofa_total"
)

cox_exposure_vars <- c(
  "vtpbw", "vtabw", "vtnchs", "pfvc", "pbwpfvc", "nchspfvc", "abwpfvc"
)

cox_required_vars <- unique(c(
  "surv_time", "mortality_event_60", cox_exposure_vars, cox_adjustment_vars
))

missing_cox_vars <- setdiff(cox_required_vars, names(cross_sectional))

stopifnot(length(missing_cox_vars) == 0L)

cox_analysis_data <- cross_sectional %>%
  mutate(.source_row = row_number(), event = as.integer(mortality_event_60)) %>%
  select(.source_row, event, all_of(cox_required_vars)) %>%
  drop_na(all_of(cox_required_vars))

cox_sample_qc <- cox_analysis_data %>%
  summarise(
    n_total_source = nrow(cross_sectional), n_analyzed = n(),
    n_excluded = n_total_source - n_analyzed, pct_analyzed = 100 * n_analyzed / n_total_source,
    n_events = sum(event == 1L), n_censored = sum(event == 0L),
    event_percent = 100 * mean(event == 1L),
    minimum_time = min(surv_time), median_time = median(surv_time), maximum_time = max(surv_time)
  )

cox_sample_by_status_qc <- cross_sectional %>%
  mutate(
    included_in_cox_analysis = complete.cases(across(all_of(cox_required_vars))),
    survival_status = case_when(
      mortality_event_60 == 0L ~ "Censored alive at day 60",
      mortality_event_60 == 1L ~ "Death by day 60",
      TRUE ~ "Missing or unexpected"
    )
  ) %>%
  count(mortality_event_60, survival_status, included_in_cox_analysis, name = "n") %>%
  group_by(mortality_event_60, survival_status) %>%
  mutate(percent_within_status = 100 * n / sum(n)) %>%
  ungroup()

cox_sample_checks <- list(
  cox_required_variables_present = length(missing_cox_vars) == 0L,
  cox_common_sample_nonempty = nrow(cox_analysis_data) > 0L,
  cox_common_sample_complete = all(complete.cases(cox_analysis_data[cox_required_vars])),
  cox_common_sample_has_events = cox_sample_qc$n_events > 0L,
  cox_common_sample_has_censored_observations = cox_sample_qc$n_censored > 0L,
  cox_common_sample_event_codes_valid = all(cox_analysis_data$event %in% c(0L, 1L)),
  cox_common_sample_times_positive = all(cox_analysis_data$surv_time > 0),
  cox_common_sample_times_within_60_days = all(cox_analysis_data$surv_time <= 60),
  cox_source_rows_complete = !any(is.na(cox_analysis_data$.source_row)),
  cox_source_rows_unique = !anyDuplicated(cox_analysis_data$.source_row)
)


stopifnot(isTRUE(all(unlist(cox_sample_checks))))

model_qc$survival$cox_analysis_sample <- list(
  required_variables = cox_required_vars,
  missing_variables = missing_cox_vars,
  sample_summary = cox_sample_qc,
  inclusion_by_status = cox_sample_by_status_qc,
  source_rows = cox_analysis_data %>% select(.source_row)
)

model_qc$checks <- c(model_qc$checks, cox_sample_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Common Cox-model sample QC added to: ", model_qc_path)

# The preceding checks require both deaths and censored observations, so the adjusted Cox models are fit directly.
cox_model <- coxph(
  Surv(surv_time, event) ~ pbwpfvc + vtpbw + age_at_admission +
    sex_category + race_category + sf10 + sofa_total,
  data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
)

# Companion model with PFVC as the scaling exposure. Both adjusted Cox models use the same complete-case sample.
cox_model_pfvc <- coxph(
  Surv(surv_time, event) ~ pfvc + vtpbw + age_at_admission +
    sex_category + race_category + sf10 + sofa_total,
  data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
)

cox_model_abwpfvc <- coxph(
  Surv(surv_time, event) ~ abwpfvc + vtabw + age_at_admission +
    sex_category + race_category + sf10 + sofa_total,
  data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
)

cox_model_nchspfvc <- coxph(
  Surv(surv_time, event) ~ nchspfvc + vtnchs + age_at_admission +
    sex_category + race_category + sf10 + sofa_total,
  data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
)

cox_models <- list(
  vtpbw_pbwpfvc = cox_model,
  vtnchs_nchspfvc = cox_model_nchspfvc,
  vtabw_abwpfvc = cox_model_abwpfvc,
  vtpbw_pfvc = cox_model_pfvc
)

cox_model_labels <- c(
  vtpbw_pbwpfvc = "VT/McLaren PBW + McLaren PBW/PFVC",
  vtnchs_nchspfvc = "VT/NCHS weight + NCHS weight/PFVC",
  vtabw_abwpfvc = "VT/ABW + ABW/PFVC",
  vtpbw_pfvc = "VT/McLaren PBW + PFVC"
)

cox_model_qc <- purrr::imap_dfr(cox_models, function(model, model_name) {
  coefficient_table <- summary(model)$coefficients
  model_matrix <- model$x
  survival_response <- model$y
  n_patients_model <- model$n
  n_events_model <- model$nevent
  n_coefficients <- length(coef(model))
  
  tibble(
    model_spec = model_name, label = unname(cox_model_labels[[model_name]]),
    formula = paste(deparse(formula(model)), collapse = " "),
    n_obs = n_patients_model, n_events = n_events_model, n_censored = n_patients_model - n_events_model,
    n_response_rows = nrow(survival_response), n_coefficients = n_coefficients,
    events_per_coefficient = n_events_model / n_coefficients, iterations = model$iter,
    model_matrix_columns = ncol(model_matrix), model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    coefficients_finite = all(is.finite(coef(model))),
    standard_errors_finite = all(is.finite(coefficient_table[, "se(coef)"])),
    maximum_absolute_coefficient = max(abs(coef(model))),
    maximum_standard_error = max(coefficient_table[, "se(coef)"]),
    log_likelihood = as.numeric(logLik(model)), AIC = AIC(model)
  )
})

cox_reference_rows <- rownames(model.frame(cox_models[[1L]]))

cox_same_analysis_rows <- all(
  purrr::map_lgl(cox_models, ~ identical(rownames(model.frame(.x)), cox_reference_rows))
)

cox_model_checks <- list(
  all_expected_cox_models_fitted = length(cox_models) == 4L,
  cox_models_use_common_n = all(cox_model_qc$n_obs == nrow(cox_analysis_data)),
  cox_response_rows_match_common_n = all(cox_model_qc$n_response_rows == nrow(cox_analysis_data)),
  cox_models_use_same_rows = cox_same_analysis_rows,
  cox_models_have_expected_events = all(cox_model_qc$n_events == cox_sample_qc$n_events),
  cox_models_have_expected_censoring = all(cox_model_qc$n_censored == cox_sample_qc$n_censored),
  cox_model_matrices_full_rank = all(cox_model_qc$full_rank),
  cox_coefficients_finite = all(cox_model_qc$coefficients_finite),
  cox_standard_errors_finite = all(cox_model_qc$standard_errors_finite),
  cox_log_likelihoods_finite = all(is.finite(cox_model_qc$log_likelihood)),
  cox_aic_values_finite = all(is.finite(cox_model_qc$AIC))
)


stopifnot(isTRUE(all(unlist(cox_model_checks))))

model_qc$survival$cox_model_diagnostics <- cox_model_qc
model_qc$survival$events_per_coefficient_below_10 <- cox_model_qc %>%
  select(model_spec, label, events_per_coefficient) %>%
  mutate(below_10 = events_per_coefficient < 10)

model_qc$checks <- c(model_qc$checks, cox_model_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Adjusted Cox-model diagnostics added to: ", model_qc_path)

# Test whether adjusted hazard ratios remain approximately constant over the 60-day follow-up period.
cox_zph_objects <- purrr::map(cox_models, function(model) {
  survival::cox.zph(model, transform = "km", terms = TRUE, singledf = FALSE)
})

cox_zph_qc <- purrr::imap_dfr(cox_zph_objects, function(zph_object, model_name) {
  zph_table <- as.data.frame(zph_object$table) %>%
    rownames_to_column("term") %>%
    as_tibble()
  
  names(zph_table) <- c("term", "chisq", "df", "p_value")
  
  zph_table %>%
    mutate(
      model_spec = model_name,
      model_label = unname(cox_model_labels[[model_name]]),
      .before = 1
    )
})

cox_zph_summary <- cox_zph_qc %>%
  mutate(
    is_global_test = term == "GLOBAL",
    proportional_hazards_flag = p_value < 0.05
  ) %>%
  group_by(model_spec, model_label) %>%
  summarise(
    n_term_tests = sum(!is_global_test),
    n_term_flags = sum(proportional_hazards_flag & !is_global_test),
    global_p_value = p_value[is_global_test][1],
    global_test_flag = global_p_value < 0.05,
    .groups = "drop"
  )

cox_zph_checks <- list(
  cox_zph_all_models_evaluated = length(cox_zph_objects) == length(cox_models),
  cox_zph_all_models_have_global_test = cox_zph_qc %>%
    group_by(model_spec) %>%
    summarise(has_global = any(term == "GLOBAL"), .groups = "drop") %>%
    pull(has_global) %>%
    all(),
  cox_zph_statistics_finite = all(is.finite(cox_zph_qc$chisq)),
  cox_zph_degrees_of_freedom_valid = all(is.finite(cox_zph_qc$df) & cox_zph_qc$df > 0),
  cox_zph_p_values_valid = all(cox_zph_qc$p_value >= 0 & cox_zph_qc$p_value <= 1)
)


stopifnot(isTRUE(all(unlist(cox_zph_checks))))

model_qc$survival$proportional_hazards <- list(
  tests = cox_zph_qc,
  summary = cox_zph_summary,
  interpretation = paste0(
    "Schoenfeld-residual tests evaluate whether each adjusted hazard ratio is constant over time. ",
    "P values below 0.05 are assumption flags rather than computational failures."
  )
)

model_qc$checks <- c(model_qc$checks, cox_zph_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Adjusted Cox proportional-hazards diagnostics added to: ", model_qc_path)

# =============================================================================
# Sensitivity analysis: Cox models stratified by race
# =============================================================================

# Race-stratified sensitivity models allow separate baseline hazards by race while estimating
# common coefficients for the focal exposures and remaining covariates. Their interpretation
# uses the proportional-hazards diagnostics generated from the current analysis sample.

cox_race_stratum_qc <- cox_analysis_data %>%
  group_by(race_category) %>%
  summarise(
    n = n(), n_events = sum(event == 1L), n_censored = sum(event == 0L),
    event_percent = 100 * mean(event == 1L), .groups = "drop"
  )

cox_race_stratum_checks <- list(
  all_race_strata_nonempty = all(cox_race_stratum_qc$n > 0L),
  race_stratum_counts_match_cox_n = sum(cox_race_stratum_qc$n) == nrow(cox_analysis_data),
  race_stratum_events_match_cox_events = sum(cox_race_stratum_qc$n_events) == cox_sample_qc$n_events,
  race_stratum_censoring_matches_cox_censoring = sum(cox_race_stratum_qc$n_censored) == cox_sample_qc$n_censored,
  at_least_two_race_strata_present = nrow(cox_race_stratum_qc) >= 2L
)


stopifnot(isTRUE(all(unlist(cox_race_stratum_checks))))

cox_model_race_stratified <- coxph(
  Surv(surv_time, event) ~ pbwpfvc + vtpbw + age_at_admission +
    sex_category + sf10 + sofa_total + strata(race_category),
  data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
)

cox_model_pfvc_race_stratified <- coxph(
  Surv(surv_time, event) ~ pfvc + vtpbw + age_at_admission +
    sex_category + sf10 + sofa_total + strata(race_category),
  data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
)

cox_model_abwpfvc_race_stratified <- coxph(
  Surv(surv_time, event) ~ abwpfvc + vtabw + age_at_admission +
    sex_category + sf10 + sofa_total + strata(race_category),
  data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
)

cox_model_nchspfvc_race_stratified <- coxph(
  Surv(surv_time, event) ~ nchspfvc + vtnchs + age_at_admission +
    sex_category + sf10 + sofa_total + strata(race_category),
  data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
)

cox_race_stratified_models <- list(
  vtpbw_pbwpfvc = cox_model_race_stratified,
  vtnchs_nchspfvc = cox_model_nchspfvc_race_stratified,
  vtabw_abwpfvc = cox_model_abwpfvc_race_stratified,
  vtpbw_pfvc = cox_model_pfvc_race_stratified
)

cox_race_stratified_qc <- purrr::imap_dfr(cox_race_stratified_models, function(model, model_name) {
  coefficient_table <- summary(model)$coefficients
  model_matrix <- model$x
  n_patients_model <- model$n
  n_events_model <- model$nevent
  
  tibble(
    model_spec = model_name, label = unname(cox_model_labels[[model_name]]),
    formula = paste(deparse(formula(model)), collapse = " "),
    n_obs = n_patients_model, n_events = n_events_model, n_censored = n_patients_model - n_events_model,
    n_coefficients = length(coef(model)), events_per_coefficient = n_events_model / length(coef(model)),
    iterations = model$iter, model_matrix_columns = ncol(model_matrix),
    model_matrix_rank = qr(model_matrix)$rank, full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    coefficients_finite = all(is.finite(coef(model))),
    standard_errors_finite = all(is.finite(coefficient_table[, "se(coef)"])),
    maximum_absolute_coefficient = max(abs(coef(model))),
    maximum_standard_error = max(coefficient_table[, "se(coef)"]),
    log_likelihood = as.numeric(logLik(model)), AIC = AIC(model)
  )
})

cox_race_stratified_checks <- list(
  all_expected_race_stratified_models_fitted = length(cox_race_stratified_models) == 4L,
  race_stratified_models_use_common_n = all(cox_race_stratified_qc$n_obs == nrow(cox_analysis_data)),
  race_stratified_models_have_expected_events = all(cox_race_stratified_qc$n_events == cox_sample_qc$n_events),
  race_stratified_models_have_expected_censoring = all(cox_race_stratified_qc$n_censored == cox_sample_qc$n_censored),
  race_stratified_model_matrices_full_rank = all(cox_race_stratified_qc$full_rank),
  race_stratified_coefficients_finite = all(cox_race_stratified_qc$coefficients_finite),
  race_stratified_standard_errors_finite = all(cox_race_stratified_qc$standard_errors_finite),
  race_stratified_log_likelihoods_finite = all(is.finite(cox_race_stratified_qc$log_likelihood)),
  race_stratified_aic_values_finite = all(is.finite(cox_race_stratified_qc$AIC))
)


stopifnot(isTRUE(all(unlist(cox_race_stratified_checks))))

# Compare focal exposure estimates between the primary and race-stratified models.
cox_focal_terms <- list(
  vtpbw_pbwpfvc = c("vtpbw", "pbwpfvc"),
  vtnchs_nchspfvc = c("vtnchs", "nchspfvc"),
  vtabw_abwpfvc = c("vtabw", "abwpfvc"),
  vtpbw_pfvc = c("vtpbw", "pfvc")
)

cox_primary_focal_results <- purrr::imap_dfr(cox_models, function(model, model_name) {
  broom::tidy(model, conf.int = TRUE) %>%
    filter(term %in% cox_focal_terms[[model_name]]) %>%
    transmute(
      model_spec = model_name, term,
      primary_log_hr = estimate, primary_standard_error = std.error,
      primary_conf_low_log = conf.low, primary_conf_high_log = conf.high,
      primary_p_value = p.value
    )
})

cox_stratified_focal_results <- purrr::imap_dfr(cox_race_stratified_models, function(model, model_name) {
  broom::tidy(model, conf.int = TRUE) %>%
    filter(term %in% cox_focal_terms[[model_name]]) %>%
    transmute(
      model_spec = model_name, term,
      stratified_log_hr = estimate, stratified_standard_error = std.error,
      stratified_conf_low_log = conf.low, stratified_conf_high_log = conf.high,
      stratified_p_value = p.value
    )
})

cox_race_stratified_comparison <- cox_primary_focal_results %>%
  inner_join(cox_stratified_focal_results, by = c("model_spec", "term")) %>%
  mutate(
    model_label = unname(cox_model_labels[model_spec]),
    primary_hr = exp(primary_log_hr),
    primary_conf_low = exp(primary_conf_low_log),
    primary_conf_high = exp(primary_conf_high_log),
    stratified_hr = exp(stratified_log_hr),
    stratified_conf_low = exp(stratified_conf_low_log),
    stratified_conf_high = exp(stratified_conf_high_log),
    change_in_log_hr = stratified_log_hr - primary_log_hr,
    ratio_of_hazard_ratios = exp(change_in_log_hr)
  ) %>%
  select(
    model_spec, model_label, term,
    primary_hr, primary_conf_low, primary_conf_high, primary_p_value,
    stratified_hr, stratified_conf_low, stratified_conf_high, stratified_p_value,
    change_in_log_hr, ratio_of_hazard_ratios
  )

cox_race_stratified_comparison_checks <- list(
  all_focal_terms_compared = nrow(cox_race_stratified_comparison) == 8L,
  primary_focal_estimates_finite = all(is.finite(c(
    cox_race_stratified_comparison$primary_hr,
    cox_race_stratified_comparison$primary_conf_low,
    cox_race_stratified_comparison$primary_conf_high
  ))),
  stratified_focal_estimates_finite = all(is.finite(c(
    cox_race_stratified_comparison$stratified_hr,
    cox_race_stratified_comparison$stratified_conf_low,
    cox_race_stratified_comparison$stratified_conf_high
  ))),
  hazard_ratio_comparisons_finite = all(is.finite(cox_race_stratified_comparison$ratio_of_hazard_ratios))
)


stopifnot(isTRUE(all(unlist(cox_race_stratified_comparison_checks))))

# Reassess proportional hazards after allowing the baseline hazard to differ by race.
cox_race_stratified_zph_objects <- purrr::map(cox_race_stratified_models, function(model) {
  survival::cox.zph(model, transform = "km", terms = TRUE, singledf = FALSE)
})

cox_race_stratified_zph_qc <- purrr::imap_dfr(
  cox_race_stratified_zph_objects,
  function(zph_object, model_name) {
    zph_table <- as.data.frame(zph_object$table) %>%
      rownames_to_column("term") %>%
      as_tibble()
    
    names(zph_table) <- c("term", "chisq", "df", "p_value")
    
    zph_table %>%
      mutate(
        model_spec = model_name,
        model_label = unname(cox_model_labels[[model_name]]),
        .before = 1
      )
  }
)

cox_race_stratified_zph_summary <- cox_race_stratified_zph_qc %>%
  mutate(
    is_global_test = term == "GLOBAL",
    proportional_hazards_flag = p_value < 0.05
  ) %>%
  group_by(model_spec, model_label) %>%
  summarise(
    n_term_tests = sum(!is_global_test),
    n_term_flags = sum(proportional_hazards_flag & !is_global_test),
    global_p_value = p_value[is_global_test][1],
    global_test_flag = global_p_value < 0.05,
    .groups = "drop"
  )

cox_race_stratified_zph_checks <- list(
  all_race_stratified_models_evaluated = length(cox_race_stratified_zph_objects) == 4L,
  all_stratified_models_have_global_test = cox_race_stratified_zph_qc %>%
    group_by(model_spec) %>%
    summarise(has_global = any(term == "GLOBAL"), .groups = "drop") %>%
    pull(has_global) %>%
    all(),
  stratified_zph_statistics_finite = all(is.finite(cox_race_stratified_zph_qc$chisq)),
  stratified_zph_degrees_of_freedom_valid = all(
    is.finite(cox_race_stratified_zph_qc$df) & cox_race_stratified_zph_qc$df > 0
  ),
  stratified_zph_p_values_valid = all(
    cox_race_stratified_zph_qc$p_value >= 0 & cox_race_stratified_zph_qc$p_value <= 1
  )
)


stopifnot(isTRUE(all(unlist(cox_race_stratified_zph_checks))))

# Render a sensitivity table. Race has no coefficient because it defines the baseline-hazard strata.
cox_sensitivity_labels <- c(
  vtpbw = "VT/McLaren PBW (mL/kg)",
  vtnchs = "VT/NCHS weight (mL/kg)",
  vtabw = "VT/ABW (mL/kg)",
  pfvc = "Predicted FVC (L)",
  pbwpfvc = "McLaren PBW/PFVC (kg/L)",
  nchspfvc = "NCHS weight/PFVC (kg/L)",
  abwpfvc = "ABW/PFVC (kg/L)",
  age_at_admission = "Age at admission (per year)",
  sex_category = "Sex",
  sf10 = "SF ratio (per 10)",
  sofa_total = "pSOFA total"
)

cox_race_stratified_tables <- purrr::map(cox_race_stratified_models, function(model) {
  tbl_regression(model, exponentiate = TRUE) %>%
    modify_table_body(function(table_body) {
      table_body %>%
        mutate(
          replacement_label = unname(cox_sensitivity_labels[variable]),
          label = if_else(row_type == "label" & !is.na(replacement_label), replacement_label, label)
        ) %>%
        select(-replacement_label)
    }) %>%
    bold_p()
})

cox_race_stratified_table <- tbl_merge(
  cox_race_stratified_tables, tab_spanner = unname(cox_model_labels), quiet = TRUE
) %>%
  modify_caption(paste0("**Race-stratified Cox sensitivity models for 60-day mortality — ", site_name, "**"))

cox_race_stratified_table_gt <- cox_race_stratified_table %>%
  as_gt() %>%
  tab_source_note(
    source_note = paste0(
      "Models are stratified by race category, allowing a separate baseline hazard for each race group. ",
      "Race therefore has no estimated hazard ratio. All models use the same complete-case sample: N = ",
      nrow(cox_analysis_data), ", including ", cox_sample_qc$n_events, " deaths."
    )
  )

cox_race_stratified_output_path <- file.path(
  table_dir, paste0("cox_race_stratified_sensitivity_", site_name, ".html")
)

gt::gtsave(cox_race_stratified_table_gt, cox_race_stratified_output_path)

cox_race_stratified_output_exists <- file.exists(cox_race_stratified_output_path)
cox_race_stratified_output_size <- if (cox_race_stratified_output_exists) {
  file.info(cox_race_stratified_output_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

cox_race_stratified_output_checks <- list(
  race_stratified_output_exists = cox_race_stratified_output_exists,
  race_stratified_output_nonempty = cox_race_stratified_output_exists &&
    is.finite(cox_race_stratified_output_size) &&
    cox_race_stratified_output_size > 0
)

model_qc$survival$race_stratified_sensitivity <- list(
  race_stratum_summary = cox_race_stratum_qc,
  model_diagnostics = cox_race_stratified_qc,
  focal_estimate_comparison = cox_race_stratified_comparison,
  proportional_hazards = list(
    tests = cox_race_stratified_zph_qc,
    summary = cox_race_stratified_zph_summary
  ),
  output = list(
    path = cox_race_stratified_output_path,
    size_bytes = cox_race_stratified_output_size
  )
)

model_qc$checks <- c(
  model_qc$checks,
  cox_race_stratum_checks,
  cox_race_stratified_checks,
  cox_race_stratified_comparison_checks,
  cox_race_stratified_zph_checks,
  cox_race_stratified_output_checks
)

model_qc$metadata$updated_at <- Sys.time()

duplicate_model_qc_check_names <- unique(names(model_qc$checks)[duplicated(names(model_qc$checks))])

if (length(duplicate_model_qc_check_names) > 0L) {
  stop(
    "model_qc contains duplicate check name(s): ",
    paste(duplicate_model_qc_check_names, collapse = ", ")
  )
}

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Race-stratified Cox sensitivity analysis saved to: ", cox_race_stratified_output_path)

message(
  "Adjusted Cox models validated: ", nrow(cox_analysis_data),
  " patients and ", cox_sample_qc$n_events, " deaths."
)

# =============================================================================
# Kaplan-Meier analyses by weight/PFVC ratio tercile
# =============================================================================

km_ratio_specs <- list(
  mclaren = list(
    ratio = "pbwpfvc",
    label = "McLaren PBW/PFVC",
    output_stub = "mclaren_pbwpfvc",
    palette = c("#CC5555", "#266CAE", "#D47E0E")
  ),
  nchs = list(
    ratio = "nchspfvc",
    label = "NCHS weight/PFVC",
    output_stub = "nchs_weight_pfvc",
    palette = c("#009E73", "#0072B2", "#D55E00")
  ),
  actual = list(
    ratio = "abwpfvc",
    label = "ABW/PFVC",
    output_stub = "abw_pfvc",
    palette = c("#CC79A7", "#0072B2", "#E69F00")
  )
)

km_results <- purrr::imap(km_ratio_specs, function(spec, method_name) {
  ratio_name <- spec$ratio
  tercile_name <- paste0(ratio_name, "_tercile")
  
  km_data <- surv_data %>%
    filter(is.finite(.data[[ratio_name]])) %>%
    mutate(
      !!tercile_name := factor(
        ntile(.data[[ratio_name]], 3L),
        levels = 1:3,
        labels = c("T1 (Low)", "T2 (Mid)", "T3 (High)")
      )
    )
  
  km_n <- nrow(km_data)
  km_n_deaths <- sum(km_data$event == 1L)
  km_n_censored <- sum(km_data$event == 0L)
  
  km_summary <- km_data %>%
    group_by(.data[[tercile_name]]) %>%
    summarise(
      n = n(),
      n_deaths = sum(event == 1L),
      n_censored = sum(event == 0L),
      mortality_percent = 100 * mean(event == 1L),
      minimum_ratio = min(.data[[ratio_name]]),
      maximum_ratio = max(.data[[ratio_name]]),
      .groups = "drop"
    ) %>%
    rename(tercile = all_of(tercile_name))
  
  km_formula <- as.formula(
    paste0("survival::Surv(surv_time, event) ~ ", tercile_name)
  )
  
  fit <- do.call(
    survival::survfit,
    list(
      formula = km_formula,
      data = km_data
    )
  )
  
  fit$call$formula <- km_formula
  
  logrank <- do.call(
    survival::survdiff,
    list(
      formula = km_formula,
      data = km_data
    )
  )
  
  logrank_df <- length(logrank$n) - 1L
  logrank_p_value <- pchisq(logrank$chisq, df = logrank_df, lower.tail = FALSE)
  
  checks <- list(
    ratio_values_complete_and_finite =
      all(is.finite(km_data[[ratio_name]])),
    analysis_sample_nonempty =
      km_n > 0L,
    analysis_sample_has_events =
      km_n_deaths > 0L,
    analysis_sample_has_censored_observations =
      km_n_censored > 0L,
    tercile_assignment_complete =
      !any(is.na(km_data[[tercile_name]])),
    tercile_counts_match_analysis_n =
      sum(km_summary$n) == km_n,
    tercile_deaths_match_analysis_deaths =
      sum(km_summary$n_deaths) == km_n_deaths,
    tercile_censoring_matches_analysis_censoring =
      sum(km_summary$n_censored) == km_n_censored,
    all_terciles_nonempty =
      all(km_summary$n > 0L),
    fit_has_three_strata =
      length(fit$strata) == 3L,
    logrank_statistic_finite =
      is.finite(logrank$chisq),
    logrank_df_valid =
      logrank_df == 2L,
    logrank_p_value_valid =
      is.finite(logrank_p_value) &&
      logrank_p_value >= 0 &&
      logrank_p_value <= 1
  )
  
  stopifnot(isTRUE(all(unlist(checks))))
  
  plot <- ggsurvplot(
    fit,
    data = km_data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    palette = spec$palette,
    xlab = "Days from admission",
    ylab = "Survival probability",
    title = paste0("Kaplan-Meier survival by ", spec$label, " tercile"),
    legend.title = paste0(spec$label, " tercile"),
    legend.labs = c("Lowest tercile (T1)", "Middle tercile (T2)", "Highest tercile (T3)"),
    ggtheme = theme_minimal(),
    tables.theme = theme_cleantable()
  )
  
  output_path <- file.path(
    figure_dir,
    paste0("km_curves_", spec$output_stub, "_", site_name, ".pdf")
  )
  
  km_arranged_plot <- survminer::arrange_ggsurvplots(
    list(plot),
    print = FALSE,
    ncol = 1, nrow = 1,
    surv.plot.height = 0.75,
    risk.table.height = 0.25
  )
  
  if (length(km_arranged_plot) != 1L) {
    stop("Expected exactly one arranged Kaplan-Meier page.")
  }
  
  if (!grid::is.grob(km_arranged_plot[[1L]])) {
    stop("The arranged Kaplan-Meier page is not a drawable grid object.")
  }
  
  save_script04_plot(
    filename = output_path, plot = km_arranged_plot[[1L]],
    device = "pdf", width = 10, height = 8, units = "in", bg = "white"
  )
  
  output_exists <- file.exists(output_path)
  output_size <- if (output_exists) {
    file.info(output_path, extra_cols = FALSE)$size[[1]]
  } else {
    NA_real_
  }
  
  stopifnot(
    output_exists,
    is.finite(output_size),
    output_size > 0
  )
  
  list(
    method = method_name,
    ratio = ratio_name,
    label = spec$label,
    sample = tibble(
      n_survival_endpoint = nrow(surv_data),
      n_analyzed = km_n,
      n_excluded_missing_or_nonfinite_ratio = nrow(surv_data) - km_n,
      n_deaths = km_n_deaths,
      n_censored = km_n_censored
    ),
    summary = km_summary,
    fit = fit,
    logrank = tibble(
      chisq = unname(logrank$chisq),
      df = logrank_df,
      p_value = logrank_p_value
    ),
    checks = tibble(
      check = names(checks),
      passed = as.logical(unlist(checks, use.names = FALSE))
    ),
    plot = plot,
    output = tibble(
      path = output_path,
      exists = output_exists,
      size_bytes = output_size
    )
  )
})

km_sample_qc <- purrr::imap_dfr(km_results, function(result, method_name) {
  result$sample %>%
    mutate(
      method = method_name,
      ratio = result$ratio,
      label = result$label,
      .before = 1
    )
})

km_method_summary <- purrr::imap_dfr(km_results, function(result, method_name) {
  result$summary %>%
    mutate(
      method = method_name,
      ratio = result$ratio,
      label = result$label,
      .before = 1
    )
})

km_logrank_qc <- purrr::imap_dfr(km_results, function(result, method_name) {
  result$logrank %>%
    mutate(
      method = method_name,
      ratio = result$ratio,
      label = result$label,
      .before = 1
    )
})

km_output_qc <- purrr::imap_dfr(km_results, function(result, method_name) {
  result$output %>%
    mutate(
      method = method_name,
      ratio = result$ratio,
      label = result$label,
      .before = 1
    )
})

km_checks_table <- purrr::imap_dfr(km_results, function(result, method_name) {
  result$checks %>%
    mutate(
      method = method_name,
      ratio = result$ratio,
      label = result$label,
      .before = 1
    )
})


stopifnot(
  nrow(km_sample_qc) == 3L,
  nrow(km_output_qc) == 3L,
  all(km_sample_qc$n_analyzed > 0L),
  all(km_sample_qc$n_deaths > 0L),
  all(km_sample_qc$n_censored > 0L),
  isTRUE(all(km_checks_table$passed)),
  all(km_output_qc$exists),
  all(is.finite(km_output_qc$size_bytes) & km_output_qc$size_bytes > 0)
)

model_qc$survival$kaplan_meier <- list(
  sample_summary = km_sample_qc,
  method_summary = km_method_summary,
  logrank_tests = km_logrank_qc,
  checks = km_checks_table,
  outputs = km_output_qc
)

# Preserve these objects as the McLaren-primary aliases for downstream compatibility.
km_fit <- km_results$mclaren$fit
km_plot <- km_results$mclaren$plot
km_output_path <- km_results$mclaren$output$path

model_qc$checks <- c(
  model_qc$checks,
  list(
    km_all_method_checks_passed = isTRUE(all(km_checks_table$passed)),
    km_all_outputs_exist = all(km_output_qc$exists),
    km_all_outputs_nonempty = all(
      is.finite(km_output_qc$size_bytes) & km_output_qc$size_bytes > 0
    )
  )
)

model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Independent Kaplan-Meier analyses saved for McLaren PBW, NCHS weight, and ABW.")

# =============================================================================
# Cox-model text output
# =============================================================================

cox_summary_output_path <- file.path(figure_dir, paste0("cox_model_summary_", site_name, ".txt"))

cox_summary_lines <- c(
  "=== Cox model: VT/McLaren PBW + McLaren PBW/PFVC ===",
  capture.output(print(summary(cox_model))),
  "",
  "=== Cox model: VT/NCHS weight + NCHS weight/PFVC ===",
  capture.output(print(summary(cox_model_nchspfvc))),
  "",
  "=== Cox model: VT/ABW + ABW/PFVC ===",
  capture.output(print(summary(cox_model_abwpfvc))),
  "",
  "=== Cox model: VT/McLaren PBW + PFVC ===",
  capture.output(print(summary(cox_model_pfvc))),
  "",
  "=== Proportional-hazards diagnostics ===",
  capture.output(print(cox_zph_qc, n = Inf, width = Inf)),
  "",
  "=== Kaplan-Meier log-rank tests by weight method ===",
  capture.output(print(km_logrank_qc, n = Inf, width = Inf))
)

writeLines(cox_summary_lines, con = cox_summary_output_path)

cox_summary_output_exists <- file.exists(cox_summary_output_path)

cox_summary_output_size <- if (cox_summary_output_exists) {
  file.info(cox_summary_output_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

survival_output_checks <- list(
  km_outputs_exist = all(km_output_qc$exists),
  km_outputs_nonempty = all(
    is.finite(km_output_qc$size_bytes) & km_output_qc$size_bytes > 0
  ),
  cox_summary_output_exists = cox_summary_output_exists,
  cox_summary_output_nonempty = cox_summary_output_exists &&
    is.finite(cox_summary_output_size) &&
    cox_summary_output_size > 0
)

model_qc$survival$output <- list(
  km_outputs = km_output_qc,
  cox_summary_path = cox_summary_output_path,
  cox_summary_size_bytes = cox_summary_output_size
)

model_qc$checks <- c(model_qc$checks, survival_output_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message(
  "Kaplan-Meier curves saved to: ",
  paste(km_output_qc$path, collapse = "; ")
)
message("Cox-model summary saved to: ", cox_summary_output_path)

# =============================================================================
# 4f2. Demographic-bias models (outcome ~ demographics only)
# =============================================================================

# How each dosing, driving-pressure, and elastance metric varies by demographics (the algorithmic-bias view). Age is
# modeled per year; height and SF ratio are scaled per 10 units. Dosing and elastance-normalized outcomes are z-scored
# within site, so coefficients are expressed in outcome SD units. Driving pressure remains on the raw cmH2O scale and
# additionally adjusts for BMI; 60-day mortality is logistic (OR). Reference categories are male and white.

demo_required_vars <- c(
  "age_at_admission", "sex_category", "race_category", "height_cm", "sf_ratio", "sofa_total", "bmi",
  "pbw", "pbw_nchs", "abw", "pfvc", "vtpbw", "vtnchs", "vtabw", "vtpfvc",
  "pbwpfvc", "nchspfvc", "abwpfvc", "ers_pbw", "ers_pfvc", "mortality_event_60", "dp"
)

missing_demo_vars <- setdiff(demo_required_vars, names(cross_sectional))

stopifnot(length(missing_demo_vars) == 0L)

demo_data <- cross_sectional %>%
  mutate(
    height10 = height_cm / 10,
    sf10 = sf_ratio / 10,
    vtpbw_z = as.numeric(scale(vtpbw)),
    vtabw_z = as.numeric(scale(vtabw)),
    vtnchs_z = as.numeric(scale(vtnchs)),
    vtpfvc_z = as.numeric(scale(vtpfvc)),
    abwpfvc_z = as.numeric(scale(abwpfvc)),
    pbwpfvc_z = as.numeric(scale(pbwpfvc)),
    nchspfvc_z = as.numeric(scale(nchspfvc)),
    ers_pbw_z = as.numeric(scale(ers_pbw)),
    ers_pfvc_z = as.numeric(scale(ers_pfvc))
  )

demo_covars <- "age_at_admission + sex_category + race_category + height10 + sf10 + sofa_total"
demo_covars_bmi <- paste(demo_covars, "+ bmi")

demo_numeric_vars <- c(
  "age_at_admission", "height_cm", "height10", "sf_ratio", "sf10", "sofa_total", "bmi",
  "pbw", "pbw_nchs", "abw", "pfvc", "vtpbw", "vtnchs", "vtabw", "vtpfvc",
  "pbwpfvc", "nchspfvc", "abwpfvc", "ers_pbw", "ers_pfvc", "mortality_event_60", "dp"
)

demo_input_qc <- purrr::map_dfr(demo_numeric_vars, function(v) {
  x <- demo_data[[v]]
  finite_values <- x[is.finite(x)]
  
  tibble(
    variable = v, class = paste(class(x), collapse = "/"), n_total = length(x),
    n_nonmissing = sum(!is.na(x)), n_missing = sum(is.na(x)), pct_missing = 100 * mean(is.na(x)),
    n_unique_nonmissing = n_distinct(x, na.rm = TRUE), n_infinite = sum(is.infinite(x), na.rm = TRUE),
    minimum = if (length(finite_values) > 0L) min(finite_values) else NA_real_,
    maximum = if (length(finite_values) > 0L) max(finite_values) else NA_real_
  )
})

demo_standardization_qc <- tibble(
  original_variable = c(
    "vtpbw", "vtnchs", "vtabw", "vtpfvc", "pbwpfvc", "nchspfvc", "abwpfvc", "ers_pbw", "ers_pfvc"
  ),
  standardized_variable = c(
    "vtpbw_z", "vtnchs_z", "vtabw_z", "vtpfvc_z", "pbwpfvc_z", "nchspfvc_z", "abwpfvc_z",
    "ers_pbw_z", "ers_pfvc_z"
  )
) %>%
  mutate(
    n_nonmissing = purrr::map_int(standardized_variable, ~ sum(!is.na(demo_data[[.x]]))),
    mean_standardized = purrr::map_dbl(standardized_variable, ~ mean(demo_data[[.x]], na.rm = TRUE)),
    sd_standardized = purrr::map_dbl(standardized_variable, ~ sd(demo_data[[.x]], na.rm = TRUE))
  )

height10_compared <- is.finite(demo_data$height10) & is.finite(demo_data$height_cm)
sf10_demo_compared <- is.finite(demo_data$sf10) & is.finite(demo_data$sf_ratio)

height10_max_absolute_difference <- if (any(height10_compared)) {
  max(abs(demo_data$height10[height10_compared] - demo_data$height_cm[height10_compared] / 10))
} else {
  NA_real_
}

sf10_demo_max_absolute_difference <- if (any(sf10_demo_compared)) {
  max(abs(demo_data$sf10[sf10_demo_compared] - demo_data$sf_ratio[sf10_demo_compared] / 10))
} else {
  NA_real_
}

demo_input_checks <- list(
  demographic_variables_present = length(missing_demo_vars) == 0L,
  demographic_numeric_variables_numeric = all(purrr::map_lgl(demo_data[demo_numeric_vars], is.numeric)),
  demographic_numeric_variables_have_no_infinite_values = !any(demo_input_qc$n_infinite > 0L),
  demographic_numeric_variables_have_variation = all(demo_input_qc$n_unique_nonmissing >= 2L),
  demographic_sex_levels_valid = identical(levels(demo_data$sex_category), c("Male", "Female")),
  demographic_race_levels_valid = identical(levels(demo_data$race_category), c("WHITE", "BLACK", "OTHER")),
  height10_derivation_valid = is.finite(height10_max_absolute_difference) &&
    height10_max_absolute_difference <= 1e-12,
  demographic_sf10_derivation_valid = is.finite(sf10_demo_max_absolute_difference) &&
    sf10_demo_max_absolute_difference <= 1e-12,
  standardized_means_are_zero = all(abs(demo_standardization_qc$mean_standardized) <= 1e-12),
  standardized_sds_are_one = all(abs(demo_standardization_qc$sd_standardized - 1) <= 1e-12),
  standardized_nonmissing_counts_preserved = all(
    demo_standardization_qc$n_nonmissing ==
      purrr::map_int(demo_standardization_qc$original_variable, ~ sum(!is.na(demo_data[[.x]])))
  )
)


stopifnot(isTRUE(all(unlist(demo_input_checks))))

model_qc$demographic_bias <- list(
  required_variables = demo_required_vars,
  missing_variables = missing_demo_vars,
  input_summary = demo_input_qc,
  standardization = demo_standardization_qc,
  scale_derivations = tibble(
    derived_variable = c("height10", "sf10"),
    max_absolute_difference = c(height10_max_absolute_difference, sf10_demo_max_absolute_difference)
  )
)

model_qc$checks <- c(model_qc$checks, demo_input_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Demographic-bias input QC added to: ", model_qc_path)

# Each entry contains a fitted model and metadata for the unified long-format table.
demo_models <- list()

# Z-scored linear outcomes. The elastance-normalized outcomes are driving-pressure derivatives, so they are adjusted for
# BMI to match the original paper; VT/PBW and VT/PFVC are not.
demo_z_outcomes <- list(
  vtpbw = list(z_var = "vtpbw_z", label = "VT/McLaren PBW", bmi_adjusted = FALSE),
  vtnchs = list(z_var = "vtnchs_z", label = "VT/NCHS weight", bmi_adjusted = FALSE),
  vtabw = list(z_var = "vtabw_z", label = "VT/ABW", bmi_adjusted = FALSE),
  vtpfvc = list(z_var = "vtpfvc_z", label = "VT/PFVC (%)", bmi_adjusted = FALSE),
  pbwpfvc = list(z_var = "pbwpfvc_z", label = "McLaren PBW/PFVC", bmi_adjusted = FALSE),
  nchspfvc = list(z_var = "nchspfvc_z", label = "NCHS weight/PFVC", bmi_adjusted = FALSE),
  abwpfvc = list(z_var = "abwpfvc_z", label = "ABW/PFVC", bmi_adjusted = FALSE),
  ers_pbw = list(z_var = "ers_pbw_z", label = "Ers x PBW", bmi_adjusted = TRUE),
  ers_pfvc = list(z_var = "ers_pfvc_z", label = "Ers x PFVC", bmi_adjusted = TRUE)
)

for (v in names(demo_z_outcomes)) {
  outcome_spec <- demo_z_outcomes[[v]]
  cov_v <- if (outcome_spec$bmi_adjusted) demo_covars_bmi else demo_covars
  fstr <- paste(outcome_spec$z_var, "~", cov_v)
  
  demo_models[[outcome_spec$label]] <- list(
    model = lm(as.formula(fstr), data = demo_data),
    type = "Beta", family = "linear", formula = fstr,
    original_outcome = v, standardized_outcome = outcome_spec$z_var
  )
}

# Sixty-day mortality logistic model. BMI is not included in this demographic model.
if (has_mortality_variation) {
  fstr <- paste("mortality_event_60 ~", demo_covars)
  
  demo_models[["60-day Mortality"]] <- list(
    model = glm(as.formula(fstr), data = demo_data, family = binomial),
    type = "OR", family = "logistic", formula = fstr
  )
}

# Driving-pressure estimate on its raw scale, with BMI adjustment and a minimum sample-size guard.
# Plateau-derived static estimates and PIP-derived fallback estimates retain their upstream provenance.
if (sum(!is.na(demo_data$dp)) >= 30L) {
  fstr <- paste("dp ~", demo_covars_bmi)
  
  demo_models[["Driving pressure"]] <- list(
    model = lm(as.formula(fstr), data = demo_data),
    type = "Beta", family = "linear", formula = fstr
  )
} else {
  message("Skipping driving-pressure demographic-bias model: < 30 observations")
}

# =============================================================================
# Model QC: demographic-bias model samples and fits
# =============================================================================

demo_model_qc <- purrr::imap_dfr(demo_models, function(model_info, outcome_name) {
  model <- model_info$model
  model_frame <- model.frame(model)
  model_matrix <- model.matrix(model)
  coefficient_table <- summary(model)$coefficients
  outcome_values <- model.response(model_frame)
  is_logistic <- identical(model_info$family, "logistic")
  aliased_coefficients <- is.na(coef(model))
  n_parameters_excluding_intercept <- model$rank - 1L
  
  tibble(
    outcome = outcome_name, model_family = model_info$family, estimate_type = model_info$type,
    formula = model_info$formula, n_source = nrow(demo_data), n_obs = nrow(model_frame),
    n_excluded = n_source - n_obs, pct_analyzed = 100 * n_obs / n_source,
    n_events = if (is_logistic) sum(outcome_values == 1L) else NA_integer_,
    n_nonevents = if (is_logistic) sum(outcome_values == 0L) else NA_integer_,
    n_parameters_excluding_intercept = n_parameters_excluding_intercept,
    events_per_parameter = if (is_logistic) {
      sum(outcome_values == 1L) / n_parameters_excluding_intercept
    } else {
      NA_real_
    },
    residual_df = df.residual(model), model_matrix_columns = ncol(model_matrix),
    model_matrix_rank = qr(model_matrix)$rank, full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    n_aliased_coefficients = sum(aliased_coefficients),
    coefficients_finite = all(is.finite(coef(model)[!aliased_coefficients])),
    standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
    fitted_values_finite = all(is.finite(fitted(model))),
    residuals_finite = all(is.finite(residuals(model))),
    converged = if (is_logistic) isTRUE(model$converged) else TRUE,
    boundary = if (is_logistic) isTRUE(model$boundary) else FALSE,
    minimum_fitted_value = min(fitted(model)), maximum_fitted_value = max(fitted(model)),
    n_extreme_fitted_probabilities = if (is_logistic) {
      sum(fitted(model) < 1e-6 | fitted(model) > 1 - 1e-6)
    } else {
      NA_integer_
    }
  )
})

demo_model_sample_qc <- demo_model_qc %>%
  select(
    outcome, model_family, estimate_type, n_source, n_obs, n_excluded, pct_analyzed,
    n_events, n_nonevents, n_parameters_excluding_intercept, events_per_parameter
  )

demo_logistic_qc <- demo_model_qc %>%
  filter(model_family == "logistic")

demo_model_checks <- list(
  demographic_models_fitted = length(demo_models) > 0L,
  demographic_model_matrices_full_rank = all(demo_model_qc$full_rank),
  demographic_models_have_no_aliased_coefficients = all(demo_model_qc$n_aliased_coefficients == 0L),
  demographic_model_coefficients_finite = all(demo_model_qc$coefficients_finite),
  demographic_model_standard_errors_finite = all(demo_model_qc$standard_errors_finite),
  demographic_model_fitted_values_finite = all(demo_model_qc$fitted_values_finite),
  demographic_model_residuals_finite = all(demo_model_qc$residuals_finite),
  demographic_models_have_positive_residual_df = all(demo_model_qc$residual_df > 0L),
  demographic_logistic_models_converged = nrow(demo_logistic_qc) == 0L || all(demo_logistic_qc$converged),
  no_demographic_logistic_model_at_boundary = nrow(demo_logistic_qc) == 0L || !any(demo_logistic_qc$boundary),
  no_extreme_demographic_logistic_probabilities = nrow(demo_logistic_qc) == 0L ||
    all(demo_logistic_qc$n_extreme_fitted_probabilities == 0L),
  demographic_models_have_nonempty_samples = all(demo_model_qc$n_obs > 0L)
)


stopifnot(isTRUE(all(unlist(demo_model_checks))))

model_qc$demographic_bias$model_diagnostics <- demo_model_qc
model_qc$demographic_bias$sample_summary <- demo_model_sample_qc
model_qc$demographic_bias$sparse_event_warning <- demo_logistic_qc %>%
  select(outcome, n_obs, n_events, n_parameters_excluding_intercept, events_per_parameter) %>%
  mutate(events_per_parameter_below_10 = events_per_parameter < 10)

model_qc$checks <- c(model_qc$checks, demo_model_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Demographic-bias model diagnostics added to: ", model_qc_path)

# Combined table: rows are covariates and columns are outcomes.
demo_star <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

demo_term_labels <- c(
  age_at_admission = "Age at Admission (per year)",
  sex_categoryFemale = "Female vs male",
  race_categoryOTHER = "Other vs white",
  race_categoryBLACK = "Black vs white",
  height10 = "Height (per 10 cm)",
  sf10 = "SF ratio (per 10)",
  sofa_total = "pSOFA",
  bmi = "BMI"
)

demo_outcome_order <- c(
  "VT/McLaren PBW", "VT/NCHS weight", "VT/ABW", "VT/PFVC (%)",
  "McLaren PBW/PFVC", "NCHS weight/PFVC", "ABW/PFVC",
  "60-day Mortality", "Driving pressure", "Ers x PBW", "Ers x PFVC"
)

demo_model_results <- purrr::imap_dfr(demo_models, function(model_info, outcome) {
  tidy_results <- if (identical(model_info$family, "linear")) {
    tidy_lm_hc3(model_info$model)
  } else {
    broom::tidy(model_info$model, conf.int = TRUE, exponentiate = TRUE)
  }
  
  tidy_results %>%
    filter(term != "(Intercept)") %>%
    mutate(
      outcome = outcome, estimate_type = model_info$type,
      inference = if_else(model_info$family == "linear", "HC3 robust", "Model-based"),
      .before = 1
    )
})

demo_cells <- demo_model_results %>%
  transmute(
    outcome,
    term,
    cell = sprintf("%.2f [%.2f, %.2f]%s", estimate, conf.low, conf.high, demo_star(p.value))
  )

demo_n_row <- tibble(
  term_label = "N",
  outcome = names(demo_models),
  cell = map_chr(demo_models, ~ format(stats::nobs(.x$model), big.mark = ","))
)

demo_table <- demo_cells %>%
  mutate(term_label = recode(term, !!!demo_term_labels)) %>%
  select(term_label, outcome, cell) %>%
  bind_rows(demo_n_row) %>%
  mutate(
    term_label = factor(term_label, levels = c(unname(demo_term_labels), "N")),
    outcome = factor(outcome, levels = demo_outcome_order)
  ) %>%
  arrange(term_label, outcome) %>%
  pivot_wider(names_from = outcome, values_from = cell) %>%
  arrange(term_label)

demo_gt <- demo_table %>%
  mutate(term_label = as.character(term_label)) %>%
  gt::gt(rowname_col = "term_label") %>%
  gt::tab_header(
    title = "Demographic variation in dosing, driving pressure, and elastance",
    subtitle = paste0(
      site_name,
      " — z-scored outcomes (SD units) except 60-day mortality (OR) and driving pressure (cmH2O); ",
      "linear-model intervals use HC3 robust standard errors; SF ratio per 10 units; reference = male, white"
    )
  ) %>%
  gt::sub_missing(missing_text = "")

demo_output_path <- file.path(table_dir, paste0("table_demographic_bias_", site_name, ".html"))

gt::gtsave(demo_gt, demo_output_path)

demo_output_exists <- file.exists(demo_output_path)
demo_output_size <- if (demo_output_exists) {
  file.info(demo_output_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

demo_output_checks <- list(
  demographic_bias_output_exists = demo_output_exists,
  demographic_bias_output_nonempty = demo_output_exists &&
    is.finite(demo_output_size) &&
    demo_output_size > 0,
  demographic_bias_results_nonempty = nrow(demo_model_results) > 0L,
  demographic_bias_estimates_finite = all(is.finite(demo_model_results$estimate)),
  demographic_bias_confidence_intervals_finite = all(
    is.finite(demo_model_results$conf.low) & is.finite(demo_model_results$conf.high)
  ),
  demographic_bias_p_values_valid = all(
    demo_model_results$p.value >= 0 & demo_model_results$p.value <= 1
  )
)

model_qc$demographic_bias$coefficient_results <- demo_model_results
model_qc$demographic_bias$rendered_table <- demo_table
model_qc$demographic_bias$output <- list(
  path = demo_output_path,
  exists = demo_output_exists,
  size_bytes = demo_output_size
)

model_qc$checks <- c(model_qc$checks, demo_output_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Demographic-bias table written (", length(demo_models), " outcome models): ", demo_output_path)

# =============================================================================
# 4f3. Predicted FVC vs predicted body weight
# =============================================================================

# PFVC is regressed on PBW and demographic variables in the broad eligible cohort. This evaluates whether demographic
# information explains variation in predicted FVC beyond that represented by PBW.

broad_pfvc_path <- file.path(output_dir, "analysis_broad_pfvc.parquet")

stopifnot(file.exists(broad_pfvc_path))

broad_pfvc <- read_parquet(broad_pfvc_path)

broad_pfvc_required_vars <- c(
  "pfvc", "pbw", "pbw_nchs", "abw", "age_at_admission", "sex_category", "race_category"
)
missing_broad_pfvc_vars <- setdiff(broad_pfvc_required_vars, names(broad_pfvc))

stopifnot(length(missing_broad_pfvc_vars) == 0L)

broad_pfvc_sex_values <- sort(unique(na.omit(as.character(broad_pfvc$sex_category))))
broad_pfvc_race_values <- sort(unique(na.omit(as.character(broad_pfvc$race_category))))

broad_pfvc <- broad_pfvc %>%
  mutate(
    sex_category = factor(sex_category, levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER"))
  )

broad_pfvc_input_qc <- purrr::map_dfr(broad_pfvc_required_vars, function(v) {
  x <- broad_pfvc[[v]]
  finite_values <- if (is.numeric(x)) x[is.finite(x)] else numeric(0)
  
  tibble(
    variable = v, class = paste(class(x), collapse = "/"), n_total = length(x),
    n_nonmissing = sum(!is.na(x)), n_missing = sum(is.na(x)), pct_missing = 100 * mean(is.na(x)),
    n_unique_nonmissing = n_distinct(x, na.rm = TRUE),
    n_infinite = if (is.numeric(x)) sum(is.infinite(x), na.rm = TRUE) else NA_integer_,
    minimum = if (length(finite_values) > 0L) min(finite_values) else NA_real_,
    maximum = if (length(finite_values) > 0L) max(finite_values) else NA_real_
  )
})

broad_pfvc_analysis_data <- broad_pfvc %>%
  mutate(.source_row = row_number()) %>%
  select(.source_row, all_of(broad_pfvc_required_vars)) %>%
  drop_na(all_of(broad_pfvc_required_vars))

broad_pfvc_sample_qc <- tibble(
  n_source = nrow(broad_pfvc), n_analyzed = nrow(broad_pfvc_analysis_data),
  n_excluded = n_source - n_analyzed, pct_analyzed = 100 * n_analyzed / n_source
)

broad_pfvc_input_checks <- list(
  broad_pfvc_file_exists = file.exists(broad_pfvc_path),
  broad_pfvc_required_variables_present = length(missing_broad_pfvc_vars) == 0L,
  broad_pfvc_numeric_variables_numeric = all(
    purrr::map_lgl(broad_pfvc[c("pfvc", "pbw", "pbw_nchs", "abw", "age_at_admission")], is.numeric)
  ),
  broad_pfvc_numeric_variables_have_no_infinite_values = !any(
    broad_pfvc_input_qc$n_infinite > 0L, na.rm = TRUE
  ),
  broad_pfvc_variables_have_variation = all(broad_pfvc_input_qc$n_unique_nonmissing >= 2L),
  broad_pfvc_sex_values_valid = all(broad_pfvc_sex_values %in% c("Female", "Male")),
  broad_pfvc_race_values_valid = all(broad_pfvc_race_values %in% c("BLACK", "OTHER", "WHITE")),
  broad_pfvc_analysis_sample_nonempty = nrow(broad_pfvc_analysis_data) > 0L,
  broad_pfvc_analysis_sample_complete = all(complete.cases(broad_pfvc_analysis_data[broad_pfvc_required_vars])),
  broad_pfvc_positive = all(broad_pfvc_analysis_data$pfvc > 0),
  broad_pbw_positive = all(broad_pfvc_analysis_data$pbw > 0),
  broad_nchs_weight_positive = all(broad_pfvc_analysis_data$pbw_nchs > 0),
  broad_actual_weight_positive = all(broad_pfvc_analysis_data$abw > 0),
  broad_pfvc_source_rows_unique = !anyDuplicated(broad_pfvc_analysis_data$.source_row)
)


stopifnot(isTRUE(all(unlist(broad_pfvc_input_checks))))

model_qc$pfvc_vs_pbw <- list(
  input_file = broad_pfvc_path,
  required_variables = broad_pfvc_required_vars,
  missing_variables = missing_broad_pfvc_vars,
  input_summary = broad_pfvc_input_qc,
  sample_summary = broad_pfvc_sample_qc
)

model_qc$checks <- c(model_qc$checks, broad_pfvc_input_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("PFVC-versus-PBW input QC added to: ", model_qc_path)

pfvc_vs_pbw_model <- lm(
  pfvc ~ pbw + age_at_admission + sex_category + race_category,
  data = broad_pfvc_analysis_data, x = TRUE, y = TRUE, model = TRUE
)

pfvc_vs_nchs_model <- lm(
  pfvc ~ pbw_nchs + age_at_admission + sex_category + race_category,
  data = broad_pfvc_analysis_data, x = TRUE, y = TRUE, model = TRUE
)

pfvc_vs_abw_model <- lm(
  pfvc ~ abw + age_at_admission + sex_category + race_category,
  data = broad_pfvc_analysis_data, x = TRUE, y = TRUE, model = TRUE
)

pfvc_vs_weight_models <- list(
  mclaren = pfvc_vs_pbw_model,
  nchs = pfvc_vs_nchs_model,
  actual = pfvc_vs_abw_model
)

pfvc_vs_weight_labels <- c(
  mclaren = "McLaren PBW",
  nchs = "NCHS weight",
  actual = "Actual body weight"
)

pfvc_vs_weight_terms <- c(
  mclaren = "pbw",
  nchs = "pbw_nchs",
  actual = "abw"
)

pfvc_vs_weight_model_qc <- purrr::imap_dfr(pfvc_vs_weight_models, function(model, method_name) {
  model_matrix <- model.matrix(model)
  coefficient_table <- summary(model)$coefficients
  aliased_coefficients <- is.na(coef(model))
  
  tibble(
    method = method_name,
    method_label = unname(pfvc_vs_weight_labels[[method_name]]),
    weight_term = unname(pfvc_vs_weight_terms[[method_name]]),
    formula = paste(deparse(formula(model)), collapse = " "),
    n_obs = nobs(model),
    residual_df = df.residual(model),
    model_matrix_columns = ncol(model_matrix),
    model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    n_aliased_coefficients = sum(aliased_coefficients),
    coefficients_finite = all(is.finite(coef(model)[!aliased_coefficients])),
    standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
    fitted_values_finite = all(is.finite(fitted(model))),
    residuals_finite = all(is.finite(residuals(model))),
    r_squared = summary(model)$r.squared,
    adjusted_r_squared = summary(model)$adj.r.squared,
    AIC = AIC(model)
  )
})

pfvc_vs_weight_hc3_results <- purrr::imap_dfr(pfvc_vs_weight_models, function(model, method_name) {
  tidy_lm_hc3(model) %>%
    mutate(
      method = method_name,
      method_label = unname(pfvc_vs_weight_labels[[method_name]]),
      weight_term = unname(pfvc_vs_weight_terms[[method_name]]),
      n_obs = nobs(model),
      .before = 1
    )
})

pfvc_vs_weight_checks <- list(
  all_weight_models_fitted = length(pfvc_vs_weight_models) == 3L,
  all_weight_models_use_common_n = all(
    pfvc_vs_weight_model_qc$n_obs == nrow(broad_pfvc_analysis_data)
  ),
  all_weight_model_matrices_full_rank = all(pfvc_vs_weight_model_qc$full_rank),
  no_weight_models_have_aliased_coefficients = all(
    pfvc_vs_weight_model_qc$n_aliased_coefficients == 0L
  ),
  all_weight_model_coefficients_finite = all(pfvc_vs_weight_model_qc$coefficients_finite),
  all_weight_model_standard_errors_finite = all(pfvc_vs_weight_model_qc$standard_errors_finite),
  all_weight_model_fitted_values_finite = all(pfvc_vs_weight_model_qc$fitted_values_finite),
  all_weight_model_residuals_finite = all(pfvc_vs_weight_model_qc$residuals_finite),
  all_weight_models_have_positive_residual_df = all(pfvc_vs_weight_model_qc$residual_df > 0L),
  all_weight_model_aic_values_finite = all(is.finite(pfvc_vs_weight_model_qc$AIC)),
  all_weight_model_hc3_results_finite = all(is.finite(c(
    pfvc_vs_weight_hc3_results$estimate,
    pfvc_vs_weight_hc3_results$std.error,
    pfvc_vs_weight_hc3_results$conf.low,
    pfvc_vs_weight_hc3_results$conf.high,
    pfvc_vs_weight_hc3_results$p.value
  )))
)


stopifnot(isTRUE(all(unlist(pfvc_vs_weight_checks))))

# =============================================================================
# Model QC: predicted FVC versus predicted body weight
# =============================================================================

pfvc_vs_pbw_model_matrix <- model.matrix(pfvc_vs_pbw_model)
pfvc_vs_pbw_coefficient_table <- summary(pfvc_vs_pbw_model)$coefficients
pfvc_vs_pbw_aliased <- is.na(coef(pfvc_vs_pbw_model))
pfvc_vs_pbw_residuals <- residuals(pfvc_vs_pbw_model)
pfvc_vs_pbw_studentized_residuals <- rstudent(pfvc_vs_pbw_model)
pfvc_vs_pbw_cooks_distance <- cooks.distance(pfvc_vs_pbw_model)
pfvc_vs_pbw_leverage <- hatvalues(pfvc_vs_pbw_model)
pfvc_vs_pbw_n <- nobs(pfvc_vs_pbw_model)
pfvc_vs_pbw_p <- pfvc_vs_pbw_model$rank

pfvc_vs_pbw_residual_z <- as.numeric(scale(pfvc_vs_pbw_residuals))
pfvc_vs_pbw_shapiro <- safe_shapiro_test(pfvc_vs_pbw_residuals)
pfvc_vs_pbw_bp <- lmtest::bptest(pfvc_vs_pbw_model)
pfvc_vs_pbw_reset <- lmtest::resettest(pfvc_vs_pbw_model, power = 2:3, type = "fitted")

pfvc_vs_pbw_model_qc <- tibble(
  formula = paste(deparse(formula(pfvc_vs_pbw_model)), collapse = " "),
  n_source = nrow(broad_pfvc), n_obs = pfvc_vs_pbw_n,
  n_excluded = n_source - n_obs, pct_analyzed = 100 * n_obs / n_source,
  n_parameters = pfvc_vs_pbw_model$rank, residual_df = df.residual(pfvc_vs_pbw_model),
  model_matrix_columns = ncol(pfvc_vs_pbw_model_matrix),
  model_matrix_rank = qr(pfvc_vs_pbw_model_matrix)$rank,
  full_rank = qr(pfvc_vs_pbw_model_matrix)$rank == ncol(pfvc_vs_pbw_model_matrix),
  n_aliased_coefficients = sum(pfvc_vs_pbw_aliased),
  coefficients_finite = all(is.finite(coef(pfvc_vs_pbw_model)[!pfvc_vs_pbw_aliased])),
  standard_errors_finite = all(is.finite(pfvc_vs_pbw_coefficient_table[, "Std. Error"])),
  fitted_values_finite = all(is.finite(fitted(pfvc_vs_pbw_model))),
  residuals_finite = all(is.finite(pfvc_vs_pbw_residuals)),
  residual_standard_error = sigma(pfvc_vs_pbw_model),
  r_squared = summary(pfvc_vs_pbw_model)$r.squared,
  adjusted_r_squared = summary(pfvc_vs_pbw_model)$adj.r.squared,
  residual_skewness = mean(pfvc_vs_pbw_residual_z^3),
  residual_excess_kurtosis = mean(pfvc_vs_pbw_residual_z^4) - 3,
  shapiro_p_value = pfvc_vs_pbw_shapiro$p.value,
  breusch_pagan_p_value = pfvc_vs_pbw_bp$p.value,
  reset_p_value = pfvc_vs_pbw_reset$p.value,
  maximum_absolute_studentized_residual = max(abs(pfvc_vs_pbw_studentized_residuals)),
  n_absolute_studentized_residual_gt_3 = sum(abs(pfvc_vs_pbw_studentized_residuals) > 3),
  cooks_threshold = 4 / pfvc_vs_pbw_n,
  maximum_cooks_distance = max(pfvc_vs_pbw_cooks_distance),
  n_cooks_distance_above_threshold = sum(pfvc_vs_pbw_cooks_distance > 4 / pfvc_vs_pbw_n),
  leverage_threshold = 2 * pfvc_vs_pbw_p / pfvc_vs_pbw_n,
  maximum_leverage = max(pfvc_vs_pbw_leverage),
  n_leverage_above_threshold = sum(pfvc_vs_pbw_leverage > 2 * pfvc_vs_pbw_p / pfvc_vs_pbw_n)
)

pfvc_vs_pbw_model_checks <- list(
  pfvc_vs_pbw_model_uses_expected_n = pfvc_vs_pbw_n == nrow(broad_pfvc_analysis_data),
  pfvc_vs_pbw_model_matrix_full_rank = pfvc_vs_pbw_model_qc$full_rank,
  pfvc_vs_pbw_model_has_no_aliased_coefficients = pfvc_vs_pbw_model_qc$n_aliased_coefficients == 0L,
  pfvc_vs_pbw_coefficients_finite = pfvc_vs_pbw_model_qc$coefficients_finite,
  pfvc_vs_pbw_standard_errors_finite = pfvc_vs_pbw_model_qc$standard_errors_finite,
  pfvc_vs_pbw_fitted_values_finite = pfvc_vs_pbw_model_qc$fitted_values_finite,
  pfvc_vs_pbw_residuals_finite = pfvc_vs_pbw_model_qc$residuals_finite,
  pfvc_vs_pbw_positive_residual_df = pfvc_vs_pbw_model_qc$residual_df > 0L,
  pfvc_vs_pbw_diagnostics_finite = all(is.finite(c(
    pfvc_vs_pbw_model_qc$residual_skewness,
    pfvc_vs_pbw_model_qc$residual_excess_kurtosis,
    pfvc_vs_pbw_model_qc$shapiro_p_value,
    pfvc_vs_pbw_model_qc$breusch_pagan_p_value,
    pfvc_vs_pbw_model_qc$reset_p_value,
    pfvc_vs_pbw_model_qc$maximum_absolute_studentized_residual,
    pfvc_vs_pbw_model_qc$maximum_cooks_distance,
    pfvc_vs_pbw_model_qc$maximum_leverage
  )))
)

pfvc_vs_pbw_diagnostic_flags <- tibble(
  nonnormal_residual_test = pfvc_vs_pbw_model_qc$shapiro_p_value < 0.05,
  heteroskedasticity_test = pfvc_vs_pbw_model_qc$breusch_pagan_p_value < 0.05,
  specification_test = pfvc_vs_pbw_model_qc$reset_p_value < 0.05,
  large_studentized_residuals_present = pfvc_vs_pbw_model_qc$n_absolute_studentized_residual_gt_3 > 0L,
  influential_observations_present = pfvc_vs_pbw_model_qc$n_cooks_distance_above_threshold > 0L,
  high_leverage_observations_present = pfvc_vs_pbw_model_qc$n_leverage_above_threshold > 0L
)


stopifnot(isTRUE(all(unlist(pfvc_vs_pbw_model_checks))))

model_qc$pfvc_vs_pbw$model_diagnostics <- pfvc_vs_pbw_model_qc
model_qc$pfvc_vs_pbw$diagnostic_flags <- pfvc_vs_pbw_diagnostic_flags
model_qc$pfvc_vs_pbw$interpretation <- paste0(
  "PFVC and PBW are predicted quantities derived partly from demographic variables. This model characterizes ",
  "relationships among prediction formulas and should not be interpreted as a model of measured lung size."
)

model_qc$checks <- c(model_qc$checks, pfvc_vs_pbw_model_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("PFVC-versus-PBW model diagnostics added to: ", model_qc_path)

pfvc_vs_pbw_formula <- "pfvc ~ pbw + age_at_admission + sex_category + race_category"

pfvc_vs_pbw_labels <- c(
  pbw = "Predicted body weight (kg)",
  age_at_admission = "Age at admission (per year)",
  sex_category = "Sex",
  race_category = "Race"
)

pfvc_vs_pbw_results_hc3 <- tidy_lm_hc3(pfvc_vs_pbw_model)

pfvc_vs_pbw_table <- tbl_regression(pfvc_vs_pbw_model, tidy_fun = tidy_lm_hc3) %>%
  modify_table_body(function(table_body) {
    table_body %>%
      mutate(
        replacement_label = unname(pfvc_vs_pbw_labels[variable]),
        label = if_else(row_type == "label" & !is.na(replacement_label), replacement_label, label)
      ) %>%
      select(-replacement_label)
  }) %>%
  bold_p() %>%
  modify_caption(paste0(
    "**Predicted FVC versus predicted body weight in the broad cohort — ",
    site_name, "**"
  ))

pfvc_vs_pbw_table_gt <- pfvc_vs_pbw_table %>%
  as_gt() %>%
  tab_source_note(
    source_note = paste0(
      "N = ", nobs(pfvc_vs_pbw_model),
      ". Coefficients are unstandardized linear-regression estimates. Confidence intervals and P values use HC3 ",
      "heteroskedasticity-robust standard errors. PFVC and PBW are predicted quantities derived partly from demographic ",
      "variables; estimates characterize relationships among prediction formulas rather than measured lung function."
    )
  )

pfvc_vs_pbw_output_path <- file.path(table_dir, paste0("table_pfvc_vs_pbw_", site_name, ".html"))

gt::gtsave(pfvc_vs_pbw_table_gt, pfvc_vs_pbw_output_path)

pfvc_vs_pbw_output_exists <- file.exists(pfvc_vs_pbw_output_path)
pfvc_vs_pbw_output_size <- if (pfvc_vs_pbw_output_exists) {
  file.info(pfvc_vs_pbw_output_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

pfvc_vs_pbw_output_checks <- list(
  pfvc_vs_pbw_output_exists = pfvc_vs_pbw_output_exists,
  pfvc_vs_pbw_output_nonempty = pfvc_vs_pbw_output_exists &&
    is.finite(pfvc_vs_pbw_output_size) &&
    pfvc_vs_pbw_output_size > 0,
  pfvc_vs_pbw_hc3_estimates_finite = all(is.finite(pfvc_vs_pbw_results_hc3$estimate)),
  pfvc_vs_pbw_hc3_standard_errors_finite = all(is.finite(pfvc_vs_pbw_results_hc3$std.error)),
  pfvc_vs_pbw_hc3_confidence_intervals_finite = all(
    is.finite(pfvc_vs_pbw_results_hc3$conf.low) & is.finite(pfvc_vs_pbw_results_hc3$conf.high)
  ),
  pfvc_vs_pbw_hc3_p_values_valid = all(
    pfvc_vs_pbw_results_hc3$p.value >= 0 & pfvc_vs_pbw_results_hc3$p.value <= 1
  )
)

model_qc$pfvc_vs_pbw$HC3_coefficient_results <- pfvc_vs_pbw_results_hc3
model_qc$pfvc_vs_pbw$output <- list(
  path = pfvc_vs_pbw_output_path,
  size_bytes = pfvc_vs_pbw_output_size
)

model_qc$checks <- c(model_qc$checks, pfvc_vs_pbw_output_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message(
  "PFVC-versus-PBW table written using HC3 inference (N = ",
  nobs(pfvc_vs_pbw_model), "): ", pfvc_vs_pbw_output_path
)

pfvc_vs_weight_regression_labels <- c(
  pbw = "McLaren predicted body weight (kg)",
  pbw_nchs = "NCHS Standards weight (kg)",
  abw = "Actual body weight (kg)",
  age_at_admission = "Age at admission (per year)",
  sex_category = "Sex",
  race_category = "Race"
)

pfvc_vs_weight_tables <- purrr::map(pfvc_vs_weight_models, function(model) {
  tbl_regression(model, tidy_fun = tidy_lm_hc3) %>%
    modify_table_body(function(table_body) {
      table_body %>%
        mutate(
          replacement_label = unname(pfvc_vs_weight_regression_labels[variable]),
          label = if_else(row_type == "label" & !is.na(replacement_label), replacement_label, label)
        ) %>%
        select(-replacement_label)
    }) %>%
    bold_p()
})

pfvc_vs_weight_table <- tbl_merge(
  pfvc_vs_weight_tables,
  tab_spanner = unname(pfvc_vs_weight_labels),
  quiet = TRUE
) %>%
  modify_caption(paste0(
    "**Predicted FVC versus weight-method estimates in the broad cohort — ",
    site_name, "**"
  ))

pfvc_vs_weight_table_gt <- pfvc_vs_weight_table %>%
  as_gt() %>%
  tab_source_note(
    source_note = paste0(
      "All three models use the same complete-case sample: N = ",
      nrow(broad_pfvc_analysis_data),
      ". Coefficients are unstandardized linear-regression estimates. ",
      "Confidence intervals and P values use HC3 heteroskedasticity-robust standard errors."
    )
  )

pfvc_vs_weight_output_path <- file.path(
  table_dir,
  paste0("table_pfvc_vs_weight_methods_", site_name, ".html")
)

gt::gtsave(pfvc_vs_weight_table_gt, pfvc_vs_weight_output_path)

pfvc_vs_weight_output_exists <- file.exists(pfvc_vs_weight_output_path)
pfvc_vs_weight_output_size <- if (pfvc_vs_weight_output_exists) {
  file.info(pfvc_vs_weight_output_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

pfvc_vs_weight_output_checks <- list(
  pfvc_vs_weight_output_exists = pfvc_vs_weight_output_exists,
  pfvc_vs_weight_output_nonempty = pfvc_vs_weight_output_exists &&
    is.finite(pfvc_vs_weight_output_size) &&
    pfvc_vs_weight_output_size > 0
)

model_qc$pfvc_vs_pbw$weight_method_comparison <- list(
  model_diagnostics = pfvc_vs_weight_model_qc,
  HC3_coefficient_results = pfvc_vs_weight_hc3_results,
  output = list(
    path = pfvc_vs_weight_output_path,
    size_bytes = pfvc_vs_weight_output_size
  )
)

model_qc$checks <- c(
  model_qc$checks,
  pfvc_vs_weight_checks,
  pfvc_vs_weight_output_checks
)

model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

# =============================================================================
# Sensitivity analysis: nonlinear PFVC-versus-PBW functional form
# =============================================================================
# Natural splines allow PFVC to vary nonlinearly with PBW and age while retaining sex and race
# as categorical covariates. This sensitivity analysis evaluates functional form; diagnostic
# results are generated from the current analysis sample.

pfvc_vs_pbw_spline_model <- lm(
  pfvc ~ ns(pbw, df = 4) + ns(age_at_admission, df = 4) + sex_category + race_category,
  data = broad_pfvc_analysis_data, x = TRUE, y = TRUE, model = TRUE
)

pfvc_vs_pbw_spline_matrix <- model.matrix(pfvc_vs_pbw_spline_model)
pfvc_vs_pbw_spline_coefficients <- summary(pfvc_vs_pbw_spline_model)$coefficients
pfvc_vs_pbw_spline_aliased <- is.na(coef(pfvc_vs_pbw_spline_model))
pfvc_vs_pbw_spline_residuals <- residuals(pfvc_vs_pbw_spline_model)
pfvc_vs_pbw_spline_studentized <- rstudent(pfvc_vs_pbw_spline_model)
pfvc_vs_pbw_spline_cooks <- cooks.distance(pfvc_vs_pbw_spline_model)
pfvc_vs_pbw_spline_leverage <- hatvalues(pfvc_vs_pbw_spline_model)
pfvc_vs_pbw_spline_residual_z <- as.numeric(scale(pfvc_vs_pbw_spline_residuals))

pfvc_vs_pbw_spline_shapiro <- safe_shapiro_test(pfvc_vs_pbw_spline_residuals)
pfvc_vs_pbw_spline_bp <- lmtest::bptest(pfvc_vs_pbw_spline_model)
pfvc_vs_pbw_spline_reset <- lmtest::resettest(pfvc_vs_pbw_spline_model, power = 2:3, type = "fitted")

pfvc_vs_pbw_spline_qc <- tibble(
  formula = paste(deparse(formula(pfvc_vs_pbw_spline_model)), collapse = " "),
  n_obs = nobs(pfvc_vs_pbw_spline_model), n_parameters = pfvc_vs_pbw_spline_model$rank,
  residual_df = df.residual(pfvc_vs_pbw_spline_model),
  model_matrix_columns = ncol(pfvc_vs_pbw_spline_matrix),
  model_matrix_rank = qr(pfvc_vs_pbw_spline_matrix)$rank,
  full_rank = qr(pfvc_vs_pbw_spline_matrix)$rank == ncol(pfvc_vs_pbw_spline_matrix),
  n_aliased_coefficients = sum(pfvc_vs_pbw_spline_aliased),
  coefficients_finite = all(is.finite(coef(pfvc_vs_pbw_spline_model)[!pfvc_vs_pbw_spline_aliased])),
  standard_errors_finite = all(is.finite(pfvc_vs_pbw_spline_coefficients[, "Std. Error"])),
  residual_standard_error = sigma(pfvc_vs_pbw_spline_model),
  r_squared = summary(pfvc_vs_pbw_spline_model)$r.squared,
  adjusted_r_squared = summary(pfvc_vs_pbw_spline_model)$adj.r.squared,
  AIC = AIC(pfvc_vs_pbw_spline_model),
  residual_skewness = mean(pfvc_vs_pbw_spline_residual_z^3),
  residual_excess_kurtosis = mean(pfvc_vs_pbw_spline_residual_z^4) - 3,
  shapiro_p_value = pfvc_vs_pbw_spline_shapiro$p.value,
  breusch_pagan_p_value = pfvc_vs_pbw_spline_bp$p.value,
  reset_p_value = pfvc_vs_pbw_spline_reset$p.value,
  maximum_absolute_studentized_residual = max(abs(pfvc_vs_pbw_spline_studentized)),
  n_absolute_studentized_residual_gt_3 = sum(abs(pfvc_vs_pbw_spline_studentized) > 3),
  maximum_cooks_distance = max(pfvc_vs_pbw_spline_cooks),
  n_cooks_distance_above_threshold = sum(
    pfvc_vs_pbw_spline_cooks > 4 / nobs(pfvc_vs_pbw_spline_model)
  ),
  maximum_leverage = max(pfvc_vs_pbw_spline_leverage),
  n_leverage_above_threshold = sum(
    pfvc_vs_pbw_spline_leverage >
      2 * pfvc_vs_pbw_spline_model$rank / nobs(pfvc_vs_pbw_spline_model)
  )
)

pfvc_vs_pbw_model_comparison <- bind_rows(
  pfvc_vs_pbw_model_qc %>%
    transmute(
      model = "Linear PBW and age", n_obs, n_parameters, residual_df,
      residual_standard_error, r_squared, adjusted_r_squared, AIC = AIC(pfvc_vs_pbw_model),
      residual_skewness, residual_excess_kurtosis, shapiro_p_value,
      breusch_pagan_p_value, reset_p_value,
      maximum_absolute_studentized_residual, maximum_cooks_distance, maximum_leverage
    ),
  pfvc_vs_pbw_spline_qc %>%
    transmute(
      model = "Natural splines for PBW and age", n_obs, n_parameters, residual_df,
      residual_standard_error, r_squared, adjusted_r_squared, AIC,
      residual_skewness, residual_excess_kurtosis, shapiro_p_value,
      breusch_pagan_p_value, reset_p_value,
      maximum_absolute_studentized_residual, maximum_cooks_distance, maximum_leverage
    )
) %>%
  mutate(delta_AIC_vs_linear = AIC - AIC[model == "Linear PBW and age"])

pfvc_vs_pbw_nested_comparison <- anova(pfvc_vs_pbw_model, pfvc_vs_pbw_spline_model) %>%
  as.data.frame() %>%
  rownames_to_column("model_number") %>%
  as_tibble()

pfvc_vs_pbw_spline_checks <- list(
  spline_model_uses_primary_sample = nobs(pfvc_vs_pbw_spline_model) == nobs(pfvc_vs_pbw_model),
  spline_model_matrix_full_rank = pfvc_vs_pbw_spline_qc$full_rank,
  spline_model_has_no_aliased_coefficients = pfvc_vs_pbw_spline_qc$n_aliased_coefficients == 0L,
  spline_coefficients_finite = pfvc_vs_pbw_spline_qc$coefficients_finite,
  spline_standard_errors_finite = pfvc_vs_pbw_spline_qc$standard_errors_finite,
  spline_residual_df_positive = pfvc_vs_pbw_spline_qc$residual_df > 0L,
  spline_aic_finite = is.finite(pfvc_vs_pbw_spline_qc$AIC),
  spline_diagnostics_finite = all(is.finite(c(
    pfvc_vs_pbw_spline_qc$residual_skewness,
    pfvc_vs_pbw_spline_qc$residual_excess_kurtosis,
    pfvc_vs_pbw_spline_qc$shapiro_p_value,
    pfvc_vs_pbw_spline_qc$breusch_pagan_p_value,
    pfvc_vs_pbw_spline_qc$reset_p_value,
    pfvc_vs_pbw_spline_qc$maximum_absolute_studentized_residual,
    pfvc_vs_pbw_spline_qc$maximum_cooks_distance,
    pfvc_vs_pbw_spline_qc$maximum_leverage
  )))
)


stopifnot(isTRUE(all(unlist(pfvc_vs_pbw_spline_checks))))

model_qc$pfvc_vs_pbw$nonlinear_sensitivity <- list(
  spline_degrees_of_freedom = 4L,
  spline_model_diagnostics = pfvc_vs_pbw_spline_qc,
  linear_vs_spline_comparison = pfvc_vs_pbw_model_comparison,
  nested_model_comparison = pfvc_vs_pbw_nested_comparison,
  interpretation = paste0(
    "The sensitivity model uses four-degree-of-freedom natural splines for PBW and age. ",
    "Spline coefficients are basis-function parameters and should not be interpreted individually."
  )
)

model_qc$checks <- c(model_qc$checks, pfvc_vs_pbw_spline_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Nonlinear PFVC-versus-PBW sensitivity analysis added to: ", model_qc_path)

# =============================================================================
# Sensitivity analysis: sex-specific nonlinear PFVC-versus-PBW relationships
# =============================================================================
# PBW and PFVC are calculated using sex-specific prediction equations. This model therefore allows the nonlinear
# relationships of PFVC with PBW and age to differ by sex while retaining race as a categorical covariate.

pfvc_vs_pbw_interaction_model <- lm(
  pfvc ~ ns(pbw, df = 4) * sex_category +
    ns(age_at_admission, df = 4) * sex_category + race_category,
  data = broad_pfvc_analysis_data, x = TRUE, y = TRUE, model = TRUE
)

pfvc_vs_pbw_interaction_matrix <- model.matrix(pfvc_vs_pbw_interaction_model)
pfvc_vs_pbw_interaction_coefficients <- summary(pfvc_vs_pbw_interaction_model)$coefficients
pfvc_vs_pbw_interaction_aliased <- is.na(coef(pfvc_vs_pbw_interaction_model))
pfvc_vs_pbw_interaction_residuals <- residuals(pfvc_vs_pbw_interaction_model)
pfvc_vs_pbw_interaction_studentized <- rstudent(pfvc_vs_pbw_interaction_model)
pfvc_vs_pbw_interaction_cooks <- cooks.distance(pfvc_vs_pbw_interaction_model)
pfvc_vs_pbw_interaction_leverage <- hatvalues(pfvc_vs_pbw_interaction_model)
pfvc_vs_pbw_interaction_residual_z <- as.numeric(scale(pfvc_vs_pbw_interaction_residuals))

pfvc_vs_pbw_interaction_shapiro <- safe_shapiro_test(pfvc_vs_pbw_interaction_residuals)
pfvc_vs_pbw_interaction_bp <- lmtest::bptest(pfvc_vs_pbw_interaction_model)
pfvc_vs_pbw_interaction_reset <- lmtest::resettest(
  pfvc_vs_pbw_interaction_model, power = 2:3, type = "fitted"
)

pfvc_vs_pbw_interaction_n <- nobs(pfvc_vs_pbw_interaction_model)

pfvc_vs_pbw_interaction_qc <- tibble(
  formula = paste(deparse(formula(pfvc_vs_pbw_interaction_model)), collapse = " "),
  n_obs = pfvc_vs_pbw_interaction_n, n_parameters = pfvc_vs_pbw_interaction_model$rank,
  residual_df = df.residual(pfvc_vs_pbw_interaction_model),
  model_matrix_columns = ncol(pfvc_vs_pbw_interaction_matrix),
  model_matrix_rank = qr(pfvc_vs_pbw_interaction_matrix)$rank,
  full_rank = qr(pfvc_vs_pbw_interaction_matrix)$rank == ncol(pfvc_vs_pbw_interaction_matrix),
  n_aliased_coefficients = sum(pfvc_vs_pbw_interaction_aliased),
  coefficients_finite = all(
    is.finite(coef(pfvc_vs_pbw_interaction_model)[!pfvc_vs_pbw_interaction_aliased])
  ),
  standard_errors_finite = all(is.finite(pfvc_vs_pbw_interaction_coefficients[, "Std. Error"])),
  residual_standard_error = sigma(pfvc_vs_pbw_interaction_model),
  r_squared = summary(pfvc_vs_pbw_interaction_model)$r.squared,
  adjusted_r_squared = summary(pfvc_vs_pbw_interaction_model)$adj.r.squared,
  AIC = AIC(pfvc_vs_pbw_interaction_model),
  residual_skewness = mean(pfvc_vs_pbw_interaction_residual_z^3),
  residual_excess_kurtosis = mean(pfvc_vs_pbw_interaction_residual_z^4) - 3,
  shapiro_p_value = pfvc_vs_pbw_interaction_shapiro$p.value,
  breusch_pagan_p_value = pfvc_vs_pbw_interaction_bp$p.value,
  reset_p_value = pfvc_vs_pbw_interaction_reset$p.value,
  maximum_absolute_studentized_residual = max(abs(pfvc_vs_pbw_interaction_studentized)),
  n_absolute_studentized_residual_gt_3 = sum(abs(pfvc_vs_pbw_interaction_studentized) > 3),
  maximum_cooks_distance = max(pfvc_vs_pbw_interaction_cooks),
  n_cooks_distance_above_threshold = sum(
    pfvc_vs_pbw_interaction_cooks > 4 / pfvc_vs_pbw_interaction_n
  ),
  maximum_leverage = max(pfvc_vs_pbw_interaction_leverage),
  n_leverage_above_threshold = sum(
    pfvc_vs_pbw_interaction_leverage >
      2 * pfvc_vs_pbw_interaction_model$rank / pfvc_vs_pbw_interaction_n
  )
)

pfvc_vs_pbw_interaction_comparison <- bind_rows(
  pfvc_vs_pbw_model_qc %>%
    transmute(
      model = "Linear PBW and age", n_obs, n_parameters, residual_df,
      residual_standard_error, r_squared, adjusted_r_squared, AIC = AIC(pfvc_vs_pbw_model),
      residual_skewness, residual_excess_kurtosis, breusch_pagan_p_value, reset_p_value
    ),
  pfvc_vs_pbw_spline_qc %>%
    transmute(
      model = "Natural splines for PBW and age", n_obs, n_parameters, residual_df,
      residual_standard_error, r_squared, adjusted_r_squared, AIC,
      residual_skewness, residual_excess_kurtosis, breusch_pagan_p_value, reset_p_value
    ),
  pfvc_vs_pbw_interaction_qc %>%
    transmute(
      model = "Sex-specific natural splines for PBW and age", n_obs, n_parameters, residual_df,
      residual_standard_error, r_squared, adjusted_r_squared, AIC,
      residual_skewness, residual_excess_kurtosis, breusch_pagan_p_value, reset_p_value
    )
) %>%
  mutate(delta_AIC_vs_linear = AIC - AIC[model == "Linear PBW and age"])

pfvc_vs_pbw_interaction_checks <- list(
  interaction_model_uses_primary_sample = pfvc_vs_pbw_interaction_n == nobs(pfvc_vs_pbw_model),
  interaction_model_matrix_full_rank = pfvc_vs_pbw_interaction_qc$full_rank,
  interaction_model_has_no_aliased_coefficients = pfvc_vs_pbw_interaction_qc$n_aliased_coefficients == 0L,
  interaction_coefficients_finite = pfvc_vs_pbw_interaction_qc$coefficients_finite,
  interaction_standard_errors_finite = pfvc_vs_pbw_interaction_qc$standard_errors_finite,
  interaction_residual_df_positive = pfvc_vs_pbw_interaction_qc$residual_df > 0L,
  interaction_aic_finite = is.finite(pfvc_vs_pbw_interaction_qc$AIC),
  interaction_diagnostics_finite = all(is.finite(c(
    pfvc_vs_pbw_interaction_qc$residual_skewness,
    pfvc_vs_pbw_interaction_qc$residual_excess_kurtosis,
    pfvc_vs_pbw_interaction_qc$shapiro_p_value,
    pfvc_vs_pbw_interaction_qc$breusch_pagan_p_value,
    pfvc_vs_pbw_interaction_qc$reset_p_value,
    pfvc_vs_pbw_interaction_qc$maximum_absolute_studentized_residual,
    pfvc_vs_pbw_interaction_qc$maximum_cooks_distance,
    pfvc_vs_pbw_interaction_qc$maximum_leverage
  )))
)


stopifnot(isTRUE(all(unlist(pfvc_vs_pbw_interaction_checks))))

model_qc$pfvc_vs_pbw$sex_specific_nonlinear_sensitivity <- list(
  model_diagnostics = pfvc_vs_pbw_interaction_qc,
  model_comparison = pfvc_vs_pbw_interaction_comparison,
  interpretation = paste0(
    "The model allows the four-degree-of-freedom spline relationships of PFVC with PBW and age to differ by sex. ",
    "Spline and interaction coefficients are basis-function parameters and should not be interpreted individually."
  )
)

model_qc$checks <- c(model_qc$checks, pfvc_vs_pbw_interaction_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Sex-specific nonlinear PFVC-versus-PBW sensitivity analysis added to: ", model_qc_path)

# =============================================================================
# 4g. Unified long-format regression results table
# =============================================================================

# One row per model and term. Estimates are reported on their natural scale: OR for logistic regression, HR for Cox and
# Fine-Gray regression, and Beta for linear regression. The stable column schema supports downstream cross-site analyses.

# Convert already-tidied model results to the stable long-format schema. This allows each model family to use its
# appropriate inference: model-based logistic and Cox estimates, HC3 linear estimates, and robust Fine-Gray estimates.
assemble_model_results <- function(tidy_results, estimate_type, analysis, model_spec, model_family, formula_str, n_obs,
                                   adjustment = "adjusted") {
  required_columns <- c("term", "estimate", "conf.low", "conf.high", "std.error", "statistic", "p.value")
  missing_columns <- setdiff(required_columns, names(tidy_results))
  
  if (length(missing_columns) > 0L) {
    stop("Tidied model results are missing required column(s): ", paste(missing_columns, collapse = ", "))
  }
  
  tidy_results %>%
    transmute(
      term,
      estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      std_error = std.error,
      statistic,
      p_value = p.value,
      estimate_type = estimate_type,
      analysis = analysis,
      model_spec = model_spec,
      model_family = model_family,
      adjustment = adjustment,
      formula = formula_str,
      n_obs = as.integer(n_obs)
    )
}

model_formula_string <- function(model) {
  paste(deparse(formula(model)), collapse = " ")
}

results_long <- list()

# -----------------------------------------------------------------------------
# Primary adjusted 60-day mortality models
# -----------------------------------------------------------------------------

if (!is.null(mortality_models)) {
  results_long <- c(
    results_long,
    purrr::imap(mortality_models, function(model, model_name) {
      tidy_results <- broom::tidy(
        model,
        conf.int = TRUE,
        exponentiate = TRUE
      )
      
      assemble_model_results(
        tidy_results = tidy_results,
        estimate_type = "OR",
        analysis = "Mortality",
        model_spec = unname(exposure_labels[[model_name]]),
        model_family = "logistic",
        formula_str = model_formula_string(model),
        n_obs = nrow(model.frame(model)),
        adjustment = "adjusted"
      )
    })
  )
}

# -----------------------------------------------------------------------------
# Primary adjusted continuous-outcome models
# -----------------------------------------------------------------------------
# Use the HC3 results underlying the rendered continuous-outcome tables rather than model-based OLS standard errors.

for (outcome_name in names(continuous_outcomes)) {
  outcome_label <- continuous_outcomes[[outcome_name]]$label
  
  results_long <- c(
    results_long,
    purrr::imap(continuous_models[[outcome_name]], function(model, model_name) {
      tidy_results <- continuous_hc3_qc_combined %>%
        filter(outcome_name == !!outcome_name, model_spec == !!model_name) %>%
        select(term, estimate, conf.low, conf.high, std.error, statistic, p.value)
      
      assemble_model_results(
        tidy_results = tidy_results,
        estimate_type = "Beta",
        analysis = outcome_label,
        model_spec = unname(exposure_labels[[model_name]]),
        model_family = "linear",
        formula_str = model_formula_string(model),
        n_obs = nrow(model.frame(model)),
        adjustment = "adjusted"
      )
    })
  )
}

# -----------------------------------------------------------------------------
# Primary adjusted 28-day Fine-Gray models
# -----------------------------------------------------------------------------
# estimate_type remains "HR" for compatibility with downstream ratio-scale plotting code. These estimates are
# subdistribution hazard ratios, as identified by model_family = "finegray" and by the model formula.

results_long <- c(
  results_long,
  purrr::imap(vfd_cr_models, function(model, model_name) {
    tidy_results <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
    
    assemble_model_results(
      tidy_results = tidy_results,
      estimate_type = "HR",
      analysis = "28-day VFDs",
      model_spec = unname(exposure_labels[[model_name]]),
      model_family = "finegray",
      formula_str = model_formula_string(model),
      n_obs = attr(model, "n_patients"),
      adjustment = "adjusted"
    )
  })
)

# -----------------------------------------------------------------------------
# Primary adjusted 60-day Cox models
# -----------------------------------------------------------------------------

if (all(!purrr::map_lgl(cox_models, is.null))) {
  adjusted_cox_models <- cox_models
  adjusted_cox_labels <- cox_model_labels
  
  results_long <- c(
    results_long,
    purrr::imap(adjusted_cox_models, function(model, model_name) {
      tidy_results <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
      
      assemble_model_results(
        tidy_results = tidy_results,
        estimate_type = "HR",
        analysis = "Survival",
        model_spec = unname(adjusted_cox_labels[[model_name]]),
        model_family = "cox",
        formula_str = model_formula_string(model),
        n_obs = model$n,
        adjustment = "adjusted"
      )
    })
  )
}

# Validate the primary adjusted results before fitting minimally adjusted variants.
primary_results_qc <- bind_rows(results_long)

primary_results_checks <- list(
  primary_results_nonempty = nrow(primary_results_qc) > 0L,
  primary_mortality_results_present =
    is.null(mortality_models) ||
    any(
      primary_results_qc$analysis == "Mortality" &
        primary_results_qc$model_family == "logistic" &
        primary_results_qc$adjustment == "adjusted"
    ),
  primary_estimates_finite = all(is.finite(primary_results_qc$estimate)),
  primary_confidence_intervals_finite = all(
    is.finite(primary_results_qc$conf_low) & is.finite(primary_results_qc$conf_high)
  ),
  primary_standard_errors_finite = all(is.finite(primary_results_qc$std_error)),
  primary_statistics_finite = all(is.finite(primary_results_qc$statistic)),
  primary_p_values_valid = all(primary_results_qc$p_value >= 0 & primary_results_qc$p_value <= 1),
  primary_sample_sizes_positive = all(primary_results_qc$n_obs > 0L),
  mortality_results_use_60_day_outcome = primary_results_qc %>%
    filter(analysis == "Mortality") %>%
    summarise(valid = all(str_detect(formula, fixed("mortality_event_60")))) %>%
    pull(valid),
  continuous_results_use_hc3_rows = primary_results_qc %>%
    filter(model_family == "linear") %>%
    summarise(valid = n() == nrow(continuous_hc3_qc_combined)) %>%
    pull(valid),
  finegray_sample_sizes_are_patient_counts = primary_results_qc %>%
    filter(model_family == "finegray") %>%
    summarise(valid = all(n_obs == nrow(vfd_analysis_data))) %>%
    pull(valid),
  cox_sample_sizes_are_patient_counts = primary_results_qc %>%
    filter(model_family == "cox") %>%
    summarise(valid = all(n_obs == nrow(cox_analysis_data))) %>%
    pull(valid)
)


stopifnot(isTRUE(all(unlist(primary_results_checks))))

model_qc$unified_results <- list(
  primary_adjusted_summary = primary_results_qc %>%
    count(analysis, model_family, adjustment, n_obs, name = "n_result_rows")
)

model_qc$checks <- c(model_qc$checks, primary_results_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Primary adjusted long-format regression results validated.")

# -----------------------------------------------------------------------------
# Minimally adjusted exposure-outcome models
# -----------------------------------------------------------------------------
# Drop age, sex, and race while retaining illness severity. BMI remains included for driving-pressure-derived outcomes.
# These estimates are reported alongside the primary adjusted estimates, but the rendered tables, E-values, and evidence
# ratios above continue to use the primary adjusted models.

make_formula_unadj <- function(lhs, exposure, determinant) {
  cov <- model_covariates(determinant, adjusted = FALSE)
  rhs <- if (cov == "") exposure else paste(exposure, "+", cov)
  as.formula(paste(lhs, "~", rhs))
}

fit_vfd_finegray_minimal <- function(exposure_spec) {
  model_rhs <- paste(exposure_spec, "+", covariates_unadj)
  
  finegray_data <- vfd_analysis_data %>%
    mutate(vfd_status_f = factor(
      vfd_status, levels = c(0L, 1L, 2L),
      labels = c("censored", "extubation", "death")
    ))
  
  finegray_formula <- as.formula(
    paste("survival::Surv(vfd_time, vfd_status_f) ~ .patient_id +", model_rhs)
  )
  
  finegray_expanded <- survival::finegray(
    finegray_formula, data = finegray_data, etype = "extubation", id = .patient_id
  )
  
  model_formula <- as.formula(
    paste0("survival::Surv(fgstart, fgstop, fgstatus) ~ ", model_rhs, " + cluster(.patient_id)")
  )
  
  model <- survival::coxph(
    model_formula, data = finegray_expanded, weights = fgwt,
    robust = TRUE, x = TRUE, model = TRUE
  )
  
  attr(model, "n_patients") <- nrow(finegray_data)
  attr(model, "n_finegray_rows") <- nrow(finegray_expanded)
  attr(model, "n_liberation") <- sum(finegray_data$vfd_status == 1L)
  attr(model, "n_death") <- sum(finegray_data$vfd_status == 2L)
  attr(model, "n_censored") <- sum(finegray_data$vfd_status == 0L)
  
  model
}

if (!is.null(mortality_models)) {
  mortality_minimal_vars <- unique(c(
    mortality_outcome, exposure_vars, all.vars(as.formula(paste("~", covariates_unadj)))
  ))
  
  mortality_minimal_data <- cross_sectional %>%
    mutate(.source_row = row_number()) %>%
    select(.source_row, all_of(mortality_minimal_vars)) %>%
    drop_na(all_of(mortality_minimal_vars))
  
  mortality_models_minimal <- purrr::map(exposure_specs, function(exposure_spec) {
    glm(
      make_formula_unadj(mortality_outcome, exposure_spec, exposure_spec),
      data = mortality_minimal_data, family = binomial
    )
  })
  
  mortality_minimal_qc <- purrr::imap_dfr(mortality_models_minimal, function(model, model_name) {
    coefficient_table <- summary(model)$coefficients
    model_matrix <- model.matrix(model)
    outcome_values <- model.response(model.frame(model))
    n_parameters_excluding_intercept <- model$rank - 1L
    
    tibble(
      model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
      formula = model_formula_string(model), n_obs = nrow(model.frame(model)),
      n_events = sum(outcome_values == 1L), n_nonevents = sum(outcome_values == 0L),
      n_parameters_excluding_intercept = n_parameters_excluding_intercept,
      events_per_parameter = sum(outcome_values == 1L) / n_parameters_excluding_intercept,
      converged = isTRUE(model$converged), boundary = isTRUE(model$boundary),
      full_rank = qr(model_matrix)$rank == ncol(model_matrix),
      coefficients_finite = all(is.finite(coef(model))),
      standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
      n_extreme_fitted_probabilities = sum(fitted(model) < 1e-6 | fitted(model) > 1 - 1e-6)
    )
  })
  
  mortality_minimal_checks <- list(
    minimally_adjusted_mortality_models_fitted = length(mortality_models_minimal) == length(exposure_specs),
    minimally_adjusted_mortality_models_use_common_n = all(
      mortality_minimal_qc$n_obs == nrow(mortality_minimal_data)
    ),
    minimally_adjusted_mortality_models_have_expected_events = all(
      mortality_minimal_qc$n_events == sum(mortality_minimal_data[[mortality_outcome]] == 1L)
    ),
    minimally_adjusted_mortality_models_converged = all(mortality_minimal_qc$converged),
    no_minimally_adjusted_mortality_model_at_boundary = !any(mortality_minimal_qc$boundary),
    minimally_adjusted_mortality_model_matrices_full_rank = all(mortality_minimal_qc$full_rank),
    minimally_adjusted_mortality_coefficients_finite = all(mortality_minimal_qc$coefficients_finite),
    minimally_adjusted_mortality_standard_errors_finite = all(mortality_minimal_qc$standard_errors_finite),
    no_extreme_minimally_adjusted_mortality_probabilities = all(
      mortality_minimal_qc$n_extreme_fitted_probabilities == 0L
    )
  )
  
  
  stopifnot(isTRUE(all(unlist(mortality_minimal_checks))))
  
  results_long <- c(
    results_long,
    purrr::imap(mortality_models_minimal, function(model, model_name) {
      tidy_results <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
      
      assemble_model_results(
        tidy_results = tidy_results,
        estimate_type = "OR",
        analysis = "Mortality",
        model_spec = unname(exposure_labels[[model_name]]),
        model_family = "logistic",
        formula_str = model_formula_string(model),
        n_obs = nrow(model.frame(model)),
        adjustment = "minimally adjusted"
      )
    })
  )
  
  model_qc$unified_results$minimally_adjusted_mortality <- list(
    required_variables = mortality_minimal_vars,
    sample_n = nrow(mortality_minimal_data),
    event_n = sum(mortality_minimal_data[[mortality_outcome]] == 1L),
    model_diagnostics = mortality_minimal_qc,
    sparse_event_warning = mortality_minimal_qc %>%
      select(model_spec, exposure_label, events_per_parameter) %>%
      mutate(events_per_parameter_below_10 = events_per_parameter < 10)
  )
  
  model_qc$checks <- c(model_qc$checks, mortality_minimal_checks)
  model_qc$metadata$updated_at <- Sys.time()
  
  save_model_qc(model_qc, model_qc_path)
  
  stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))
  
  message("Minimally adjusted mortality results added to the unified results.")
}

continuous_models_minimal <- list()
continuous_hc3_minimal <- list()
continuous_minimal_qc <- list()

for (outcome_name in names(continuous_outcomes)) {
  outcome_var <- continuous_outcomes[[outcome_name]]$var
  outcome_label <- continuous_outcomes[[outcome_name]]$label
  minimal_covariates <- model_covariates(outcome_var, adjusted = FALSE)
  minimal_covariate_vars <- all.vars(as.formula(paste("~", minimal_covariates)))
  minimal_required_vars <- unique(c(outcome_var, exposure_vars, minimal_covariate_vars))
  
  continuous_minimal_data <- cross_sectional %>%
    mutate(.source_row = row_number()) %>%
    select(.source_row, all_of(minimal_required_vars)) %>%
    drop_na(all_of(minimal_required_vars))
  
  stopifnot(
    nrow(continuous_minimal_data) > 0L,
    n_distinct(continuous_minimal_data[[outcome_var]]) >= 2L
  )
  
  models_for_outcome <- purrr::map(exposure_specs, function(exposure_spec) {
    lm(
      make_formula_unadj(outcome_var, exposure_spec, outcome_var),
      data = continuous_minimal_data, x = TRUE, y = TRUE, model = TRUE
    )
  })
  
  hc3_for_outcome <- purrr::imap_dfr(models_for_outcome, function(model, model_name) {
    tidy_lm_hc3(model) %>%
      mutate(
        outcome_name = outcome_name, outcome_variable = outcome_var, outcome_label = outcome_label,
        model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
        n_obs = nrow(model.frame(model)), .before = 1
      )
  })
  
  model_qc_for_outcome <- purrr::imap_dfr(models_for_outcome, function(model, model_name) {
    model_matrix <- model.matrix(model)
    coefficient_table <- summary(model)$coefficients
    aliased_coefficients <- is.na(coef(model))
    
    tibble(
      outcome_name = outcome_name, outcome_variable = outcome_var, outcome_label = outcome_label,
      model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
      formula = model_formula_string(model), n_obs = nrow(model.frame(model)),
      residual_df = df.residual(model), model_matrix_columns = ncol(model_matrix),
      model_matrix_rank = qr(model_matrix)$rank, full_rank = qr(model_matrix)$rank == ncol(model_matrix),
      n_aliased_coefficients = sum(aliased_coefficients),
      coefficients_finite = all(is.finite(coef(model)[!aliased_coefficients])),
      standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
      fitted_values_finite = all(is.finite(fitted(model))),
      residuals_finite = all(is.finite(residuals(model)))
    )
  })
  
  same_rows_for_outcome <- all(purrr::map_lgl(
    models_for_outcome,
    ~ identical(rownames(model.frame(.x)), rownames(model.frame(models_for_outcome[[1]])))
  ))
  
  checks_for_outcome <- list(
    models_fitted = length(models_for_outcome) == length(exposure_specs),
    models_use_common_n = all(model_qc_for_outcome$n_obs == nrow(continuous_minimal_data)),
    models_use_same_rows = same_rows_for_outcome,
    model_matrices_full_rank = all(model_qc_for_outcome$full_rank),
    no_aliased_coefficients = all(model_qc_for_outcome$n_aliased_coefficients == 0L),
    coefficients_finite = all(model_qc_for_outcome$coefficients_finite),
    standard_errors_finite = all(model_qc_for_outcome$standard_errors_finite),
    fitted_values_finite = all(model_qc_for_outcome$fitted_values_finite),
    residuals_finite = all(model_qc_for_outcome$residuals_finite),
    residual_df_positive = all(model_qc_for_outcome$residual_df > 0L),
    hc3_estimates_finite = all(is.finite(hc3_for_outcome$estimate)),
    hc3_standard_errors_finite = all(is.finite(hc3_for_outcome$std.error)),
    hc3_confidence_intervals_finite = all(
      is.finite(hc3_for_outcome$conf.low) & is.finite(hc3_for_outcome$conf.high)
    ),
    hc3_p_values_valid = all(hc3_for_outcome$p.value >= 0 & hc3_for_outcome$p.value <= 1)
  )
  
  stopifnot(isTRUE(all(unlist(checks_for_outcome))))
  
  results_long <- c(
    results_long,
    purrr::imap(models_for_outcome, function(model, model_name) {
      tidy_results <- hc3_for_outcome %>%
        filter(model_spec == !!model_name) %>%
        select(term, estimate, conf.low, conf.high, std.error, statistic, p.value)
      
      assemble_model_results(
        tidy_results = tidy_results,
        estimate_type = "Beta",
        analysis = outcome_label,
        model_spec = unname(exposure_labels[[model_name]]),
        model_family = "linear",
        formula_str = model_formula_string(model),
        n_obs = nrow(model.frame(model)),
        adjustment = "minimally adjusted"
      )
    })
  )
  
  continuous_models_minimal[[outcome_name]] <- models_for_outcome
  continuous_hc3_minimal[[outcome_name]] <- hc3_for_outcome
  continuous_minimal_qc[[outcome_name]] <- list(
    sample_n = nrow(continuous_minimal_data),
    required_variables = minimal_required_vars,
    source_rows = continuous_minimal_data$.source_row,
    model_diagnostics = model_qc_for_outcome,
    checks = checks_for_outcome
  )
}

continuous_minimal_model_qc <- purrr::map_dfr(continuous_minimal_qc, ~ .x$model_diagnostics)
continuous_minimal_hc3_results <- bind_rows(continuous_hc3_minimal)

continuous_minimal_checks <- list(
  minimally_adjusted_continuous_outcomes_fitted = identical(
    names(continuous_models_minimal), names(continuous_outcomes)
  ),
  expected_minimally_adjusted_models_per_outcome = all(
    purrr::map_int(continuous_models_minimal, length) == length(exposure_specs)
  ),
  minimally_adjusted_continuous_model_matrices_full_rank = all(continuous_minimal_model_qc$full_rank),
  minimally_adjusted_continuous_coefficients_finite = all(continuous_minimal_model_qc$coefficients_finite),
  minimally_adjusted_continuous_hc3_results_finite = all(is.finite(c(
    continuous_minimal_hc3_results$estimate,
    continuous_minimal_hc3_results$std.error,
    continuous_minimal_hc3_results$conf.low,
    continuous_minimal_hc3_results$conf.high,
    continuous_minimal_hc3_results$p.value
  )))
)


stopifnot(isTRUE(all(unlist(continuous_minimal_checks))))

model_qc$unified_results$minimally_adjusted_continuous <- list(
  model_diagnostics = continuous_minimal_model_qc,
  HC3_coefficient_results = continuous_minimal_hc3_results,
  sample_summary = continuous_minimal_model_qc %>%
    count(outcome_name, outcome_label, n_obs, name = "n_models")
)

model_qc$checks <- c(model_qc$checks, continuous_minimal_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Minimally adjusted continuous-outcome results added to the unified results.")

vfd_cr_models_minimal <- purrr::map(exposure_specs, fit_vfd_finegray_minimal)

vfd_minimal_qc <- purrr::imap_dfr(vfd_cr_models_minimal, function(model, model_name) {
  coefficient_table <- summary(model)$coefficients
  model_matrix <- model$x
  
  tibble(
    model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
    formula = model_formula_string(model), n_patients = attr(model, "n_patients"),
    n_finegray_rows = attr(model, "n_finegray_rows"), n_liberation = attr(model, "n_liberation"),
    n_death = attr(model, "n_death"), n_censored = attr(model, "n_censored"),
    model_matrix_columns = ncol(model_matrix), model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    coefficients_finite = all(is.finite(coef(model))),
    robust_standard_errors_finite = all(is.finite(coefficient_table[, "robust se"]))
  )
})

vfd_minimal_checks <- list(
  minimally_adjusted_finegray_models_fitted = length(vfd_cr_models_minimal) == length(exposure_specs),
  minimally_adjusted_finegray_models_use_common_patient_n = all(
    vfd_minimal_qc$n_patients == nrow(vfd_analysis_data)
  ),
  minimally_adjusted_finegray_models_use_common_expanded_n = n_distinct(vfd_minimal_qc$n_finegray_rows) == 1L,
  minimally_adjusted_finegray_models_have_expected_liberations = all(
    vfd_minimal_qc$n_liberation == vfd_sample_qc$n_liberation
  ),
  minimally_adjusted_finegray_models_have_expected_deaths = all(
    vfd_minimal_qc$n_death == vfd_sample_qc$n_death
  ),
  minimally_adjusted_finegray_models_have_expected_censoring = all(
    vfd_minimal_qc$n_censored == vfd_sample_qc$n_status_0
  ),
  minimally_adjusted_finegray_model_matrices_full_rank = all(vfd_minimal_qc$full_rank),
  minimally_adjusted_finegray_coefficients_finite = all(vfd_minimal_qc$coefficients_finite),
  minimally_adjusted_finegray_robust_standard_errors_finite = all(
    vfd_minimal_qc$robust_standard_errors_finite
  )
)


stopifnot(isTRUE(all(unlist(vfd_minimal_checks))))

results_long <- c(
  results_long,
  purrr::imap(vfd_cr_models_minimal, function(model, model_name) {
    tidy_results <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
    
    assemble_model_results(
      tidy_results = tidy_results,
      estimate_type = "HR",
      analysis = "28-day VFDs",
      model_spec = unname(exposure_labels[[model_name]]),
      model_family = "finegray",
      formula_str = model_formula_string(model),
      n_obs = attr(model, "n_patients"),
      adjustment = "minimally adjusted"
    )
  })
)

model_qc$unified_results$minimally_adjusted_finegray <- list(
  model_diagnostics = vfd_minimal_qc,
  patient_n = nrow(vfd_analysis_data),
  liberation_n = vfd_sample_qc$n_liberation,
  death_n = vfd_sample_qc$n_death,
  censored_n = vfd_sample_qc$n_status_0
)

model_qc$checks <- c(model_qc$checks, vfd_minimal_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Minimally adjusted Fine-Gray results added to the unified results.")

if (length(cox_models) == 4L) {
  cox_unadj <- coxph(
    Surv(surv_time, event) ~ pbwpfvc + vtpbw + sf10 + sofa_total,
    data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
  )
  
  cox_pfvc_unadj <- coxph(
    Surv(surv_time, event) ~ pfvc + vtpbw + sf10 + sofa_total,
    data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
  )
  
  cox_abwpfvc_unadj <- coxph(
    Surv(surv_time, event) ~ abwpfvc + vtabw + sf10 + sofa_total,
    data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
  )
  
  cox_nchspfvc_unadj <- coxph(
    Surv(surv_time, event) ~ nchspfvc + vtnchs + sf10 + sofa_total,
    data = cox_analysis_data, ties = "efron", x = TRUE, model = TRUE, y = TRUE
  )
  
  cox_models_minimal <- list(
    vtpbw_pbwpfvc = cox_unadj,
    vtnchs_nchspfvc = cox_nchspfvc_unadj,
    vtabw_abwpfvc = cox_abwpfvc_unadj,
    vtpbw_pfvc = cox_pfvc_unadj
  )
  
  cox_minimal_qc <- purrr::imap_dfr(cox_models_minimal, function(model, model_name) {
    coefficient_table <- summary(model)$coefficients
    model_matrix <- model$x
    n_patients_model <- model$n
    n_events_model <- model$nevent
    n_coefficients <- length(coef(model))
    
    tibble(
      model_spec = model_name, model_label = unname(cox_model_labels[[model_name]]),
      formula = model_formula_string(model), n_obs = n_patients_model,
      n_events = n_events_model, n_censored = n_patients_model - n_events_model,
      n_response_rows = nrow(model$y), n_coefficients = n_coefficients,
      events_per_coefficient = n_events_model / n_coefficients, iterations = model$iter,
      model_matrix_columns = ncol(model_matrix), model_matrix_rank = qr(model_matrix)$rank,
      full_rank = qr(model_matrix)$rank == ncol(model_matrix),
      coefficients_finite = all(is.finite(coef(model))),
      standard_errors_finite = all(is.finite(coefficient_table[, "se(coef)"])),
      maximum_absolute_coefficient = max(abs(coef(model))),
      maximum_standard_error = max(coefficient_table[, "se(coef)"]),
      log_likelihood = as.numeric(logLik(model)), AIC = AIC(model)
    )
  })
  
  cox_minimal_reference_rows <- rownames(model.frame(cox_models_minimal[[1L]]))
  
  cox_minimal_same_rows <- all(
    purrr::map_lgl(
      cox_models_minimal,
      ~ identical(rownames(model.frame(.x)), cox_minimal_reference_rows)
    )
  )
  
  cox_minimal_checks <- list(
    all_minimally_adjusted_cox_models_fitted = length(cox_models_minimal) == 4L,
    minimally_adjusted_cox_models_use_common_n = all(cox_minimal_qc$n_obs == nrow(cox_analysis_data)),
    minimally_adjusted_cox_response_rows_match_common_n = all(
      cox_minimal_qc$n_response_rows == nrow(cox_analysis_data)
    ),
    minimally_adjusted_cox_models_use_same_rows = cox_minimal_same_rows,
    minimally_adjusted_cox_models_have_expected_events = all(
      cox_minimal_qc$n_events == cox_sample_qc$n_events
    ),
    minimally_adjusted_cox_models_have_expected_censoring = all(
      cox_minimal_qc$n_censored == cox_sample_qc$n_censored
    ),
    minimally_adjusted_cox_model_matrices_full_rank = all(cox_minimal_qc$full_rank),
    minimally_adjusted_cox_coefficients_finite = all(cox_minimal_qc$coefficients_finite),
    minimally_adjusted_cox_standard_errors_finite = all(cox_minimal_qc$standard_errors_finite),
    minimally_adjusted_cox_log_likelihoods_finite = all(is.finite(cox_minimal_qc$log_likelihood)),
    minimally_adjusted_cox_aic_values_finite = all(is.finite(cox_minimal_qc$AIC))
  )
  
  
  stopifnot(isTRUE(all(unlist(cox_minimal_checks))))
  
  results_long <- c(
    results_long,
    purrr::imap(cox_models_minimal, function(model, model_name) {
      tidy_results <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
      
      assemble_model_results(
        tidy_results = tidy_results,
        estimate_type = "HR",
        analysis = "Survival",
        model_spec = unname(cox_model_labels[[model_name]]),
        model_family = "cox",
        formula_str = model_formula_string(model),
        n_obs = model$n,
        adjustment = "minimally adjusted"
      )
    })
  )
  
  model_qc$unified_results$minimally_adjusted_cox <- list(
    sample_n = nrow(cox_analysis_data),
    event_n = cox_sample_qc$n_events,
    censored_n = cox_sample_qc$n_censored,
    model_diagnostics = cox_minimal_qc,
    sparse_event_warning = cox_minimal_qc %>%
      select(model_spec, model_label, events_per_coefficient) %>%
      mutate(events_per_coefficient_below_10 = events_per_coefficient < 10)
  )
  
  model_qc$checks <- c(model_qc$checks, cox_minimal_checks)
  model_qc$metadata$updated_at <- Sys.time()
  
  save_model_qc(model_qc, model_qc_path)
  
  stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))
  
  message("Minimally adjusted Cox results added to the unified results.")
}

# -----------------------------------------------------------------------------
# Demographic-bias models
# -----------------------------------------------------------------------------
# Linear models use the HC3 results reported in the rendered table. The 60-day mortality model uses model-based logistic
# odds ratios. These analyses retain model_spec = "Demographics" for compatibility with downstream scripts.

results_long <- c(
  results_long,
  purrr::imap(demo_models, function(model_info, outcome_name) {
    tidy_results <- demo_model_results %>%
      filter(outcome == !!outcome_name) %>%
      select(term, estimate, conf.low, conf.high, std.error, statistic, p.value)
    
    assemble_model_results(
      tidy_results = tidy_results,
      estimate_type = if (identical(model_info$type, "OR")) "OR" else "Beta",
      analysis = paste0("Demo bias: ", outcome_name),
      model_spec = "Demographics",
      model_family = model_info$family,
      formula_str = model_info$formula,
      n_obs = nrow(model.frame(model_info$model)),
      adjustment = "adjusted"
    )
  })
)

# -----------------------------------------------------------------------------
# Predicted FVC versus weight-method estimates in the broad cohort
# -----------------------------------------------------------------------------

results_long <- c(
  results_long,
  purrr::imap(pfvc_vs_weight_models, function(model, method_name) {
    tidy_results <- pfvc_vs_weight_hc3_results %>%
      filter(method == !!method_name) %>%
      select(term, estimate, conf.low, conf.high, std.error, statistic, p.value)
    
    assemble_model_results(
      tidy_results = tidy_results,
      estimate_type = "Beta",
      analysis = "PFVC vs weight method",
      model_spec = unname(pfvc_vs_weight_labels[[method_name]]),
      model_family = "linear",
      formula_str = model_formula_string(model),
      n_obs = nobs(model),
      adjustment = "adjusted"
    )
  })
)

# Drop intercepts and stamp the site name so the table remains self-describing after cross-site aggregation.
regression_results_long <- bind_rows(results_long) %>%
  filter(term != "(Intercept)") %>%
  mutate(site = site_name, .before = 1)

unified_results_required_columns <- c(
  "site", "term", "estimate", "conf_low", "conf_high", "std_error", "statistic", "p_value",
  "estimate_type", "analysis", "model_spec", "model_family", "adjustment", "formula", "n_obs"
)

missing_unified_results_columns <- setdiff(unified_results_required_columns, names(regression_results_long))

unified_results_group_qc <- regression_results_long %>%
  count(analysis, model_family, adjustment, n_obs, name = "n_result_rows") %>%
  arrange(analysis, adjustment, n_obs)

unified_results_duplicate_qc <- regression_results_long %>%
  count(site, analysis, model_spec, model_family, adjustment, term, name = "n_rows") %>%
  filter(n_rows > 1L)

unified_results_checks <- list(
  unified_results_required_columns_present = length(missing_unified_results_columns) == 0L,
  unified_results_nonempty = nrow(regression_results_long) > 0L,
  unified_results_have_no_intercepts = !any(regression_results_long$term == "(Intercept)"),
  unified_results_site_complete = all(
    !is.na(regression_results_long$site) & regression_results_long$site == site_name
  ),
  unified_results_estimates_finite = all(is.finite(regression_results_long$estimate)),
  unified_results_confidence_intervals_finite = all(
    is.finite(regression_results_long$conf_low) & is.finite(regression_results_long$conf_high)
  ),
  unified_results_standard_errors_finite = all(is.finite(regression_results_long$std_error)),
  unified_results_statistics_finite = all(is.finite(regression_results_long$statistic)),
  unified_results_p_values_valid = all(
    regression_results_long$p_value >= 0 & regression_results_long$p_value <= 1
  ),
  unified_results_sample_sizes_positive = all(regression_results_long$n_obs > 0L),
  unified_results_estimate_types_valid = all(regression_results_long$estimate_type %in% c("OR", "HR", "Beta")),
  unified_results_adjustment_labels_valid = all(
    regression_results_long$adjustment %in% c("adjusted", "minimally adjusted")
  ),
  unified_results_have_no_duplicate_keys = nrow(unified_results_duplicate_qc) == 0L,
  unified_mortality_formulas_use_60_day_outcome = regression_results_long %>%
    filter(analysis == "Mortality") %>%
    summarise(valid = all(str_detect(formula, fixed("mortality_event_60")))) %>%
    pull(valid),
  unified_finegray_sample_sizes_are_patient_counts = regression_results_long %>%
    filter(model_family == "finegray") %>%
    summarise(valid = all(n_obs == nrow(vfd_analysis_data))) %>%
    pull(valid),
  unified_cox_sample_sizes_are_patient_counts = regression_results_long %>%
    filter(model_family == "cox") %>%
    summarise(valid = all(n_obs == nrow(cox_analysis_data))) %>%
    pull(valid)
)


stopifnot(isTRUE(all(unlist(unified_results_checks))))

regression_results_csv_path <- file.path(table_dir,paste0("regression_results_long_", site_name, ".csv"))
regression_results_parquet_path <- file.path(table_dir, paste0("regression_results_long_", site_name, ".parquet"))

write_csv(regression_results_long, regression_results_csv_path)
write_parquet(regression_results_long, regression_results_parquet_path)

regression_results_csv_exists <- file.exists(regression_results_csv_path)
regression_results_parquet_exists <- file.exists(regression_results_parquet_path)

regression_results_csv_size <- if (regression_results_csv_exists) {
  file.info(regression_results_csv_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

regression_results_parquet_size <- if (regression_results_parquet_exists) {
  file.info(regression_results_parquet_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

unified_results_output_checks <- list(
  unified_results_csv_exists = regression_results_csv_exists,
  unified_results_csv_nonempty = regression_results_csv_exists &&
    is.finite(regression_results_csv_size) &&
    regression_results_csv_size > 0,
  unified_results_parquet_exists = regression_results_parquet_exists,
  unified_results_parquet_nonempty = regression_results_parquet_exists &&
    is.finite(regression_results_parquet_size) &&
    regression_results_parquet_size > 0
)

model_qc$unified_results$final_artifact <- list(
  row_count = nrow(regression_results_long),
  analysis_count = n_distinct(regression_results_long$analysis),
  group_summary = unified_results_group_qc,
  duplicate_keys = unified_results_duplicate_qc,
  csv_path = regression_results_csv_path,
  csv_size_bytes = regression_results_csv_size,
  parquet_path = regression_results_parquet_path,
  parquet_size_bytes = regression_results_parquet_size,
  schema = names(regression_results_long)
)

model_qc$checks <- c(model_qc$checks, unified_results_checks, unified_results_output_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message(
  "Unified regression results table: ", nrow(regression_results_long),
  " rows across ", n_distinct(regression_results_long$analysis),
  " analyses; saved to ", regression_results_csv_path
)

# =============================================================================
# 4g2. Sensitivity analyses: E-values and nonlinear age adjustment
# =============================================================================

# PBW and PFVC are constructed from demographic inputs, including age, sex, height, and race. This deterministic
# construction does not establish that these variables constitute a complete causal adjustment set or that height is
# necessarily a mediator rather than a potential common cause. The appropriate adjustment set depends on the target
# estimand and assumed causal structure.
#
# This section evaluates two forms of sensitivity:
#
#   1. Age functional form: primary models containing linear age are compared with otherwise similar models containing
#      a four-degree-of-freedom natural spline for age. This assesses sensitivity to nonlinear age adjustment but does
#      not address other forms of residual or unmeasured confounding.
#
#   2. E-values: ratio-scale estimates are summarized using point-estimate and confidence-limit E-values. An E-value
#      describes the minimum strength of association that an unmeasured confounder would need with both the exposure and
#      outcome, conditional on measured covariates, to explain away the observed association under the assumptions of
#      the E-value framework.
#
# E-values do not address selection bias, measurement error, model misspecification, positivity violations, informative
# censoring, or exposure misclassification. Odds ratios and hazard ratios are converted using the EValue package's
# common-outcome approximations. Applying the hazard-ratio E-value approximation to a Fine-Gray subdistribution hazard
# ratio is exploratory because a subdistribution hazard ratio is not a conventional risk ratio or cause-specific hazard
# ratio.
#
# Exposure estimates are rescaled to a one-SD increase using the patient-level analysis sample shared by the linear-age
# and spline-age versions of each model. This provides a common contrast within each analysis but does not make effect
# sizes directly comparable across different exposures or outcomes.

EXPOSURE_TERMS <- c(
  "vtpbw", "pbwpfvc",
  "vtnchs", "nchspfvc",
  "vtabw", "abwpfvc",
  "vtpfvc", "pfvc"
)

exposure_term_labels_ev <- c(
  vtpbw = "VT/McLaren PBW",
  pbwpfvc = "McLaren PBW/PFVC",
  vtnchs = "VT/NCHS weight",
  nchspfvc = "NCHS weight/PFVC",
  vtabw = "VT/ABW",
  abwpfvc = "ABW/PFVC",
  vtpfvc = "VT/PFVC",
  pfvc = "PFVC"
)

AGE_SPLINE_DF <- 4L

# This replaces linear age_at_admission with a four-degree-of-freedom natural spline while retaining the remaining
# primary adjustment covariates.
spline_covars <- paste0(
  "race_category + ns(age_at_admission, df = ", AGE_SPLINE_DF,
  ") + sex_category + sofa_total + sf10"
)

# Calculate a per-SD ratio estimate and 95% confidence interval from a fitted model. beta and its standard error are on
# the log-ratio scale; multiplying by the exposure SD gives the per-SD log ratio before exponentiation.
persd_ratio <- function(model, term, exposure_sd) {
  beta <- coef(model)[[term]]
  se <- sqrt(diag(vcov(model)))[[term]]
  
  list(
    estimate = exp(beta * exposure_sd),
    conf_low = exp((beta - 1.96 * se) * exposure_sd),
    conf_high = exp((beta + 1.96 * se) * exposure_sd)
  )
}

# Calculate per-SD estimates and confidence intervals for every focal exposure term in a model. The SD is the marginal,
# patient-level SD of the exposure in sd_data, so the linear-age and spline-age fits use an identical contrast.
persd_terms <- function(model, sd_data, terms = EXPOSURE_TERMS) {
  present <- intersect(terms, names(coef(model)))
  
  map_dfr(present, function(tm) {
    exposure_sd <- sd(sd_data[[tm]], na.rm = TRUE)
    pr <- persd_ratio(model, tm, exposure_sd)
    
    tibble(
      term = tm,
      exposure_sd = exposure_sd,
      estimate = pr$estimate,
      conf_low = pr$conf_low,
      conf_high = pr$conf_high
    )
  })
}

# Calculate point-estimate and confidence-limit E-values for one per-SD ratio estimate. OR and HR calculations use
# rare = FALSE for common outcomes. SHR calculations use the EValue package's HR conversion as an approximation.
evalue_for <- function(est, lo, hi, type) {
  if (
    length(est) != 1L || length(lo) != 1L || length(hi) != 1L ||
    !all(is.finite(c(est, lo, hi))) || est <= 0 || lo <= 0 || hi <= 0 || lo > hi
  ) {
    stop("E-value inputs must be finite positive scalars satisfying conf_low <= conf_high.")
  }
  
  ev_object <- suppressMessages(
    suppressWarnings(
      switch(
        type,
        OR = EValue::evalues.OR(est, lo, hi, rare = FALSE),
        HR = EValue::evalues.HR(est, lo, hi, rare = FALSE),
        SHR = EValue::evalues.HR(est, lo, hi, rare = FALSE),
        stop("Unsupported estimate type for E-value calculation: ", type)
      )
    )
  )
  
  ev_matrix <- as.matrix(ev_object)
  
  if (nrow(ev_matrix) < 1L || ncol(ev_matrix) < 1L) {
    stop("The EValue package returned an empty result.")
  }
  
  row_labels <- tolower(trimws(rownames(ev_matrix)))
  column_labels <- tolower(trimws(colnames(ev_matrix)))
  
  evalue_row <- grep("e.?value", row_labels)
  evalue_row <- if (length(evalue_row) > 0L) evalue_row[[1]] else nrow(ev_matrix)
  
  point_column <- grep("point", column_labels)
  lower_column <- grep("lower", column_labels)
  upper_column <- grep("upper", column_labels)
  
  point_column <- if (length(point_column) > 0L) point_column[[1]] else 1L
  lower_column <- if (length(lower_column) > 0L) lower_column[[1]] else min(2L, ncol(ev_matrix))
  upper_column <- if (length(upper_column) > 0L) upper_column[[1]] else min(3L, ncol(ev_matrix))
  
  evalue_point <- as.numeric(ev_matrix[evalue_row, point_column])
  
  evalue_ci <- if (lo <= 1 && hi >= 1) {
    1
  } else if (lo > 1) {
    as.numeric(ev_matrix[evalue_row, lower_column])
  } else {
    as.numeric(ev_matrix[evalue_row, upper_column])
  }
  
  if (!is.finite(evalue_ci)) {
    ci_candidates <- as.numeric(ev_matrix[evalue_row, unique(c(lower_column, upper_column))])
    ci_candidates <- ci_candidates[is.finite(ci_candidates)]
    evalue_ci <- if (length(ci_candidates) > 0L) ci_candidates[[1]] else NA_real_
  }
  
  if (!is.finite(evalue_point) || evalue_point < 1) {
    stop("The EValue package returned an invalid point-estimate E-value.")
  }
  
  if (!is.finite(evalue_ci) || evalue_ci < 1) {
    stop("The EValue package returned an invalid confidence-limit E-value.")
  }
  
  list(evalue_point = evalue_point, evalue_ci = evalue_ci)
}

# Return one row per focal exposure with linear-age and spline-age estimates side by side, plus the E-value for the
# linear-age estimate. Both models use the same sd_data, so the per-SD contrast is identical.
residual_conf_rows <- function(model_lin, model_spl, model_spec, analysis, type, sd_data) {
  lin <- persd_terms(model_lin, sd_data)
  spl <- persd_terms(model_spl, sd_data)
  
  if (nrow(lin) == 0L) {
    return(tibble())
  }
  
  if (!setequal(lin$term, spl$term)) {
    stop("Linear-age and spline-age models do not contain the same focal exposure terms for: ", model_spec)
  }
  
  if (any(!is.finite(lin$exposure_sd)) || any(lin$exposure_sd <= 0)) {
    stop("A focal exposure has an invalid SD for: ", model_spec)
  }
  
  spline_match <- match(lin$term, spl$term)
  
  ev <- purrr::pmap_dfr(
    list(lin$estimate, lin$conf_low, lin$conf_high),
    function(estimate, conf_low, conf_high) {
      ev_result <- evalue_for(est = estimate, lo = conf_low, hi = conf_high, type = type)
      
      tibble(
        evalue_point = as.numeric(ev_result$evalue_point),
        evalue_ci = as.numeric(ev_result$evalue_ci)
      )
    }
  )
  
  stopifnot(
    nrow(ev) == nrow(lin),
    all(is.finite(ev$evalue_point)),
    all(is.finite(ev$evalue_ci))
  )
  
  tibble(
    analysis = analysis,
    model_spec = model_spec,
    term = lin$term,
    term_label = unname(exposure_term_labels_ev[lin$term]),
    estimate_type = type,
    exposure_sd = lin$exposure_sd,
    est_linage = lin$estimate,
    lo_linage = lin$conf_low,
    hi_linage = lin$conf_high,
    est_splineage = spl$estimate[spline_match],
    lo_splineage = spl$conf_low[spline_match],
    hi_splineage = spl$conf_high[spline_match],
    evalue_point = ev$evalue_point,
    evalue_ci = ev$evalue_ci
  )
}

# =============================================================================
# Spline-age refits on the validated common analysis samples
# =============================================================================

# Each spline-age model uses the same patients as its corresponding linear-age model. Only the functional form of age
# changes, allowing differences between estimates to be attributed to age specification rather than sample composition.

mortality_models_spline <- if (!is.null(mortality_models)) {
  purrr::map(exposure_specs, function(exposure_spec) {
    formula_str <- paste("mortality_event_60 ~", exposure_spec, "+", spline_covars)
    
    glm(
      as.formula(formula_str), data = mortality_analysis_data,
      family = binomial, x = TRUE, y = TRUE, model = TRUE
    )
  })
} else {
  NULL
}

cox_model_spline <- if (!is.null(cox_model)) {
  coxph(
    Surv(surv_time, event) ~ pbwpfvc + vtpbw + ns(age_at_admission, df = AGE_SPLINE_DF) +
      sex_category + race_category + sf10 + sofa_total,
    data = cox_analysis_data, ties = "efron", x = TRUE, y = TRUE, model = TRUE
  )
} else {
  NULL
}

cox_model_pfvc_spline <- if (!is.null(cox_model_pfvc)) {
  coxph(
    Surv(surv_time, event) ~ pfvc + vtpbw + ns(age_at_admission, df = AGE_SPLINE_DF) +
      sex_category + race_category + sf10 + sofa_total,
    data = cox_analysis_data, ties = "efron", x = TRUE, y = TRUE, model = TRUE
  )
} else {
  NULL
}

cox_model_nchspfvc_spline <- if (!is.null(cox_model_nchspfvc)) {
  coxph(
    Surv(surv_time, event) ~ nchspfvc + vtnchs + ns(age_at_admission, df = AGE_SPLINE_DF) +
      sex_category + race_category + sf10 + sofa_total,
    data = cox_analysis_data, ties = "efron", x = TRUE, y = TRUE, model = TRUE
  )
} else {
  NULL
}

cox_model_abwpfvc_spline <- if (!is.null(cox_model_abwpfvc)) {
  coxph(
    Surv(surv_time, event) ~ abwpfvc + vtabw + ns(age_at_admission, df = AGE_SPLINE_DF) +
      sex_category + race_category + sf10 + sofa_total,
    data = cox_analysis_data, ties = "efron", x = TRUE, y = TRUE, model = TRUE
  )
} else {
  NULL
}

fit_vfd_finegray_spline <- function(exposure_spec) {
  age_spline_matrix <- splines::ns(vfd_analysis_data$age_at_admission, df = AGE_SPLINE_DF)
  age_spline_data <- as_tibble(age_spline_matrix, .name_repair = "minimal")
  names(age_spline_data) <- paste0("age_spline_", seq_len(ncol(age_spline_data)))
  
  finegray_data <- bind_cols(vfd_analysis_data, age_spline_data) %>%
    mutate(vfd_status_f = factor(
      vfd_status, levels = c(0L, 1L, 2L),
      labels = c("censored", "extubation", "death")
    ))
  
  spline_basis_vars <- names(age_spline_data)
  
  model_rhs <- paste(
    exposure_spec, "+ race_category +", paste(spline_basis_vars, collapse = " + "),
    "+ sex_category + sofa_total + sf10"
  )
  
  finegray_formula <- as.formula(
    paste("survival::Surv(vfd_time, vfd_status_f) ~ .patient_id +", model_rhs)
  )
  
  finegray_expanded <- survival::finegray(
    finegray_formula, data = finegray_data, etype = "extubation", id = .patient_id
  )
  
  model_formula <- as.formula(
    paste0("survival::Surv(fgstart, fgstop, fgstatus) ~ ", model_rhs, " + cluster(.patient_id)")
  )
  
  model <- survival::coxph(
    model_formula, data = finegray_expanded, weights = fgwt,
    robust = TRUE, x = TRUE, y = TRUE, model = TRUE
  )
  
  attr(model, "n_patients") <- nrow(finegray_data)
  attr(model, "n_finegray_rows") <- nrow(finegray_expanded)
  attr(model, "n_liberation") <- sum(finegray_data$vfd_status == 1L)
  attr(model, "n_death") <- sum(finegray_data$vfd_status == 2L)
  attr(model, "n_censored") <- sum(finegray_data$vfd_status == 0L)
  attr(model, "age_spline_df") <- AGE_SPLINE_DF
  attr(model, "age_spline_vars") <- spline_basis_vars
  
  model
}

vfd_cr_models_spline <- purrr::map(exposure_specs, fit_vfd_finegray_spline)

# =============================================================================
# Model QC: spline-age sensitivity models
# =============================================================================

mortality_spline_qc <- if (!is.null(mortality_models_spline)) {
  purrr::imap_dfr(mortality_models_spline, function(model, model_name) {
    coefficient_table <- summary(model)$coefficients
    model_matrix <- model.matrix(model)
    outcome_values <- model.response(model.frame(model))
    
    tibble(
      model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
      n_obs = nrow(model.frame(model)), n_events = sum(outcome_values == 1L),
      n_nonevents = sum(outcome_values == 0L), converged = isTRUE(model$converged),
      boundary = isTRUE(model$boundary), model_matrix_columns = ncol(model_matrix),
      model_matrix_rank = qr(model_matrix)$rank, full_rank = qr(model_matrix)$rank == ncol(model_matrix),
      coefficients_finite = all(is.finite(coef(model))),
      standard_errors_finite = all(is.finite(coefficient_table[, "Std. Error"])),
      n_extreme_fitted_probabilities = sum(fitted(model) < 1e-6 | fitted(model) > 1 - 1e-6)
    )
  })
} else {
  tibble()
}

cox_spline_models <- list(
  vtpbw_pbwpfvc = cox_model_spline,
  vtnchs_nchspfvc = cox_model_nchspfvc_spline,
  vtabw_abwpfvc = cox_model_abwpfvc_spline,
  vtpbw_pfvc = cox_model_pfvc_spline
) %>%
  purrr::compact()

cox_spline_qc <- purrr::imap_dfr(cox_spline_models, function(model, model_name) {
  coefficient_table <- summary(model)$coefficients
  model_matrix <- model$x
  
  tibble(
    model_spec = model_name, model_label = unname(cox_model_labels[[model_name]]),
    n_obs = model$n, n_events = model$nevent, n_censored = model$n - model$nevent,
    model_matrix_columns = ncol(model_matrix), model_matrix_rank = qr(model_matrix)$rank,
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    coefficients_finite = all(is.finite(coef(model))),
    standard_errors_finite = all(is.finite(coefficient_table[, "se(coef)"]))
  )
})

vfd_spline_qc <- purrr::imap_dfr(vfd_cr_models_spline, function(model, model_name) {
  coefficient_table <- summary(model)$coefficients
  model_matrix <- model$x
  
  tibble(
    model_spec = model_name, exposure_label = unname(exposure_labels[[model_name]]),
    n_patients = attr(model, "n_patients"), n_finegray_rows = attr(model, "n_finegray_rows"),
    n_liberation = attr(model, "n_liberation"), n_death = attr(model, "n_death"),
    n_censored = attr(model, "n_censored"), model_matrix_columns = ncol(model_matrix),
    model_matrix_rank = qr(model_matrix)$rank, full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    coefficients_finite = all(is.finite(coef(model))),
    robust_standard_errors_finite = all(is.finite(coefficient_table[, "robust se"]))
  )
})

spline_age_model_checks <- list(
  mortality_spline_models_use_primary_sample = is.null(mortality_models_spline) ||
    all(mortality_spline_qc$n_obs == nrow(mortality_analysis_data)),
  mortality_spline_models_have_expected_events = is.null(mortality_models_spline) ||
    all(mortality_spline_qc$n_events == mortality_sample_qc$n_events),
  mortality_spline_models_converged = is.null(mortality_models_spline) || all(mortality_spline_qc$converged),
  no_mortality_spline_model_at_boundary = is.null(mortality_models_spline) || !any(mortality_spline_qc$boundary),
  mortality_spline_model_matrices_full_rank = is.null(mortality_models_spline) || all(mortality_spline_qc$full_rank),
  mortality_spline_coefficients_finite = is.null(mortality_models_spline) ||
    all(mortality_spline_qc$coefficients_finite),
  mortality_spline_standard_errors_finite = is.null(mortality_models_spline) ||
    all(mortality_spline_qc$standard_errors_finite),
  no_extreme_mortality_spline_probabilities = is.null(mortality_models_spline) ||
    all(mortality_spline_qc$n_extreme_fitted_probabilities == 0L),
  all_cox_spline_models_fitted = length(cox_spline_models) == 4L,
  cox_spline_models_use_primary_sample = all(cox_spline_qc$n_obs == nrow(cox_analysis_data)),
  cox_spline_models_have_expected_events = all(cox_spline_qc$n_events == cox_sample_qc$n_events),
  cox_spline_model_matrices_full_rank = all(cox_spline_qc$full_rank),
  cox_spline_coefficients_finite = all(cox_spline_qc$coefficients_finite),
  cox_spline_standard_errors_finite = all(cox_spline_qc$standard_errors_finite),
  expected_finegray_spline_models_fitted = length(vfd_cr_models_spline) == length(exposure_specs),
  finegray_spline_models_use_primary_patient_sample = all(
    vfd_spline_qc$n_patients == nrow(vfd_analysis_data)
  ),
  finegray_spline_models_use_common_expanded_rows = n_distinct(vfd_spline_qc$n_finegray_rows) == 1L,
  finegray_spline_models_have_expected_liberations = all(
    vfd_spline_qc$n_liberation == vfd_sample_qc$n_liberation
  ),
  finegray_spline_models_have_expected_deaths = all(vfd_spline_qc$n_death == vfd_sample_qc$n_death),
  finegray_spline_models_have_expected_censoring = all(vfd_spline_qc$n_censored == vfd_sample_qc$n_status_0),
  finegray_spline_model_matrices_full_rank = all(vfd_spline_qc$full_rank),
  finegray_spline_coefficients_finite = all(vfd_spline_qc$coefficients_finite),
  finegray_spline_robust_standard_errors_finite = all(vfd_spline_qc$robust_standard_errors_finite)
)


stopifnot(isTRUE(all(unlist(spline_age_model_checks))))

model_qc$residual_confounding <- list(
  age_spline_df = AGE_SPLINE_DF,
  mortality_spline_model_diagnostics = mortality_spline_qc,
  cox_spline_model_diagnostics = cox_spline_qc,
  finegray_spline_model_diagnostics = vfd_spline_qc,
  interpretation = paste0(
    "Spline-age models use the same patient samples as the corresponding linear-age models. ",
    "Only the functional form of age differs."
  )
)

model_qc$checks <- c(model_qc$checks, spline_age_model_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Spline-age sensitivity-model QC added to: ", model_qc_path)

# =============================================================================
# Assemble age-functional-form and E-value sensitivity results
# =============================================================================

residual_list <- list()

# Sixty-day mortality odds ratios.
if (!is.null(mortality_models) && !is.null(mortality_models_spline)) {
  residual_list <- c(
    residual_list,
    purrr::imap(mortality_models, function(model, model_name) {
      residual_conf_rows(
        model_lin = model, model_spl = mortality_models_spline[[model_name]],
        model_spec = unname(exposure_labels[[model_name]]), analysis = "Mortality (60-day)",
        type = "OR", sd_data = mortality_analysis_data
      )
    })
  )
}

# Sixty-day death hazard ratios.
if (length(cox_models) == 4L && length(cox_spline_models) == 4L) {
  residual_list <- c(
    residual_list,
    purrr::imap(cox_models, function(model, model_name) {
      residual_conf_rows(
        model_lin = model,
        model_spl = cox_spline_models[[model_name]],
        model_spec = unname(cox_model_labels[[model_name]]),
        analysis = "Survival (60-day)",
        type = "HR",
        sd_data = cox_analysis_data
      )
    })
  )
}

# Twenty-eight-day successful-liberation subdistribution hazard ratios. E-values use the EValue package's HR conversion
# as an exploratory approximation; patient-level SDs are calculated before Fine-Gray row expansion.
residual_list <- c(
  residual_list,
  purrr::imap(vfd_cr_models, function(model, model_name) {
    residual_conf_rows(
      model_lin = model, model_spl = vfd_cr_models_spline[[model_name]],
      model_spec = unname(exposure_labels[[model_name]]), analysis = "28-day VFDs (liberation)",
      type = "SHR", sd_data = vfd_analysis_data
    )
  })
)

residual_confounding <- bind_rows(residual_list) %>%
  mutate(
    analysis = factor(
      analysis,
      levels = c("Mortality (60-day)", "Survival (60-day)", "28-day VFDs (liberation)")
    )
  ) %>%
  arrange(analysis, model_spec, term) %>%
  mutate(analysis = as.character(analysis), site = site_name, .before = 1)

residual_confounding_qc <- residual_confounding %>%
  group_by(analysis, estimate_type) %>%
  summarise(
    n_rows = n(), n_models = n_distinct(model_spec), n_terms = n_distinct(term),
    all_exposure_sds_positive = all(is.finite(exposure_sd) & exposure_sd > 0),
    all_linear_estimates_finite = all(is.finite(est_linage) & is.finite(lo_linage) & is.finite(hi_linage)),
    all_spline_estimates_finite = all(is.finite(est_splineage) & is.finite(lo_splineage) & is.finite(hi_splineage)),
    all_evalue_points_valid = all(is.finite(evalue_point) & evalue_point >= 1),
    all_evalue_ci_values_valid = all(is.finite(evalue_ci) & evalue_ci >= 1),
    .groups = "drop"
  )

residual_confounding_duplicate_qc <- residual_confounding %>%
  count(site, analysis, model_spec, term, name = "n_rows") %>%
  filter(n_rows > 1L)

residual_confounding_checks <- list(
  residual_confounding_results_nonempty = nrow(residual_confounding) > 0L,
  residual_confounding_analysis_labels_complete = !any(is.na(residual_confounding$analysis)),
  residual_confounding_term_labels_complete = !any(is.na(residual_confounding$term_label)),
  residual_confounding_estimate_types_valid = all(
    residual_confounding$estimate_type %in% c("OR", "HR", "SHR")
  ),
  residual_confounding_exposure_sds_positive = all(residual_confounding_qc$all_exposure_sds_positive),
  residual_confounding_linear_estimates_finite = all(residual_confounding_qc$all_linear_estimates_finite),
  residual_confounding_spline_estimates_finite = all(residual_confounding_qc$all_spline_estimates_finite),
  residual_confounding_evalue_points_valid = all(residual_confounding_qc$all_evalue_points_valid),
  residual_confounding_evalue_ci_values_valid = all(residual_confounding_qc$all_evalue_ci_values_valid),
  residual_confounding_has_no_duplicate_keys = nrow(residual_confounding_duplicate_qc) == 0L
)


stopifnot(isTRUE(all(unlist(residual_confounding_checks))))

model_qc$residual_confounding$result_summary <- residual_confounding_qc
model_qc$residual_confounding$duplicate_keys <- residual_confounding_duplicate_qc
model_qc$residual_confounding$cautions <- c(
  "Spline-age comparisons assess age functional form only and do not establish absence of other residual confounding.",
  "E-values do not address selection bias, measurement error, model misspecification, or informative censoring.",
  "Fine-Gray SHR E-values use the EValue package's HR conversion as an exploratory approximation."
)

model_qc$checks <- c(model_qc$checks, residual_confounding_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

# The filename is retained as `evalues_<site>` for continuity.
evalue_csv_path <- file.path(table_dir, paste0("evalues_", site_name, ".csv"))
evalue_table_path <- file.path(table_dir, paste0("table_evalues_", site_name, ".html"))

write_csv(residual_confounding, evalue_csv_path)

residual_confounding_sample_note <- tibble(
  analysis = c("Mortality (60-day)", "Survival (60-day)", "28-day VFDs (liberation)"),
  n_obs = c(
    if (!is.null(mortality_models)) nrow(mortality_analysis_data) else NA_integer_,
    if (!is.null(cox_model)) nrow(cox_analysis_data) else NA_integer_,
    nrow(vfd_analysis_data)
  )
) %>%
  filter(!is.na(n_obs))

residual_confounding_sample_text <- residual_confounding_sample_note %>%
  transmute(text = paste0(analysis, ": N = ", format(n_obs, big.mark = ","))) %>%
  pull(text) %>%
  paste(collapse = "; ")

residual_confounding_sample_checks <- list(
  residual_confounding_sample_summary_nonempty = nrow(residual_confounding_sample_note) > 0L,
  residual_confounding_sample_sizes_positive = all(residual_confounding_sample_note$n_obs > 0L),
  residual_confounding_sample_text_nonempty = nzchar(residual_confounding_sample_text)
)

stopifnot(isTRUE(all(unlist(residual_confounding_sample_checks))))

model_qc$residual_confounding$sample_sizes <- residual_confounding_sample_note
model_qc$checks <- c(model_qc$checks, residual_confounding_sample_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

# Render one row per model specification and focal exposure, grouped by analysis.
evalue_gt <- residual_confounding %>%
  transmute(
    analysis,
    Model = model_spec,
    Exposure = term_label,
    Type = estimate_type,
    `Exposure SD` = sprintf("%.2f", exposure_sd),
    `Linear-age estimate` = sprintf("%.2f (%.2f, %.2f)", est_linage, lo_linage, hi_linage),
    `Spline-age estimate` = sprintf("%.2f (%.2f, %.2f)", est_splineage, lo_splineage, hi_splineage),
    `E-value (estimate)` = sprintf("%.2f", evalue_point),
    `E-value (95% CI)` = sprintf("%.2f", evalue_ci)
  ) %>%
  gt::gt(groupname_col = "analysis") %>%
  gt::tab_header(
    title = "Sensitivity analyses for nonlinear age adjustment and unmeasured confounding",
    subtitle = paste0(
      site_name, " — ratio estimates represent a one-SD increase in each focal exposure. ",
      residual_confounding_sample_text, "."
    )
  ) %>%
  gt::cols_align("left", columns = c(Model, Exposure, Type)) %>%
  gt::cols_align(
    "center",
    columns = c(
      `Exposure SD`, `Linear-age estimate`, `Spline-age estimate`,
      `E-value (estimate)`, `E-value (95% CI)`
    )
  ) %>%
  gt::sub_missing(missing_text = "—") %>%
  gt::tab_source_note(
    source_note = paste0(
      "Linear-age and spline-age estimates differ only in the functional form used for age adjustment. ",
      "This comparison evaluates sensitivity to nonlinear age but does not rule out other residual or unmeasured confounding."
    )
  ) %>%
  gt::tab_source_note(
    source_note = paste0(
      "The E-value is the minimum association strength, on the risk-ratio scale, that an unmeasured confounder would need ",
      "with both the exposure and outcome to explain the estimate, conditional on measured covariates. An E-value of 1 for ",
      "the confidence limit indicates that the confidence interval includes the null."
    )
  ) %>%
  gt::tab_source_note(
    source_note = paste0(
      "OR = odds ratio; HR = cause-specific death hazard ratio; SHR = Fine–Gray subdistribution hazard ratio for successful ",
      "liberation. SHR E-values use the EValue package's hazard-ratio conversion as an exploratory approximation."
    )
  )

gt::gtsave(evalue_gt, evalue_table_path)

evalue_csv_exists <- file.exists(evalue_csv_path)
evalue_table_exists <- file.exists(evalue_table_path)

evalue_csv_size <- if (evalue_csv_exists) {
  file.info(evalue_csv_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

evalue_table_size <- if (evalue_table_exists) {
  file.info(evalue_table_path, extra_cols = FALSE)$size[[1]]
} else {
  NA_real_
}

evalue_output_checks <- list(
  evalue_csv_exists = evalue_csv_exists,
  evalue_csv_nonempty = evalue_csv_exists && is.finite(evalue_csv_size) && evalue_csv_size > 0,
  evalue_table_exists = evalue_table_exists,
  evalue_table_nonempty = evalue_table_exists && is.finite(evalue_table_size) && evalue_table_size > 0
)

model_qc$residual_confounding$output <- list(
  csv_path = evalue_csv_path,
  csv_size_bytes = evalue_csv_size,
  html_path = evalue_table_path,
  html_size_bytes = evalue_table_size,
  n_rows = nrow(residual_confounding),
  n_analyses = n_distinct(residual_confounding$analysis)
)

model_qc$checks <- c(model_qc$checks, evalue_output_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message(
  "Sensitivity-analysis table written: ", nrow(residual_confounding),
  " ratio-scale exposure estimates across ", n_distinct(residual_confounding$analysis),
  " analyses; linear-age, spline-age, and E-value results saved."
)

# =============================================================================
# 4h. Conditional bias diagnostic plots (algorithmDiagnostics)
# =============================================================================

# `conditional_bias_plot()` is provided by the local algorithmDiagnostics package. Its interface and defaults are retained
# unchanged. The checks below validate inputs without modifying data passed to the plotting function.

conditional_bias_required_vars <- c(
  "age_group", "race_category", "sex_category", "pfvc", "pbw", "mortality_event_60", "vtpfvc", "vtpbw",
  "pbwpfvc", "ers", "crs", "vfd_28", "vtabw", "vtnchs", "abwpfvc", "nchspfvc", "abw", "pbw_nchs"
)

missing_conditional_bias_vars <- setdiff(conditional_bias_required_vars, names(cross_sectional))

stopifnot(length(missing_conditional_bias_vars) == 0L)

# Preserve the existing package input rule: patients with missing age group are excluded from all conditional-bias plots.
diag_data <- cross_sectional %>%
  filter(!is.na(age_group))

conditional_bias_variable_qc <- purrr::map_dfr(conditional_bias_required_vars, function(v) {
  x <- diag_data[[v]]
  
  tibble(
    variable = v, class = paste(class(x), collapse = "/"), n_total = length(x),
    n_nonmissing = sum(!is.na(x)), n_missing = sum(is.na(x)), pct_missing = 100 * mean(is.na(x)),
    n_unique_nonmissing = n_distinct(x, na.rm = TRUE),
    n_infinite = if (is.numeric(x)) sum(is.infinite(x), na.rm = TRUE) else NA_integer_
  )
})

conditional_bias_group_qc <- purrr::map_dfr(
  c("race_category", "sex_category", "age_group"),
  function(grouping_variable) {
    diag_data %>%
      count(.data[[grouping_variable]], name = "n", .drop = FALSE) %>%
      rename(group_value = 1) %>%
      mutate(grouping_variable = grouping_variable, .before = 1)
  }
)

conditional_bias_mortality_qc <- diag_data %>%
  count(mortality_event_60, name = "n", .drop = FALSE)

conditional_bias_input_checks <- list(
  conditional_bias_required_variables_present = length(missing_conditional_bias_vars) == 0L,
  conditional_bias_data_nonempty = nrow(diag_data) > 0L,
  conditional_bias_age_group_complete = !any(is.na(diag_data$age_group)),
  conditional_bias_age_group_levels_valid = identical(
    levels(diag_data$age_group), c("3-5", "6-8", "9-12", "13-17")
  ),
  conditional_bias_sex_levels_valid = identical(levels(diag_data$sex_category), c("Male", "Female")),
  conditional_bias_race_levels_valid = identical(
    levels(diag_data$race_category), c("WHITE", "BLACK", "OTHER")
  ),
  conditional_bias_numeric_variables_have_no_infinite_values = !any(
    conditional_bias_variable_qc$n_infinite > 0L, na.rm = TRUE
  ),
  conditional_bias_mortality_coding_valid = all(
    na.omit(diag_data$mortality_event_60) %in% c(0L, 1L)
  ),
  conditional_bias_mortality_has_variation = n_distinct(diag_data$mortality_event_60, na.rm = TRUE) == 2L
)


stopifnot(isTRUE(all(unlist(conditional_bias_input_checks))))

model_qc$conditional_bias <- list(
  package = "algorithmDiagnostics",
  input_n = nrow(diag_data),
  excluded_for_missing_age_group = nrow(cross_sectional) - nrow(diag_data),
  required_variables = conditional_bias_required_vars,
  missing_variables = missing_conditional_bias_vars,
  variable_summary = conditional_bias_variable_qc,
  grouping_summary = conditional_bias_group_qc,
  mortality_summary = conditional_bias_mortality_qc,
  note = "The package interface and package defaults were retained unchanged."
)

model_qc$checks <- c(model_qc$checks, conditional_bias_input_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Conditional-bias input QC added to: ", model_qc_path)

# Plot 1: PBW as a prediction of PFVC, stratified by demographics.
bias_pbw_pfvc <- conditional_bias_plot(
  data = diag_data, dependent_vars = "pfvc", independent_vars = "pbw",
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "PFVC (L)", min_n = 20
)

bias_pbw_pfvc_path <- file.path(figure_dir, paste0("bias_pbw_vs_pfvc_", site_name, ".pdf"))
save_script04_plot(bias_pbw_pfvc_path, bias_pbw_pfvc, width = 14, height = 12)

# Plot 2: 60-day mortality by VT/PFVC, VT/PBW, and PBW/PFVC, stratified by demographics.
bias_mortality <- conditional_bias_plot(
  data = diag_data, dependent_vars = "mortality_event_60",
  independent_vars = c(
    "vtpbw", "vtnchs", "vtabw", "vtpfvc",
    "pbwpfvc", "nchspfvc", "abwpfvc"
  ),
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "60-day mortality",
  x_labels = c(
    "Percentile of VT/McLaren PBW",
    "Percentile of VT/NCHS weight",
    "Percentile of VT/ABW",
    "Percentile of VT/PFVC",
    "Percentile of McLaren PBW/PFVC",
    "Percentile of NCHS weight/PFVC",
    "Percentile of ABW/PFVC"
  )
)

bias_mortality_path <- file.path(figure_dir, paste0("bias_mortality_", site_name, ".pdf"))
save_script04_plot(bias_mortality_path, bias_mortality, width = 14, height = 18)

# Plot 3: Elastance by exposures, stratified by demographics.
bias_ers <- conditional_bias_plot(
  data = diag_data, dependent_vars = "ers", independent_vars = c("vtpbw", "vtnchs", "vtabw", "vtpfvc"),
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "Elastance (cmH2O/L)",
  x_labels = c(
    "Percentile of VT/McLaren PBW",
    "Percentile of VT/NCHS weight",
    "Percentile of VT/ABW",
    "Percentile of VT/PFVC"
  )
)

bias_ers_path <- file.path(figure_dir, paste0("bias_elastance_", site_name, ".pdf"))
save_script04_plot(bias_ers_path, bias_ers, width = 14, height = 18)

message("Elastance bias plots saved")

# Plot 4: Compliance by exposures, stratified by demographics.
bias_crs <- conditional_bias_plot(
  data = diag_data, dependent_vars = "crs", independent_vars = c("pbwpfvc", "nchspfvc", "abwpfvc"),
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "Compliance (mL/cmH2O)", x_labels = c(
    "Percentile of McLaren PBW/PFVC",
    "Percentile of NCHS weight/PFVC",
    "Percentile of ABW/PFVC"
  )
)

bias_crs_path <- file.path(figure_dir, paste0("bias_compliance_", site_name, ".pdf"))
save_script04_plot(bias_crs_path, bias_crs, width = 14, height = 18)

message("Compliance bias plots saved")

# Plot 5: 28-day VFDs by exposures, stratified by demographics.
bias_vfd <- conditional_bias_plot(
  data = diag_data,
  dependent_vars = "vfd_28",
  independent_vars = c("pbwpfvc", "nchspfvc", "abwpfvc"),
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "28-day VFDs",
  x_labels = c(
    "Percentile of McLaren PBW/PFVC",
    "Percentile of NCHS weight/PFVC",
    "Percentile of ABW/PFVC"
  ),
  min_n = 20
)

bias_vfd_path <- file.path(figure_dir, paste0("bias_vfd28_", site_name, ".pdf"))
save_script04_plot(bias_vfd_path, bias_vfd, width = 14, height = 18)

# Plot 6a: McLaren PBW/PFVC ratio by demographics.
bias_mclaren_pfvc <- conditional_bias_plot(
  data = diag_data,
  dependent_vars = "pbwpfvc",
  independent_vars = c("pfvc", "pbw"),
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "McLaren PBW/PFVC",
  x_labels = c("Percentile of PFVC", "Percentile of McLaren PBW")
)

bias_mclaren_pfvc_path <- file.path(
  figure_dir,
  paste0("bias_mclaren_pbwpfvc_ratio_", site_name, ".pdf")
)

save_script04_plot(bias_mclaren_pfvc_path, bias_mclaren_pfvc, width = 14, height = 18)

# Plot 6b: NCHS weight/PFVC ratio by demographics.
bias_nchs_pfvc <- conditional_bias_plot(
  data = diag_data,
  dependent_vars = "nchspfvc",
  independent_vars = c("pfvc", "pbw_nchs"),
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "NCHS weight/PFVC",
  x_labels = c("Percentile of PFVC", "Percentile of NCHS weight")
)

bias_nchs_pfvc_path <- file.path(
  figure_dir,
  paste0("bias_nchs_weight_pfvc_ratio_", site_name, ".pdf")
)

save_script04_plot(bias_nchs_pfvc_path, bias_nchs_pfvc, width = 14, height = 18)

# Plot 6c: ABW/PFVC ratio by demographics.
bias_abw_pfvc <- conditional_bias_plot(
  data = diag_data,
  dependent_vars = "abwpfvc",
  independent_vars = c("pfvc", "abw"),
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "ABW/PFVC",
  x_labels = c("Percentile of PFVC", "Percentile of ABW")
)

bias_abw_pfvc_path <- file.path(
  figure_dir,
  paste0("bias_abw_pfvc_ratio_", site_name, ".pdf")
)

save_script04_plot(bias_abw_pfvc_path, bias_abw_pfvc, width = 14, height = 18)

message("All bias diagnostic plots saved")

# =============================================================================
# Presentation-only layout of returned elastance and mortality bias plots
# =============================================================================
# The collaborator's function and original plots remain unchanged.
# Existing leaf panels are assembled in one shared grid to align outcome and distribution axes.

format_returned_bias_plot <- function(original_plot, outcome_label, n_columns) {
  demographic_labels <- c("Race", "Sex", "Age group")
  panels <- list()
  preservation_checks <- list()
  
  if (!inherits(original_plot, "patchwork") || length(original_plot) != 3L) {
    stop("Expected three demographic blocks in the returned bias plot.")
  }
  
  for (group_number in seq_len(3L)) {
    demographic_block <- original_plot[[group_number]]
    
    if (!inherits(demographic_block, "patchwork") || length(demographic_block) != 2L) {
      stop("Each demographic block must contain an outcome row and a distribution row.")
    }
    
    for (row_number in seq_len(2L)) {
      panel_row <- demographic_block[[row_number]]
      
      if (!inherits(panel_row, "patchwork") || length(panel_row) != n_columns) {
        stop("The returned bias-plot row has an unexpected number of panels.")
      }
      
      for (column_number in seq_len(n_columns)) {
        original_panel <- panel_row[[column_number]]
        expected_layers <- if (row_number == 1L) 4L else 1L
        
        if (!inherits(original_panel, "ggplot") || inherits(original_panel, "patchwork") ||
            length(original_panel$layers) != expected_layers) {
          stop("The returned bias plot contains an unexpected leaf-panel structure.")
        }
        
        show_legend <- row_number == 1L && column_number == n_columns
        
        panel <- original_panel +
          ggplot2::theme(legend.position = if (show_legend) "right" else "none")
        
        if (column_number == 1L) {
          axis_label <- if (row_number == 1L) outcome_label else demographic_labels[[group_number]]
          
          panel <- panel +
            ggplot2::labs(y = axis_label) +
            ggplot2::theme(
              axis.title.y = ggplot2::element_text(size = 11),
              axis.text.y = ggplot2::element_text(size = 10),
              axis.ticks.y = ggplot2::element_line()
            )
        } else {
          panel <- panel +
            ggplot2::labs(y = NULL) +
            ggplot2::theme(
              axis.title.y = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank()
            )
        }
        
        # Reuse the original scale collection explicitly; no scale changes are intended.
        panel$scales <- original_panel$scales
        
        panels[[length(panels) + 1L]] <- panel
        
        preservation_checks[[length(preservation_checks) + 1L]] <- tibble(
          demographic_group = demographic_labels[[group_number]],
          row_type = if (row_number == 1L) "Outcome" else "Exposure distribution",
          column_number = column_number,
          legend_displayed = show_legend,
          data_preserved = identical(panel$data, original_panel$data),
          layers_preserved = identical(panel$layers, original_panel$layers),
          mapping_preserved = identical(panel$mapping, original_panel$mapping),
          coordinates_preserved = identical(panel$coordinates, original_panel$coordinates),
          scales_preserved = identical(panel$scales, original_panel$scales)
        )
      }
    }
  }
  
  preservation_checks <- dplyr::bind_rows(preservation_checks)
  
  checks <- list(
    expected_panel_count = length(panels) == 6L * n_columns,
    exactly_three_legends = sum(preservation_checks$legend_displayed) == 3L,
    legends_only_on_rightmost_outcome_panels = all(
      preservation_checks$legend_displayed ==
        (preservation_checks$row_type == "Outcome" & preservation_checks$column_number == n_columns)
    ),
    panel_data_preserved = all(preservation_checks$data_preserved),
    panel_layers_preserved = all(preservation_checks$layers_preserved),
    panel_mappings_preserved = all(preservation_checks$mapping_preserved),
    panel_coordinates_preserved = all(preservation_checks$coordinates_preserved),
    panel_scales_preserved = all(preservation_checks$scales_preserved)
  )
  
  check_values <- vapply(checks, isTRUE, logical(1))
  
  if (!all(check_values)) {
    print(preservation_checks, n = Inf, width = Inf)
    
    stop(
      "Bias presentation preservation check(s) failed: ",
      paste(names(checks)[!check_values], collapse = ", ")
    )
  }
  
  presentation_plot <- patchwork::wrap_plots(
    plotlist = panels,
    ncol = n_columns, nrow = 6, byrow = TRUE,
    widths = rep(1, n_columns),
    heights = rep(c(4, 1.3), 3),
    guides = "keep"
  )
  
  list(
    plot = presentation_plot,
    preservation_checks = preservation_checks,
    checks = checks
  )
}

bias_presentation_specs <- list(
  elastance = list(
    original_plot = bias_ers,
    outcome_label = "Elastance (cmH2O/L)",
    n_columns = 4L,
    width = 14,
    height = 18,
    filename = paste0("bias_elastance_", site_name, "_presentation.pdf")
  ),
  mortality = list(
    original_plot = bias_mortality,
    outcome_label = "60-day mortality",
    n_columns = 7L,
    width = 18,
    height = 18,
    filename = paste0("bias_mortality_", site_name, "_presentation.pdf")
  )
)

bias_presentation_results <- purrr::map(bias_presentation_specs, function(spec) {
  format_returned_bias_plot(
    original_plot = spec$original_plot,
    outcome_label = spec$outcome_label,
    n_columns = spec$n_columns
  )
})

bias_presentation_output_qc <- purrr::imap_dfr(bias_presentation_results, function(result, plot_name) {
  spec <- bias_presentation_specs[[plot_name]]
  output_path <- file.path(figure_dir, spec$filename)
  
  message("Script 04 checkpoint | rendering presentation PDF: ", plot_name)
  save_script04_plot(
    filename = output_path, plot = result$plot,
    device = "pdf", width = spec$width, height = spec$height,
    units = "in", bg = "white"
  )
  message("Script 04 checkpoint | presentation PDF render returned: ", plot_name)
  
  output_exists <- file.exists(output_path)
  output_size <- if (output_exists) file.info(output_path, extra_cols = FALSE)$size[[1]] else NA_real_
  
  tibble(
    plot_name = plot_name,
    output_path = output_path,
    width_inches = spec$width,
    height_inches = spec$height,
    output_exists = output_exists,
    output_size_bytes = output_size,
    output_nonempty = output_exists && is.finite(output_size) && output_size > 0
  )
})

bias_presentation_checks_table <- purrr::imap_dfr(bias_presentation_results, function(result, plot_name) {
  tibble(
    plot_name = plot_name,
    check = names(result$checks),
    passed = as.logical(unlist(result$checks, use.names = FALSE))
  )
})

bias_ers_presentation <- bias_presentation_results$elastance$plot
bias_mortality_presentation <- bias_presentation_results$mortality$plot

bias_presentation_checks <- list(
  bias_presentation_panel_checks_passed = isTRUE(all(bias_presentation_checks_table$passed)),
  bias_presentation_expected_outputs_created = nrow(bias_presentation_output_qc) == 2L,
  bias_presentation_outputs_exist = all(bias_presentation_output_qc$output_exists),
  bias_presentation_outputs_nonempty = all(bias_presentation_output_qc$output_nonempty)
)

stopifnot(isTRUE(all(unlist(bias_presentation_checks))))

model_qc$conditional_bias$presentation_formatting <- list(
  configuration = list(
    method = "Post-process returned leaf panels in a shared grid; collaborator's function unchanged",
    demographic_order = c("Race", "Sex", "Age group"),
    row_order = "Each outcome row immediately followed by its demographic exposure-distribution row",
    legend_placement = "Rightmost outcome panel only within each demographic block",
    original_plot_objects_modified = FALSE,
    original_pdf_files_overwritten = FALSE
  ),
  panel_preservation = purrr::imap_dfr(bias_presentation_results, function(result, plot_name) {
    result$preservation_checks %>% mutate(plot_name = plot_name, .before = 1)
  }),
  panel_checks = bias_presentation_checks_table,
  outputs = bias_presentation_output_qc
)

# Named assignment avoids duplicate checks when this formatting block is rerun interactively.
model_qc$checks[names(bias_presentation_checks)] <- bias_presentation_checks
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)
bias_presentation_readback <- readRDS(model_qc_path)

stopifnot(
  identical(bias_presentation_readback, model_qc),
  isTRUE(all(unlist(bias_presentation_checks)))
)

message("Aligned elastance and mortality presentation copies saved.")

# Validate locally rendered package outputs without modifying the algorithmDiagnostics interface or returned objects.
conditional_bias_output_qc <- tibble(
  plot_name = c(
    "PBW versus PFVC",
    "60-day mortality",
    "Elastance",
    "Compliance",
    "28-day VFDs",
    "McLaren PBW/PFVC ratio",
    "NCHS weight/PFVC ratio",
    "ABW/PFVC ratio"
  ),
  output_path = c(
    bias_pbw_pfvc_path,
    bias_mortality_path,
    bias_ers_path,
    bias_crs_path,
    bias_vfd_path,
    bias_mclaren_pfvc_path,
    bias_nchs_pfvc_path,
    bias_abw_pfvc_path
  )
) %>%
  mutate(
    output_exists = file.exists(output_path),
    output_size_bytes = purrr::map_dbl(output_path, function(path) {
      if (file.exists(path)) file.info(path, extra_cols = FALSE)$size[[1]] else NA_real_
    }),
    output_nonempty = output_exists & is.finite(output_size_bytes) & output_size_bytes > 0
  )

conditional_bias_plot_objects <- list(
  bias_pbw_pfvc = bias_pbw_pfvc,
  bias_mortality = bias_mortality,
  bias_ers = bias_ers,
  bias_crs = bias_crs,
  bias_vfd = bias_vfd,
  bias_mclaren_pfvc = bias_mclaren_pfvc,
  bias_nchs_pfvc = bias_nchs_pfvc,
  bias_abw_pfvc = bias_abw_pfvc
)


conditional_bias_plot_specifications <- bind_rows(
  # Plot 1: PFVC by McLaren PBW
  tidyr::crossing(
    plot_name = "PBW versus PFVC",
    dependent_variable = "pfvc",
    independent_variable = "pbw"
  ) %>%
    mutate(
      dependent_variable_label = "PFVC (L)",
      independent_variable_label = "Percentile of McLaren PBW",
      explicit_min_n = 20L
    ),
  
  # Plot 2: 60-day mortality by all ventilation and ratio model types
  tibble(
    plot_name = "60-day mortality",
    dependent_variable = "mortality_event_60",
    dependent_variable_label = "60-day mortality",
    independent_variable = c(
      "vtpbw",
      "vtnchs",
      "vtabw",
      "vtpfvc",
      "pbwpfvc",
      "nchspfvc",
      "abwpfvc"
    ),
    independent_variable_label = c(
      "Percentile of VT/McLaren PBW",
      "Percentile of VT/NCHS weight",
      "Percentile of VT/ABW",
      "Percentile of VT/PFVC",
      "Percentile of McLaren PBW/PFVC",
      "Percentile of NCHS weight/PFVC",
      "Percentile of ABW/PFVC"
    ),
    explicit_min_n = NA_integer_
  ),
  
  # Plot 3: elastance by all four ventilation model types
  tibble(
    plot_name = "Elastance",
    dependent_variable = "ers",
    dependent_variable_label = "Elastance (cmH2O/L)",
    independent_variable = c(
      "vtpbw",
      "vtnchs",
      "vtabw",
      "vtpfvc"
    ),
    independent_variable_label = c(
      "Percentile of VT/McLaren PBW",
      "Percentile of VT/NCHS weight",
      "Percentile of VT/ABW",
      "Percentile of VT/PFVC"
    ),
    explicit_min_n = NA_integer_
  ),
  
  # Plot 4: compliance by all three weight/PFVC ratio model types
  tibble(
    plot_name = "Compliance",
    dependent_variable = "crs",
    dependent_variable_label = "Compliance (mL/cmH2O)",
    independent_variable = c(
      "pbwpfvc",
      "nchspfvc",
      "abwpfvc"
    ),
    independent_variable_label = c(
      "Percentile of McLaren PBW/PFVC",
      "Percentile of NCHS weight/PFVC",
      "Percentile of ABW/PFVC"
    ),
    explicit_min_n = NA_integer_
  ),
  
  # Plot 5: VFDs by all three weight/PFVC ratio model types
  tibble(
    plot_name = "28-day VFDs",
    dependent_variable = "vfd_28",
    dependent_variable_label = "28-day VFDs",
    independent_variable = c(
      "pbwpfvc",
      "nchspfvc",
      "abwpfvc"
    ),
    independent_variable_label = c(
      "Percentile of McLaren PBW/PFVC",
      "Percentile of NCHS weight/PFVC",
      "Percentile of ABW/PFVC"
    ),
    explicit_min_n = 20L
  ),
  
  # Plot 6a: McLaren ratio by PFVC and McLaren PBW
  tibble(
    plot_name = "McLaren PBW/PFVC ratio",
    dependent_variable = "pbwpfvc",
    dependent_variable_label = "McLaren PBW/PFVC",
    independent_variable = c("pfvc", "pbw"),
    independent_variable_label = c(
      "Percentile of PFVC",
      "Percentile of McLaren PBW"
    ),
    explicit_min_n = NA_integer_
  ),
  
  # Plot 6b: NCHS ratio by PFVC and NCHS weight
  tibble(
    plot_name = "NCHS weight/PFVC ratio",
    dependent_variable = "nchspfvc",
    dependent_variable_label = "NCHS weight/PFVC",
    independent_variable = c("pfvc", "pbw_nchs"),
    independent_variable_label = c(
      "Percentile of PFVC",
      "Percentile of NCHS weight"
    ),
    explicit_min_n = NA_integer_
  ),
  
  # Plot 6c: ABW ratio by PFVC and ABW
  tibble(
    plot_name = "ABW/PFVC ratio",
    dependent_variable = "abwpfvc",
    dependent_variable_label = "ABW/PFVC",
    independent_variable = c("pfvc", "abw"),
    independent_variable_label = c(
      "Percentile of PFVC",
      "Percentile of ABW"
    ),
    explicit_min_n = NA_integer_
  )
) %>%
  mutate(
    grouping_variables = "race_category + sex_category + age_group",
    package_default_min_n = is.na(explicit_min_n)
  ) %>%
  select(
    plot_name,
    dependent_variable,
    dependent_variable_label,
    independent_variable,
    independent_variable_label,
    grouping_variables,
    explicit_min_n,
    package_default_min_n
  )

conditional_bias_expected_specification_count <- 24L

conditional_bias_expected_pairs <- bind_rows(
  tibble(
    plot_name = "PBW versus PFVC",
    dependent_variable = "pfvc",
    independent_variable = "pbw"
  ),
  tibble(
    plot_name = "60-day mortality",
    dependent_variable = "mortality_event_60",
    independent_variable = c(
      "vtpbw", "vtnchs", "vtabw", "vtpfvc",
      "pbwpfvc", "nchspfvc", "abwpfvc"
    )
  ),
  tibble(
    plot_name = "Elastance",
    dependent_variable = "ers",
    independent_variable = c(
      "vtpbw", "vtnchs", "vtabw", "vtpfvc"
    )
  ),
  tibble(
    plot_name = "Compliance",
    dependent_variable = "crs",
    independent_variable = c(
      "pbwpfvc", "nchspfvc", "abwpfvc"
    )
  ),
  tibble(
    plot_name = "28-day VFDs",
    dependent_variable = "vfd_28",
    independent_variable = c(
      "pbwpfvc", "nchspfvc", "abwpfvc"
    )
  ),
  tibble(
    plot_name = "McLaren PBW/PFVC ratio",
    dependent_variable = "pbwpfvc",
    independent_variable = c("pfvc", "pbw")
  ),
  tibble(
    plot_name = "NCHS weight/PFVC ratio",
    dependent_variable = "nchspfvc",
    independent_variable = c("pfvc", "pbw_nchs")
  ),
  tibble(
    plot_name = "ABW/PFVC ratio",
    dependent_variable = "abwpfvc",
    independent_variable = c("pfvc", "abw")
  )
)

conditional_bias_output_checks <- list(
  conditional_bias_expected_plot_count =
    nrow(conditional_bias_output_qc) == 8L,
  
  conditional_bias_plot_objects_created =
    length(conditional_bias_plot_objects) == 8L &&
    all(!purrr::map_lgl(conditional_bias_plot_objects, is.null)),
  
  conditional_bias_expected_specification_count =
    nrow(conditional_bias_plot_specifications) ==
    conditional_bias_expected_specification_count,
  
  conditional_bias_all_model_types_documented =
    nrow(
      anti_join(
        conditional_bias_expected_pairs,
        conditional_bias_plot_specifications %>%
          select(plot_name, dependent_variable, independent_variable),
        by = c(
          "plot_name",
          "dependent_variable",
          "independent_variable"
        )
      )
    ) == 0L,
  
  conditional_bias_no_unexpected_model_types =
    nrow(
      anti_join(
        conditional_bias_plot_specifications %>%
          select(plot_name, dependent_variable, independent_variable),
        conditional_bias_expected_pairs,
        by = c(
          "plot_name",
          "dependent_variable",
          "independent_variable"
        )
      )
    ) == 0L,
  
  conditional_bias_specification_pairs_unique =
    conditional_bias_plot_specifications %>%
    count(
      plot_name,
      dependent_variable,
      independent_variable,
      name = "n"
    ) %>%
    summarise(valid = all(n == 1L)) %>%
    pull(valid),
  
  conditional_bias_specification_plot_names_valid =
    setequal(
      unique(conditional_bias_plot_specifications$plot_name),
      conditional_bias_output_qc$plot_name
    ),
  
  conditional_bias_outputs_exist =
    all(conditional_bias_output_qc$output_exists),
  
  conditional_bias_outputs_nonempty =
    all(conditional_bias_output_qc$output_nonempty)
)


stopifnot(isTRUE(all(unlist(conditional_bias_output_checks))))

model_qc$conditional_bias$outputs <- conditional_bias_output_qc
model_qc$conditional_bias$plot_specifications <- conditional_bias_plot_specifications
model_qc$conditional_bias$interpretation <- paste0(
  "Conditional-bias plots are descriptive local diagnostics. The algorithmDiagnostics package interface and defaults ",
  "were retained unchanged. Sparse outcome events and small demographic-percentile cells may limit visual stability."
)

model_qc$checks <- c(model_qc$checks, conditional_bias_output_checks)
model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Conditional-bias plot outputs validated and added to: ", model_qc_path)

# Cross-site conditional-bias aggregation and pooled plotting are deferred.
# This script does not produce a per-percentile federated conditional-bias dataset.
# Any future export requires an explicitly approved schema and disclosure policy.
# Plotting minimum-sample settings do not establish permission to export small cells.

# =============================================================================
# 4i. Inclusion CONSORT diagram + PBW/PFVC-by-demographics figure (site QC)
# =============================================================================

# =============================================================================
# CONSORT input QC
# =============================================================================
# The cohort-flow function is retained unchanged. These checks validate the attrition log before passing it to
# render_consort().

attrition_path <- file.path(config$table_dir,paste0("attrition_log_", site_name, ".csv"))

stopifnot(file.exists(attrition_path))

attrition <- read_csv(attrition_path, show_col_types = FALSE)

attrition_required_vars <- "n_remaining"
missing_attrition_vars <- setdiff(attrition_required_vars, names(attrition))

stopifnot(length(missing_attrition_vars) == 0L)

attrition_site_valid <- if ("site" %in% names(attrition)) {
  all(!is.na(attrition$site) & attrition$site == site_name)
} else {
  TRUE
}

attrition_qc <- tibble(
  n_steps = nrow(attrition), first_step_n = first(attrition$n_remaining),
  final_step_n = last(attrition$n_remaining), cross_sectional_n = nrow(cross_sectional),
  n_missing_counts = sum(is.na(attrition$n_remaining)),
  n_negative_counts = sum(attrition$n_remaining < 0, na.rm = TRUE),
  n_increases_between_steps = sum(diff(attrition$n_remaining) > 0, na.rm = TRUE),
  final_count_matches_cross_sectional = last(attrition$n_remaining) == nrow(cross_sectional),
  site_valid = attrition_site_valid
)

attrition_checks <- list(
  attrition_file_exists = file.exists(attrition_path),
  attrition_required_variables_present = length(missing_attrition_vars) == 0L,
  attrition_has_six_steps = nrow(attrition) == 6L,
  attrition_counts_numeric = is.numeric(attrition$n_remaining),
  attrition_counts_complete = !any(is.na(attrition$n_remaining)),
  attrition_counts_nonnegative = all(attrition$n_remaining >= 0),
  attrition_counts_nonincreasing = all(diff(attrition$n_remaining) <= 0),
  attrition_final_count_matches_cross_sectional = last(attrition$n_remaining) == nrow(cross_sectional),
  attrition_site_valid = attrition_site_valid
)


stopifnot(isTRUE(all(unlist(attrition_checks))))

# =============================================================================
# CONSORT diagram
# =============================================================================

consort_fig <- render_consort(attrition, title = paste0("Cohort inclusion - ", site_name))
consort_output_path <- file.path(figure_dir, paste0("consort_diagram_", site_name, ".pdf"))

save_script04_plot(consort_output_path, consort_fig, width = 9, height = 11)

message("CONSORT diagram saved")

# =============================================================================
# PBW/PFVC distribution input QC
# =============================================================================
# This site-level figure uses local row-level data and is not included in federated exports.

distribution_required_vars <- c(
  "pbwpfvc", "nchspfvc", "abwpfvc", "sex_category", "race_category",
  "age_at_admission", "height_cm"
)

missing_distribution_vars <- setdiff(distribution_required_vars, names(cross_sectional))

stopifnot(length(missing_distribution_vars) == 0L)

dist_base <- cross_sectional %>%
  filter(
    is.finite(pbwpfvc),
    is.finite(nchspfvc),
    is.finite(abwpfvc)
  )

dist_long <- dist_base %>%
  select(
    hospitalization_id, sex_category, race_category, age_at_admission, height_cm,
    pbwpfvc, nchspfvc, abwpfvc
  ) %>%
  pivot_longer(
    cols = c(pbwpfvc, nchspfvc, abwpfvc),
    names_to = "weight_method",
    values_to = "weight_pfvc_ratio"
  ) %>%
  mutate(
    weight_method = factor(
      weight_method,
      levels = c("pbwpfvc", "nchspfvc", "abwpfvc"),
      labels = c("McLaren PBW/PFVC", "NCHS weight/PFVC", "ABW/PFVC")
    )
  )

distribution_variable_qc <- purrr::map_dfr(distribution_required_vars, function(v) {
  x <- dist_base[[v]]
  finite_values <- if (is.numeric(x)) x[is.finite(x)] else numeric(0)
  
  tibble(
    variable = v, class = paste(class(x), collapse = "/"), n_total = length(x),
    n_nonmissing = sum(!is.na(x)), n_missing = sum(is.na(x)), pct_missing = 100 * mean(is.na(x)),
    n_unique_nonmissing = n_distinct(x, na.rm = TRUE),
    n_infinite = if (is.numeric(x)) sum(is.infinite(x), na.rm = TRUE) else NA_integer_,
    minimum = if (length(finite_values) > 0L) min(finite_values) else NA_real_,
    maximum = if (length(finite_values) > 0L) max(finite_values) else NA_real_
  )
})

distribution_group_qc <- bind_rows(
  dist_base %>%
    count(sex_category, name = "n", .drop = FALSE) %>%
    transmute(grouping_variable = "sex_category", group_value = as.character(sex_category), n),
  dist_base %>%
    count(race_category, name = "n", .drop = FALSE) %>%
    transmute(grouping_variable = "race_category", group_value = as.character(race_category), n)
)

distribution_input_checks <- list(
  distribution_required_variables_present = length(missing_distribution_vars) == 0L,
  distribution_data_nonempty = nrow(dist_base) > 0L,
  distribution_pbwpfvc_complete = !any(is.na(dist_base$pbwpfvc)),
  distribution_pbwpfvc_finite = all(is.finite(dist_base$pbwpfvc)),
  distribution_pbwpfvc_positive = all(dist_base$pbwpfvc > 0),
  distribution_nchspfvc_complete = !any(is.na(dist_base$nchspfvc)),
  distribution_nchspfvc_finite = all(is.finite(dist_base$nchspfvc)),
  distribution_nchspfvc_positive = all(dist_base$nchspfvc > 0),
  distribution_abwpfvc_complete = !any(is.na(dist_base$abwpfvc)),
  distribution_abwpfvc_finite = all(is.finite(dist_base$abwpfvc)),
  distribution_abwpfvc_positive = all(dist_base$abwpfvc > 0),
  distribution_sex_complete = !any(is.na(dist_base$sex_category)),
  distribution_race_complete = !any(is.na(dist_base$race_category)),
  distribution_age_complete_and_finite = all(is.finite(dist_base$age_at_admission)),
  distribution_height_complete_and_finite = all(is.finite(dist_base$height_cm)),
  distribution_sex_levels_valid = identical(levels(dist_base$sex_category), c("Male", "Female")),
  distribution_race_levels_valid = identical(levels(dist_base$race_category), c("WHITE", "BLACK", "OTHER"))
)


stopifnot(isTRUE(all(unlist(distribution_input_checks))))

# =============================================================================
# PBW/PFVC distribution figure
# =============================================================================

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

p_sex <- ggplot(dist_long, aes(sex_category, weight_pfvc_ratio, fill = sex_category)) +
  geom_violin(alpha = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
  facet_wrap(~ weight_method, scales = "free_y") +
  scale_fill_manual(values = okabe, guide = "none") +
  labs(x = "Sex", y = "Weight/PFVC ratio (kg/L)") +
  theme_minimal(base_size = 11)

p_race <- ggplot(dist_long, aes(race_category, weight_pfvc_ratio, fill = race_category)) +
  geom_violin(alpha = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
  facet_wrap(~ weight_method, scales = "free_y") +
  scale_fill_manual(values = okabe, guide = "none") +
  labs(x = "Race", y = "Weight/PFVC ratio (kg/L)") +
  theme_minimal(base_size = 11)

p_age <- ggplot(dist_long, aes(age_at_admission, weight_pfvc_ratio)) +
  geom_point(alpha = 0.15, colour = "#D55E00") +
  geom_smooth(method = "loess", colour = "#0072B2", se = TRUE) +
  facet_wrap(~ weight_method, scales = "free_y") +
  labs(x = "Age (years)", y = "Weight/PFVC ratio (kg/L)") +
  theme_minimal(base_size = 11)

p_height <- ggplot(dist_long, aes(height_cm, weight_pfvc_ratio)) +
  geom_point(alpha = 0.15, colour = "#D55E00") +
  geom_smooth(method = "loess", colour = "#0072B2", se = TRUE) +
  facet_wrap(~ weight_method, scales = "free_y") +
  labs(x = "Height (cm)", y = "Weight/PFVC ratio (kg/L)") +
  theme_minimal(base_size = 11)

dist_fig <- (p_sex | p_race) / (p_age | p_height) +
  plot_annotation(
    title = paste0("Weight/PFVC ratios by demographics - ", site_name),
    theme = theme(plot.title = element_text(face = "bold"))
  )

distribution_output_path <- file.path(
  figure_dir,
  paste0("distribution_weight_methods_pfvc_", site_name, ".pdf")
)

save_script04_plot(distribution_output_path, dist_fig, width = 14, height = 14)

message("Weight/PFVC distribution figure saved")

# =============================================================================
# Final site-QC output validation
# =============================================================================

site_qc_output_paths <- tibble(
  output = c("CONSORT diagram", "Weight/PFVC distribution figure"),
  path = c(consort_output_path, distribution_output_path)
) %>%
  mutate(
    exists = file.exists(path),
    size_bytes = purrr::map_dbl(path, function(output_path) {
      if (file.exists(output_path)) {
        file.info(output_path, extra_cols = FALSE)$size[[1]]
      } else {
        NA_real_
      }
    }),
    nonempty = exists & is.finite(size_bytes) & size_bytes > 0
  )

site_qc_plot_checks <- list(
  consort_plot_object_created = !is.null(consort_fig),
  distribution_plot_object_created = !is.null(dist_fig),
  site_qc_output_files_exist = all(site_qc_output_paths$exists),
  site_qc_output_files_nonempty = all(site_qc_output_paths$nonempty)
)


stopifnot(isTRUE(all(unlist(site_qc_plot_checks))))

# Preserve the CSV's column values and stable metadata while excluding the
# readr problems attribute that does not retain strict identity on serialization.
attrition_for_diagnostic_rds <- attrition
attr(attrition_for_diagnostic_rds, "problems") <- NULL

stopifnot(
  identical(names(attrition_for_diagnostic_rds), names(attrition)),
  identical(class(attrition_for_diagnostic_rds), class(attrition)),
  all(vapply(
    names(attrition),
    function(column) identical(
      attrition_for_diagnostic_rds[[column]],
      attrition[[column]]
    ),
    logical(1)
  ))
)

model_qc$site_qc_figures <- list(
  attrition = list(
    input_path = attrition_path,
    summary = attrition_qc,
    data = attrition_for_diagnostic_rds
  ),
  pbwpfvc_distribution = list(
    input_n = nrow(dist_base),
    variable_summary = distribution_variable_qc,
    grouping_summary = distribution_group_qc
  ),
  outputs = site_qc_output_paths
)

model_qc$checks <- c(
  model_qc$checks,
  attrition_checks,
  distribution_input_checks,
  site_qc_plot_checks
)

model_qc$metadata$updated_at <- Sys.time()

save_model_qc(model_qc, model_qc_path)

stopifnot(file.exists(model_qc_path), isTRUE(all(unlist(model_qc$checks))))

message("Final site-QC figures validated and added to: ", model_qc_path)

# =============================================================================
# Script 04 completion
# =============================================================================
# Diagnostic report rendering is a separate workflow step.
# Reports consume the saved Table 1 and model diagnostic RDS files without refitting
# models. Script 04 does not render, copy, or delete diagnostic report files.

message("Script 04 complete.")
message("  Table artifacts: ", table_dir)
message("  Figure artifacts: ", figure_dir)
message("  Text reports: ", report_dir)
message("  Local Table 1 diagnostic RDS: ", table1_qc_path)
message("  Local model diagnostic RDS: ", model_qc_path)
message("  Diagnostic PDF rendering is deferred to the separate reporting workflow.")
