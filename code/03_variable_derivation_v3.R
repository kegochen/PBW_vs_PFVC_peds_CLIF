# =============================================================================
# Script 03: Variable Derivation
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================

rm(list = ls())

library(tidyverse)
library(arrow)
library(here)
library(data.table)
library(lubridate)
library(rspiro)

# =============================================================================
# Load validated repository configuration
# =============================================================================
# Script 00 will supply PBWPFVC_REPO_ROOT to each clean subprocess.
# For isolated execution, run from the repository root or one of its subdirectories.

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
diagnostic_dir <- config$cohort_diagnostic_dir
diagnostic_rds_path <- file.path(diagnostic_dir, "cohort_diagnostics.rds")

# Script 03 consumes existing Script 01/02 intermediates and diagnostics.
required_script03_directories <- c(output_dir, diagnostic_dir)
missing_script03_directories <- required_script03_directories[
  !dir.exists(required_script03_directories)
]

if (length(missing_script03_directories) > 0L) {
  stop(
    "Script 03 prerequisite directory(s) not found: ",
    paste(missing_script03_directories, collapse = "; "),
    ". Confirm configuration and complete Scripts 01/02 before executing Script 03."
  )
}

if (!file.exists(diagnostic_rds_path)) {
  stop("The combined Script 01/02 diagnostic RDS was not found: ", diagnostic_rds_path)
}

script03_utility_paths <- file.path(
  config$utils_dir,
  c(
    "process_sofa_scores.R",
    "standardize_pressor_dose.R",
    "attrition_log.R"
  )
)

missing_script03_utilities <- script03_utility_paths[
  !file.exists(script03_utility_paths)
]

if (length(missing_script03_utilities) > 0L) {
  stop(
    "Script 03 repository utility file(s) not found: ",
    paste(missing_script03_utilities, collapse = "; ")
  )
}

for (utility_path in script03_utility_paths) {
  source(utility_path)
}

rm(utility_path, script03_utility_paths, missing_script03_utilities)

# Ensure the configured final-output directory exists.
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(final_dir)) {
  stop("Script 03 could not create the configured final directory: ", final_dir)
}

message("Script 03 starting | site: ", site_name)
message("Intermediate input/output directory: ", output_dir)
message("Cohort diagnostic directory: ", diagnostic_dir)
message("Final output directory: ", final_dir)

# Fixed bin edges must match across sites so site-specific histograms can be combined.
PBWPFVC_BIN_EDGES <- seq(0, 50, by = 1)
AGE_BIN_EDGES <- c(3, 6, 9, 12, 15, Inf)
HEIGHT_BIN_EDGES <- c(80, 92, 102, 112, 122, 132, 142, 152, 162, 172, 182, Inf)
DIST_MIN_CELL <- 10

# =============================================================================
# Load CDC growth-reference tables
# =============================================================================
# Shared repository references are used at every site.
# Preserve the existing CSV reader and downstream interpolation methods.

stature_reference_path <- file.path(config$growth_reference_dir, "statage.csv")
weight_reference_path <- file.path(config$growth_reference_dir, "wtage.csv")

growth_reference_paths <- c(stature_reference_path, weight_reference_path)
missing_growth_reference_paths <- growth_reference_paths[
  !file.exists(growth_reference_paths)
]

if (length(missing_growth_reference_paths) > 0L) {
  stop(
    "Required growth-reference file(s) not found: ",
    paste(missing_growth_reference_paths, collapse = "; ")
  )
}

cdc_stature_for_age <- read.csv(stature_reference_path)
cdc_weight_for_age <- read.csv(weight_reference_path)

required_growth_reference_columns <- c("Sex", "Agemos", "M")

for (reference_name in c("cdc_stature_for_age", "cdc_weight_for_age")) {
  reference_table <- get(reference_name, inherits = FALSE)
  missing_reference_columns <- setdiff(
    required_growth_reference_columns,
    names(reference_table)
  )
  
  if (length(missing_reference_columns) > 0L) {
    stop(
      reference_name,
      " is missing required column(s): ",
      paste(missing_reference_columns, collapse = ", ")
    )
  }
}

rm(
  growth_reference_paths,
  missing_growth_reference_paths,
  reference_name,
  reference_table,
  missing_reference_columns
)

# =============================================================================
# McLaren predicted-body-weight method
# =============================================================================
# 1. Identify the age at which the child's height is the median stature.
# 2. Use that age to identify the corresponding median weight.

calculate_pbw_mcclaren <- function(actual_height, sex, stature_table = cdc_stature_for_age,
                                   weight_table = cdc_weight_for_age) {
  pbw <- rep(NA_real_, length(actual_height))
  
  for (s in c(1L, 2L)) {
    idx <- which(sex == s)
    
    if (length(idx) == 0L) {
      next
    }
    
    stature_sex <- stature_table %>%
      filter(Sex == s)
    
    weight_sex <- weight_table %>%
      filter(Sex == s)
    
    height_age <- approx(
      x = stature_sex$M, y = stature_sex$Agemos, xout = actual_height[idx], rule = 2
    )$y
    
    pbw[idx] <- approx(x = weight_sex$Agemos, y = weight_sex$M, xout = height_age, rule = 2)$y
  }
  
  pbw
}

# =============================================================================
# NCHS Standards predicted-weight method
# =============================================================================
# The NCHS Standards weight is the sex-specific median weight for chronological age.
# The existing CDC/NCHS weight-for-age reference table provides age in months and median weight in M.

calculate_pbw_nchs <- function(age_years, sex, weight_table = cdc_weight_for_age) {
  pbw_nchs <- rep(NA_real_, length(age_years))
  age_months <- age_years * 12
  
  for (s in c(1L, 2L)) {
    idx <- which(sex == s & !is.na(age_months))
    
    if (length(idx) == 0L) {
      next
    }
    
    weight_sex <- weight_table %>%
      filter(Sex == s) %>%
      arrange(Agemos)
    
    pbw_nchs[idx] <- approx(
      x = weight_sex$Agemos, y = weight_sex$M, xout = age_months[idx], rule = 2
    )$y
  }
  
  pbw_nchs
}

# =============================================================================
# Load cleaned intermediate data
# =============================================================================

required_input_files <- c(
  "cohort_hospitalization_ids.rds", "resp_support_waterfall_clean.parquet",
  "post_intubation_tidal_volume_rows.parquet", "cohort_demographics.parquet", "cohort_vitals_clean.parquet",
  "cohort_labs_clean.parquet", "cohort_heights_clean.parquet", "cohort_weights.parquet",
  "cohort_meds.parquet", "cohort_assessments.parquet", "attrition_log_partial.csv"
)

missing_input_files <- required_input_files[!file.exists(file.path(output_dir, required_input_files))]

if (length(missing_input_files) > 0L) {
  stop("Script 03 is missing required input file(s): ", paste(missing_input_files, collapse = ", "))
}

cohort_ids <- readRDS(file.path(output_dir, "cohort_hospitalization_ids.rds"))
resp_waterfall <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet"))
post_intubation_timepoints <- read_parquet(file.path(output_dir, "post_intubation_tidal_volume_rows.parquet"))
cohort_demographics <- read_parquet(file.path(output_dir, "cohort_demographics.parquet"))
cohort_vitals <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet"))
cohort_labs <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet"))
cohort_heights <- read_parquet(file.path(output_dir, "cohort_heights_clean.parquet"))
cohort_weights <- read_parquet(file.path(output_dir, "cohort_weights.parquet"))
cohort_meds <- read_parquet(file.path(output_dir, "cohort_meds.parquet"))
cohort_assessments <- read_parquet(file.path(output_dir, "cohort_assessments.parquet"))

required_post_intubation_columns <- c(
  ".resp_source_row_id", "intubation_event_id", "hospitalization_id", "intubation_dttm",
  "intubation_window_end_dttm", "recorded_dttm", "support_resolved",
  "intubation_identification_method", "tidal_volume_source", "tidal_volume_value"
)

missing_post_intubation_columns <- setdiff(required_post_intubation_columns, names(post_intubation_timepoints))

if (length(missing_post_intubation_columns) > 0L) {
  stop(
    "post_intubation_tidal_volume_rows.parquet is missing required column(s): ",
    paste(missing_post_intubation_columns, collapse = ", ")
  )
}

allowed_intubation_identification_methods <- c(
  "Respiratory transition",
  "Initially invasive; first documented invasive support in ICU"
)

post_intubation_input_summary <- post_intubation_timepoints %>%
  summarise(
    n_rows = n(),
    n_hospitalizations = n_distinct(hospitalization_id),
    n_intubation_events = n_distinct(intubation_event_id),
    n_set_rows = sum(tidal_volume_source == "set", na.rm = TRUE),
    n_observed_rows = sum(tidal_volume_source == "observed", na.rm = TRUE),
    n_rows_outside_cohort = sum(!hospitalization_id %in% cohort_ids),
    n_missing_identifiers = sum(
      is.na(.resp_source_row_id) | is.na(intubation_event_id) | is.na(hospitalization_id)
    ),
    n_missing_identification_method = sum(is.na(intubation_identification_method)),
    n_invalid_identification_method = sum(
      !is.na(intubation_identification_method) &
        !intubation_identification_method %in% allowed_intubation_identification_methods
    ),
    n_missing_datetimes = sum(
      is.na(intubation_dttm) | is.na(intubation_window_end_dttm) | is.na(recorded_dttm)
    ),
    n_noninvasive_rows = sum(is.na(support_resolved) | support_resolved != "invasive"),
    n_invalid_sources = sum(is.na(tidal_volume_source) | !tidal_volume_source %in% c("set", "observed")),
    n_missing_or_nonpositive_values = sum(is.na(tidal_volume_value) | tidal_volume_value <= 0),
    n_duplicate_source_rows = sum(duplicated(.resp_source_row_id)),
    n_rows_before_intubation = sum(recorded_dttm < intubation_dttm, na.rm = TRUE),
    n_rows_after_window = sum(recorded_dttm > intubation_window_end_dttm, na.rm = TRUE),
    n_window_end_mismatches = sum(intubation_window_end_dttm != intubation_dttm + hours(1), na.rm = TRUE)
  )

post_intubation_method_summary <- post_intubation_timepoints %>%
  count(intubation_identification_method, name = "n_rows") %>%
  mutate(percent = round(100 * n_rows / sum(n_rows), 1))

post_intubation_method_event_consistency <- post_intubation_timepoints %>%
  group_by(intubation_event_id) %>%
  summarise(
    n_identification_methods = n_distinct(intubation_identification_method, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  count(n_identification_methods, name = "n_intubation_events") %>%
  arrange(n_identification_methods)

post_intubation_input_checks <- list(
  required_columns_present = length(missing_post_intubation_columns) == 0L,
  input_nonempty = nrow(post_intubation_timepoints) > 0L,
  all_cohort_hospitalizations_represented =
    all(cohort_ids %in% post_intubation_timepoints$hospitalization_id),
  no_rows_outside_cohort = post_intubation_input_summary$n_rows_outside_cohort == 0L,
  identifiers_complete = post_intubation_input_summary$n_missing_identifiers == 0L,
  identification_methods_complete =
    post_intubation_input_summary$n_missing_identification_method == 0L,
  identification_methods_valid =
    post_intubation_input_summary$n_invalid_identification_method == 0L,
  one_identification_method_per_event = all(
    post_intubation_method_event_consistency$n_identification_methods == 1L
  ),
  datetimes_complete = post_intubation_input_summary$n_missing_datetimes == 0L,
  all_rows_confirmed_invasive = post_intubation_input_summary$n_noninvasive_rows == 0L,
  sources_valid = post_intubation_input_summary$n_invalid_sources == 0L,
  values_positive = post_intubation_input_summary$n_missing_or_nonpositive_values == 0L,
  source_rows_unique = post_intubation_input_summary$n_duplicate_source_rows == 0L,
  rows_within_inclusive_window =
    post_intubation_input_summary$n_rows_before_intubation == 0L &&
    post_intubation_input_summary$n_rows_after_window == 0L,
  window_end_matches_one_hour = post_intubation_input_summary$n_window_end_mismatches == 0L,
  source_counts_reconcile =
    post_intubation_input_summary$n_set_rows + post_intubation_input_summary$n_observed_rows ==
    post_intubation_input_summary$n_rows
)

post_intubation_input_checks_table <- tibble(
  check = names(post_intubation_input_checks),
  passed = as.logical(unlist(post_intubation_input_checks, use.names = FALSE))
)

message("---- Script 03 post-intubation input summary ----")
print(post_intubation_input_summary, width = Inf)
print(post_intubation_method_summary, n = Inf, width = Inf)
print(post_intubation_method_event_consistency, n = Inf, width = Inf)
print(post_intubation_input_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(post_intubation_input_checks_table$passed)))

if (!".resp_source_row_id" %in% names(resp_waterfall)) {
  stop("resp_support_waterfall_clean.parquet is missing .resp_source_row_id.")
}

waterfall_source_key_qc <- resp_waterfall %>%
  filter(!is.na(.resp_source_row_id)) %>%
  summarise(
    n_real_source_rows = n(),
    n_distinct_source_rows = n_distinct(.resp_source_row_id),
    n_duplicate_source_rows = sum(duplicated(.resp_source_row_id))
  )

if (waterfall_source_key_qc$n_duplicate_source_rows > 0L) {
  stop("The cleaned respiratory waterfall contains duplicate nonmissing respiratory source-row identifiers.")
}

required_height_provenance_columns <- c("height_source", "height_gap_days")
missing_height_provenance_columns <- setdiff(required_height_provenance_columns, names(cohort_heights))

if (length(missing_height_provenance_columns) > 0L) {
  stop(
    "cohort_heights_clean.parquet is missing required height-provenance column(s): ",
    paste(missing_height_provenance_columns, collapse = ", ")
  )
}

height_provenance_input_summary <- cohort_heights %>%
  summarise(
    n_height_rows = n(),
    n_missing_height_source = sum(is.na(height_source)),
    n_donor_source_rows = sum(height_source == "Another encounter within window", na.rm = TRUE),
    n_donor_source_rows_missing_gap = sum(
      height_source == "Another encounter within window" & is.na(height_gap_days),
      na.rm = TRUE
    ),
    n_non_donor_rows_with_missing_or_nonzero_gap = sum(
      height_source != "Another encounter within window" & (is.na(height_gap_days) | height_gap_days != 0),
      na.rm = TRUE
    )
  )

height_provenance_input_checks <- list(
  required_columns_present = length(missing_height_provenance_columns) == 0L,
  no_missing_height_source = height_provenance_input_summary$n_missing_height_source == 0L,
  donor_rows_have_gap = height_provenance_input_summary$n_donor_source_rows_missing_gap == 0L,
  non_donor_rows_have_zero_gap =
    height_provenance_input_summary$n_non_donor_rows_with_missing_or_nonzero_gap == 0L
)

height_provenance_input_checks_table <- tibble(
  check = names(height_provenance_input_checks),
  passed = as.logical(unlist(height_provenance_input_checks, use.names = FALSE))
)

message("---- Script 03 height-provenance input summary ----")
print(height_provenance_input_summary, width = Inf)
print(height_provenance_input_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(height_provenance_input_checks_table$passed)))

# Prepare intermediate measurement timestamps for downstream time matching.
cohort_vitals <- cohort_vitals %>%
  mutate(recorded_dttm = ymd_hms(recorded_dttm))

cohort_labs <- cohort_labs %>%
  mutate(lab_result_dttm = ymd_hms(lab_result_dttm))

cohort_meds <- cohort_meds %>%
  mutate(admin_dttm = ymd_hms(admin_dttm))

cohort_assessments <- cohort_assessments %>%
  mutate(recorded_dttm = ymd_hms(recorded_dttm))

message("Loaded all cleaned intermediate data")

required_weight_comparison_columns <- c(
  "hospitalization_id", "weight_kg", "weight_source", "weight_outside_expected_range"
)

missing_weight_comparison_columns <- setdiff(required_weight_comparison_columns, names(cohort_weights))

if (length(missing_weight_comparison_columns) > 0L) {
  stop(
    "cohort_weights.parquet is missing required weight-comparison column(s): ",
    paste(missing_weight_comparison_columns, collapse = ", ")
  )
}

if (anyDuplicated(cohort_weights$hospitalization_id)) {
  stop("cohort_weights.parquet contains duplicate hospitalization IDs.")
}

if (any(
  !is.na(cohort_weights$weight_kg) &
  (!is.finite(cohort_weights$weight_kg) | cohort_weights$weight_kg <= 0)
)) {
  stop("cohort_weights.parquet contains a nonfinite or nonpositive nonmissing weight_kg value.")
}

# =============================================================================
# QC: input and height-join integrity
# =============================================================================

height_join_integrity <- tibble(
  diagnostic = c("Index hospitalization IDs", "Duplicate cohort IDs",
                 "Duplicate demographic hospitalization IDs", "Duplicate cleaned-height hospitalization IDs",
                 "Cleaned-height rows outside index cohort"),
  n = c(length(cohort_ids), sum(duplicated(cohort_ids)),
        sum(duplicated(cohort_demographics$hospitalization_id)),
        sum(duplicated(cohort_heights$hospitalization_id)),
        sum(!cohort_heights$hospitalization_id %in% cohort_ids))
)

message("---- Height-join integrity ----")
print(height_join_integrity, n = Inf)

if (anyNA(cohort_ids)) {
  stop("The index cohort contains missing hospitalization IDs.")
}

if (anyDuplicated(cohort_ids)) {
  stop("The index cohort contains duplicate hospitalization IDs.")
}

if (anyDuplicated(cohort_demographics$hospitalization_id) || anyDuplicated(cohort_heights$hospitalization_id)) {
  stop("Duplicate hospitalization IDs would make the demographics-height join non-one-to-one.")
}

height_join_diagnostic <- tibble(hospitalization_id = cohort_ids) %>%
  left_join(
    cohort_heights %>% select(hospitalization_id, height_cm, height_source, height_gap_days),
    by = "hospitalization_id", relationship = "one-to-one"
  ) %>%
  mutate(
    has_usable_cleaned_height = !is.na(height_cm),
    height_join_status = if_else(has_usable_cleaned_height, "Retained; usable cleaned height",
                                 "Excluded; no usable cleaned height")
  )

height_join_summary <- height_join_diagnostic %>%
  count(height_join_status, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

height_join_funnel <- tibble(
  filter_step = c("Index cohort entering Script 03", "Usable cleaned height available"),
  n_hospitalizations = c(length(cohort_ids), sum(height_join_diagnostic$has_usable_cleaned_height))
) %>%
  mutate(n_removed = lag(n_hospitalizations) - n_hospitalizations)

message("---- Retention after cleaned-height join ----")
print(height_join_summary, n = Inf, width = Inf)
print(height_join_funnel, n = Inf)

# =============================================================================
# 3a. Calculate PBW and PFVC
# =============================================================================

pbw_pfvc_calculated <- cohort_demographics %>%
  select(hospitalization_id, patient_id, age_at_admission, sex_category, race_category, deceased) %>%
  inner_join(cohort_heights, by = "hospitalization_id", relationship = "one-to-one") %>%
  left_join(
    cohort_weights %>%
      select(
        hospitalization_id,
        abw = weight_kg,
        abw_source = weight_source,
        abw_outside_expected_range = weight_outside_expected_range
      ),
    by = "hospitalization_id",
    relationship = "one-to-one"
  ) %>%
  filter(!is.na(height_cm), height_cm >= 80, height_cm <= 200) %>%
  mutate(
    sex_numeric = case_when(
      sex_category == "Male" ~ 1L,
      sex_category == "Female" ~ 2L,
      TRUE ~ NA_integer_
    ),
    race_numeric = case_when(
      race_category == "WHITE" ~ 1L,
      race_category == "BLACK" ~ 2L,
      TRUE ~ 5L
    ),
    pbw = calculate_pbw_mcclaren(actual_height = height_cm, sex = sex_numeric),
    pbw_nchs = calculate_pbw_nchs(age_years = age_at_admission, sex = sex_numeric),
    pfvc = pred_GLI(
      age = age_at_admission, height = height_cm / 100, gender = sex_numeric, ethnicity = race_numeric,
      param = "FVC"
    ),
    abwpfvc = if_else(!is.na(abw) & !is.na(pfvc) & pfvc > 0, abw / pfvc, NA_real_),
    nchspfvc = if_else(!is.na(pbw_nchs) & !is.na(pfvc) & pfvc > 0, pbw_nchs / pfvc, NA_real_)
  )

# =============================================================================
# QC: PBW/PFVC calculation retention
# =============================================================================

pbw_pfvc_calculation_diagnostic <- pbw_pfvc_calculated %>%
  mutate(
    pbw_pfvc_status = case_when(
      is.na(sex_numeric) ~ "Excluded; sex missing or unsupported",
      is.na(pbw) & is.na(pfvc) ~ "Excluded; both PBW and PFVC missing",
      is.na(pbw) ~ "Excluded; PBW missing",
      is.na(pfvc) ~ "Excluded; PFVC missing",
      pfvc <= 0 ~ "Excluded; PFVC nonpositive",
      TRUE ~ "Retained; PBW and PFVC available"
    )
  )

pbw_pfvc_calculation_summary <- pbw_pfvc_calculation_diagnostic %>%
  count(pbw_pfvc_status, name = "n_hospitalizations") %>%
  mutate(
    percent_of_hospitalizations_with_height = round(100 * n_hospitalizations / sum(n_hospitalizations), 1)
  ) %>%
  arrange(desc(n_hospitalizations))

pbw_pfvc_calculation_funnel <- tibble(
  filter_step = c("Hospitalizations entering PBW/PFVC calculation", "PBW and positive PFVC available"),
  n_hospitalizations = c(
    nrow(pbw_pfvc_calculated),
    sum(!is.na(pbw_pfvc_calculated$pbw) & !is.na(pbw_pfvc_calculated$pfvc) & pbw_pfvc_calculated$pfvc > 0)
  )
) %>%
  mutate(n_removed = lag(n_hospitalizations) - n_hospitalizations)

message("---- PBW/PFVC calculation retention ----")
print(pbw_pfvc_calculation_summary, n = Inf, width = Inf)
print(pbw_pfvc_calculation_funnel, n = Inf)

pbw_pfvc_data <- pbw_pfvc_calculated %>%
  filter(!is.na(pbw), !is.na(pfvc), pfvc > 0) %>%
  mutate(
    sex_category = factor(sex_category, levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER"))
  )

if (any(!is.na(pbw_pfvc_data$pbw_nchs) & (!is.finite(pbw_pfvc_data$pbw_nchs) | pbw_pfvc_data$pbw_nchs <= 0))) {
  stop("At least one nonmissing NCHS predicted weight is nonfinite or nonpositive.")
}

if (any(!is.na(pbw_pfvc_data$abwpfvc) & (!is.finite(pbw_pfvc_data$abwpfvc) | pbw_pfvc_data$abwpfvc <= 0))) {
  stop("At least one nonmissing ABW/PFVC value is nonfinite or nonpositive.")
}

if (any(!is.na(pbw_pfvc_data$nchspfvc) & (!is.finite(pbw_pfvc_data$nchspfvc) | pbw_pfvc_data$nchspfvc <= 0))) {
  stop("At least one nonmissing NCHS/PFVC value is nonfinite or nonpositive.")
}

if (anyDuplicated(pbw_pfvc_data$hospitalization_id)) {
  stop("The retained PBW/PFVC dataset contains duplicate hospitalization IDs.")
}

missing_pbw_pfvc_provenance_columns <- setdiff(required_height_provenance_columns, names(pbw_pfvc_data))

if (length(missing_pbw_pfvc_provenance_columns) > 0L) {
  stop(
    "pbw_pfvc_data is missing required height-provenance column(s): ",
    paste(missing_pbw_pfvc_provenance_columns, collapse = ", ")
  )
}

pbw_pfvc_height_provenance_summary <- pbw_pfvc_data %>%
  count(height_source, name = "n_hospitalizations", .drop = FALSE) %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

pbw_pfvc_height_provenance_checks <- list(
  required_columns_present = TRUE,
  no_missing_height_source = !any(is.na(pbw_pfvc_data$height_source)),
  donor_rows_have_gap = all(
    pbw_pfvc_data$height_source != "Another encounter within window" | !is.na(pbw_pfvc_data$height_gap_days)
  ),
  non_donor_rows_have_zero_gap = all(
    pbw_pfvc_data$height_source == "Another encounter within window" |
      (!is.na(pbw_pfvc_data$height_gap_days) & pbw_pfvc_data$height_gap_days == 0)
  )
)

pbw_pfvc_height_provenance_checks_table <- tibble(
  check = names(pbw_pfvc_height_provenance_checks),
  passed = as.logical(unlist(pbw_pfvc_height_provenance_checks, use.names = FALSE))
)

message("---- PBW/PFVC height provenance ----")
print(pbw_pfvc_height_provenance_summary, n = Inf, width = Inf)
print(pbw_pfvc_height_provenance_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(pbw_pfvc_height_provenance_checks_table$passed)))

message("PBW/PFVC calculated: ", nrow(pbw_pfvc_data), " hospitalizations")

write_parquet(pbw_pfvc_data, file.path(output_dir, "analysis_broad_pfvc.parquet"))

# =============================================================================
# 3b. Calculate oxygenation metrics
# =============================================================================

OI_HYPOXEMIA_THRESHOLD <- 4
OSI_HYPOXEMIA_THRESHOLD <- 5
PF_HYPOXEMIA_THRESHOLD <- 300
SF_HYPOXEMIA_THRESHOLD <- 250

OXYGENATION_COMPONENT_MATCH_WINDOW_HOURS <- 4
OXYGENATION_TIMEPOINT_MATCH_WINDOW_HOURS <- 5

# =============================================================================
# Prepare FiO2 records
# =============================================================================

fio2_data <- resp_waterfall %>%
  filter(!is.na(fio2_set)) %>%
  transmute(
    hospitalization_id, fio2_set, fio2_dttm = as.numeric(recorded_dttm),
    fio2_source_dttm_num = as.numeric(recorded_dttm)
  )

fio2_duplicate_keys <- fio2_data %>%
  count(hospitalization_id, fio2_dttm, name = "n_rows") %>%
  filter(n_rows > 1)

if (nrow(fio2_duplicate_keys) > 0L) {
  stop("FiO2 data contain duplicate hospitalization-time keys.")
}

fio2_dt <- as.data.table(fio2_data)
setkey(fio2_dt, hospitalization_id, fio2_dttm)

fio2_availability_summary <- tibble(hospitalization_id = cohort_ids) %>%
  left_join(
    resp_waterfall %>%
      group_by(hospitalization_id) %>%
      summarise(n_waterfall_rows = n(), n_usable_fio2_rows = sum(!is.na(fio2_set)), .groups = "drop"),
    by = "hospitalization_id",
    relationship = "one-to-one"
  ) %>%
  mutate(
    across(c(n_waterfall_rows, n_usable_fio2_rows), ~ replace_na(.x, 0L)),
    fio2_status = case_when(
      n_waterfall_rows == 0 ~ "No respiratory-waterfall rows",
      n_usable_fio2_rows == 0 ~ "No usable FiO2",
      TRUE ~ "Usable FiO2 available"
    )
  ) %>%
  count(fio2_status, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

message("---- FiO2 availability ----")
print(fio2_availability_summary, n = Inf, width = Inf)

# =============================================================================
# Calculate SF ratio at each eligible SpO2 measurement
# =============================================================================

spo2_before_range_filter <- cohort_vitals %>%
  filter(vital_category == "spo2", !is.na(vital_value)) %>%
  select(hospitalization_id, recorded_dttm, vital_value)

spo2_data <- spo2_before_range_filter %>%
  filter(vital_value >= 80, vital_value <= 97) %>%
  transmute(
    hospitalization_id, spo2_row_id = row_number(), spo2_value = vital_value,
    spo2_dttm_num = as.numeric(recorded_dttm), spo2_measurement_dttm_num = as.numeric(recorded_dttm)
  )

spo2_filter_summary <- tibble(
  measure = c("Nonmissing SpO2 rows before 80-97 filter", "SpO2 rows retained from 80-97",
              "SpO2 rows removed below 80", "SpO2 rows removed above 97"),
  n = c(nrow(spo2_before_range_filter), nrow(spo2_data),
        sum(spo2_before_range_filter$vital_value < 80, na.rm = TRUE),
        sum(spo2_before_range_filter$vital_value > 97, na.rm = TRUE))
)

message("---- SpO2 filtering summary ----")
print(spo2_filter_summary, n = Inf, width = Inf)

spo2_dt <- as.data.table(spo2_data)
setkey(spo2_dt, hospitalization_id, spo2_dttm_num)

sf_joined_raw <- fio2_dt[
  spo2_dt,
  roll = OXYGENATION_COMPONENT_MATCH_WINDOW_HOURS * 3600,
  on = .(hospitalization_id, fio2_dttm = spo2_dttm_num)
] %>%
  as_tibble() %>%
  mutate(
    matched_fio2 = !is.na(fio2_set),
    fio2_lag_hours = if_else(
      matched_fio2 & !is.na(fio2_source_dttm_num) & !is.na(spo2_measurement_dttm_num),
      (spo2_measurement_dttm_num - fio2_source_dttm_num) / 3600,
      NA_real_
    )
  )

if (nrow(sf_joined_raw) != nrow(spo2_data) || anyDuplicated(sf_joined_raw$spo2_row_id)) {
  stop("The SpO2-to-FiO2 rolling join changed or duplicated SpO2 measurement rows.")
}

sf_match_summary <- tibble(
  measure = c("Eligible SpO2 measurements", "SpO2 measurements matched to FiO2 within 4 hours",
              "SpO2 measurements without a FiO2 match", "Hospitalizations with at least one SF match"),
  n = c(n_distinct(sf_joined_raw$spo2_row_id),
        n_distinct(sf_joined_raw$spo2_row_id[sf_joined_raw$matched_fio2]),
        n_distinct(sf_joined_raw$spo2_row_id) -
          n_distinct(sf_joined_raw$spo2_row_id[sf_joined_raw$matched_fio2]),
        n_distinct(sf_joined_raw$hospitalization_id[sf_joined_raw$matched_fio2]))
)

message("---- SpO2-to-FiO2 matching ----")
print(sf_match_summary, n = Inf, width = Inf)

sf_joined <- sf_joined_raw %>%
  filter(matched_fio2, !is.na(spo2_value)) %>%
  mutate(sf_ratio = spo2_value / fio2_set, sf_dttm = spo2_measurement_dttm_num)

sf_duplicate_key_detail <- sf_joined %>%
  group_by(hospitalization_id, sf_dttm) %>%
  summarise(
    n_rows = n(), n_distinct_spo2 = n_distinct(spo2_value), n_distinct_fio2 = n_distinct(fio2_set),
    .groups = "drop"
  ) %>%
  filter(n_rows > 1)

if (nrow(sf_duplicate_key_detail) > 0L && any(sf_duplicate_key_detail$n_distinct_fio2 > 1)) {
  stop("At least one duplicate SF timestamp contains conflicting FiO2 values.")
}

sf_joined_collapsed <- sf_joined %>%
  group_by(hospitalization_id, sf_dttm) %>%
  summarise(
    spo2_value = median(spo2_value, na.rm = TRUE), fio2_set = first(fio2_set),
    sf_ratio = median(sf_ratio, na.rm = TRUE), .groups = "drop"
  )

if (anyDuplicated(sf_joined_collapsed[c("hospitalization_id", "sf_dttm")])) {
  stop("Duplicate SF hospitalization-time keys remain after aggregation.")
}

# =============================================================================
# Calculate PF ratio at each arterial PaO2 measurement
# =============================================================================

pao2_data <- cohort_labs %>%
  filter(lab_category == "po2_arterial", !is.na(lab_value_numeric)) %>%
  transmute(
    hospitalization_id, pao2_row_id = row_number(), pao2_value = lab_value_numeric,
    pao2_dttm_num = as.numeric(lab_result_dttm), pao2_measurement_dttm_num = as.numeric(lab_result_dttm)
  )

pao2_dt <- as.data.table(pao2_data)
setkey(pao2_dt, hospitalization_id, pao2_dttm_num)

pf_joined_raw <- fio2_dt[
  pao2_dt,
  roll = OXYGENATION_COMPONENT_MATCH_WINDOW_HOURS * 3600,
  on = .(hospitalization_id, fio2_dttm = pao2_dttm_num)
] %>%
  as_tibble() %>%
  mutate(
    matched_fio2 = !is.na(fio2_set),
    fio2_lag_hours = if_else(
      matched_fio2 & !is.na(fio2_source_dttm_num) & !is.na(pao2_measurement_dttm_num),
      (pao2_measurement_dttm_num - fio2_source_dttm_num) / 3600,
      NA_real_
    )
  )

if (nrow(pf_joined_raw) != nrow(pao2_data) || anyDuplicated(pf_joined_raw$pao2_row_id)) {
  stop("The PaO2-to-FiO2 rolling join changed or duplicated PaO2 measurement rows.")
}

pf_match_summary <- tibble(
  measure = c("Numeric arterial PaO2 measurements", "PaO2 measurements matched to FiO2 within 4 hours",
              "PaO2 measurements without a FiO2 match", "Hospitalizations with at least one PF match"),
  n = c(n_distinct(pf_joined_raw$pao2_row_id),
        n_distinct(pf_joined_raw$pao2_row_id[pf_joined_raw$matched_fio2]),
        n_distinct(pf_joined_raw$pao2_row_id) -
          n_distinct(pf_joined_raw$pao2_row_id[pf_joined_raw$matched_fio2]),
        n_distinct(pf_joined_raw$hospitalization_id[pf_joined_raw$matched_fio2]))
)

message("---- PaO2-to-FiO2 matching ----")
print(pf_match_summary, n = Inf, width = Inf)

pf_joined <- pf_joined_raw %>%
  filter(matched_fio2, !is.na(pao2_value), is.finite(pao2_value), pao2_value > 0) %>%
  mutate(pf_ratio = pao2_value / fio2_set, pf_dttm = pao2_measurement_dttm_num)

pf_duplicate_key_detail <- pf_joined %>%
  group_by(hospitalization_id, pf_dttm) %>%
  summarise(
    n_rows = n(), n_distinct_pao2 = n_distinct(pao2_value), n_distinct_fio2 = n_distinct(fio2_set),
    .groups = "drop"
  ) %>%
  filter(n_rows > 1)

if (nrow(pf_duplicate_key_detail) > 0L && any(pf_duplicate_key_detail$n_distinct_fio2 > 1)) {
  stop("At least one duplicate PF timestamp contains conflicting FiO2 values.")
}

pf_joined_collapsed <- pf_joined %>%
  group_by(hospitalization_id, pf_dttm) %>%
  summarise(
    pao2_value = median(pao2_value, na.rm = TRUE), fio2_set = first(fio2_set),
    pf_ratio = median(pf_ratio, na.rm = TRUE), .groups = "drop"
  )

if (anyDuplicated(pf_joined_collapsed[c("hospitalization_id", "pf_dttm")])) {
  stop("Duplicate PF hospitalization-time keys remain after aggregation.")
}

# =============================================================================
# Calculate OI at arterial PaO2 measurement times
# =============================================================================

map_data <- resp_waterfall %>%
  filter(!is.na(map_set)) %>%
  transmute(
    hospitalization_id, map_set, map_dttm = as.numeric(recorded_dttm),
    map_source_dttm_num = as.numeric(recorded_dttm)
  )

map_duplicate_keys <- map_data %>%
  count(hospitalization_id, map_dttm, name = "n_rows") %>%
  filter(n_rows > 1)

if (nrow(map_duplicate_keys) > 0L) {
  stop("Mean-airway-pressure data contain duplicate hospitalization-time keys.")
}

map_dt <- as.data.table(map_data)
setkey(map_dt, hospitalization_id, map_dttm)

oi_map_joined_raw <- map_dt[
  pao2_dt,
  roll = OXYGENATION_COMPONENT_MATCH_WINDOW_HOURS * 3600,
  on = .(hospitalization_id, map_dttm = pao2_dttm_num)
] %>%
  as_tibble() %>%
  mutate(
    matched_map = !is.na(map_set),
    map_lag_hours = if_else(
      matched_map & !is.na(map_source_dttm_num) & !is.na(pao2_measurement_dttm_num),
      (pao2_measurement_dttm_num - map_source_dttm_num) / 3600,
      NA_real_
    )
  )

if (nrow(oi_map_joined_raw) != nrow(pao2_data) || anyDuplicated(oi_map_joined_raw$pao2_row_id)) {
  stop("The PaO2-to-MAP rolling join changed or duplicated PaO2 measurement rows.")
}

oi_step1 <- oi_map_joined_raw %>%
  filter(matched_map, !is.na(pao2_value)) %>%
  mutate(pao2_dttm_num = pao2_measurement_dttm_num)

oi_step1_dt <- as.data.table(oi_step1)
setkey(oi_step1_dt, hospitalization_id, pao2_dttm_num)

oi_fio2_joined_raw <- fio2_dt[
  oi_step1_dt,
  roll = OXYGENATION_COMPONENT_MATCH_WINDOW_HOURS * 3600,
  on = .(hospitalization_id, fio2_dttm = pao2_dttm_num)
] %>%
  as_tibble() %>%
  mutate(matched_fio2 = !is.na(fio2_set))

if (nrow(oi_fio2_joined_raw) != nrow(oi_step1) || anyDuplicated(oi_fio2_joined_raw$pao2_row_id)) {
  stop("The OI FiO2 rolling join changed or duplicated PaO2/MAP rows.")
}

oi_joined <- oi_fio2_joined_raw %>%
  filter(matched_fio2, !is.na(map_set), !is.na(pao2_value), is.finite(pao2_value), pao2_value > 0) %>%
  mutate(oi = fio2_set * map_set / pao2_value * 100, oi_dttm = pao2_measurement_dttm_num)

oi_duplicate_key_detail <- oi_joined %>%
  group_by(hospitalization_id, oi_dttm) %>%
  summarise(n_rows = n(), n_distinct_oi = n_distinct(oi), .groups = "drop") %>%
  filter(n_rows > 1)

if (nrow(oi_duplicate_key_detail) > 0L && any(oi_duplicate_key_detail$n_distinct_oi > 1)) {
  stop("At least one duplicate OI timestamp contains conflicting OI values.")
}

oi_joined_collapsed <- oi_joined %>%
  group_by(hospitalization_id, oi_dttm) %>%
  summarise(oi = first(oi), .groups = "drop")

if (anyDuplicated(oi_joined_collapsed[c("hospitalization_id", "oi_dttm")])) {
  stop("Duplicate OI hospitalization-time keys remain after aggregation.")
}

oi_match_summary <- tibble(
  measure = c("PaO2 measurements", "PaO2 measurements matched to MAP", "Complete OI measurements",
              "Hospitalizations with at least one OI measurement"),
  n = c(n_distinct(pao2_data$pao2_row_id), n_distinct(oi_step1$pao2_row_id),
        n_distinct(oi_joined$pao2_row_id), n_distinct(oi_joined$hospitalization_id))
)

message("---- OI matching summary ----")
print(oi_match_summary, n = Inf, width = Inf)

# =============================================================================
# Calculate OSI at eligible SpO2 measurement times
# =============================================================================

osi_map_joined_raw <- map_dt[
  spo2_dt,
  roll = OXYGENATION_COMPONENT_MATCH_WINDOW_HOURS * 3600,
  on = .(hospitalization_id, map_dttm = spo2_dttm_num)
] %>%
  as_tibble() %>%
  mutate(matched_map = !is.na(map_set))

if (nrow(osi_map_joined_raw) != nrow(spo2_data) || anyDuplicated(osi_map_joined_raw$spo2_row_id)) {
  stop("The SpO2-to-MAP rolling join changed or duplicated SpO2 measurement rows.")
}

osi_step1 <- osi_map_joined_raw %>%
  filter(matched_map, !is.na(spo2_value)) %>%
  mutate(spo2_dttm_num = spo2_measurement_dttm_num)

osi_step1_dt <- as.data.table(osi_step1)
setkey(osi_step1_dt, hospitalization_id, spo2_dttm_num)

osi_fio2_joined_raw <- fio2_dt[
  osi_step1_dt,
  roll = OXYGENATION_COMPONENT_MATCH_WINDOW_HOURS * 3600,
  on = .(hospitalization_id, fio2_dttm = spo2_dttm_num)
] %>%
  as_tibble() %>%
  mutate(matched_fio2 = !is.na(fio2_set))

if (nrow(osi_fio2_joined_raw) != nrow(osi_step1) || anyDuplicated(osi_fio2_joined_raw$spo2_row_id)) {
  stop("The OSI FiO2 rolling join changed or duplicated SpO2/MAP rows.")
}

osi_joined <- osi_fio2_joined_raw %>%
  filter(matched_fio2, !is.na(map_set), !is.na(spo2_value), is.finite(spo2_value), spo2_value > 0) %>%
  mutate(osi = fio2_set * map_set / spo2_value * 100, osi_dttm = spo2_measurement_dttm_num)

osi_duplicate_key_detail <- osi_joined %>%
  group_by(hospitalization_id, osi_dttm) %>%
  summarise(n_rows = n(), n_distinct_osi = n_distinct(osi), .groups = "drop") %>%
  filter(n_rows > 1)

if (nrow(osi_duplicate_key_detail) > 0L && any(osi_duplicate_key_detail$n_distinct_osi > 1)) {
  stop("At least one duplicate OSI timestamp contains conflicting OSI values.")
}

osi_joined_collapsed <- osi_joined %>%
  group_by(hospitalization_id, osi_dttm) %>%
  summarise(osi = first(osi), .groups = "drop")

if (anyDuplicated(osi_joined_collapsed[c("hospitalization_id", "osi_dttm")])) {
  stop("Duplicate OSI hospitalization-time keys remain after aggregation.")
}

osi_match_summary <- tibble(
  measure = c("Eligible SpO2 measurements", "SpO2 measurements matched to MAP", "Complete OSI measurements",
              "Hospitalizations with at least one OSI measurement"),
  n = c(n_distinct(spo2_data$spo2_row_id), n_distinct(osi_step1$spo2_row_id),
        n_distinct(osi_joined$spo2_row_id), n_distinct(osi_joined$hospitalization_id))
)

message("---- OSI matching summary ----")
print(osi_match_summary, n = Inf, width = Inf)

# =============================================================================
# QC: oxygenation-metric validity and distribution
# =============================================================================

summarise_oxygenation_metric <- function(data, value_column, metric_name, threshold, direction) {
  values <- data[[value_column]]
  finite_values <- values[!is.na(values) & is.finite(values)]
  positive_values <- finite_values[finite_values > 0]
  
  threshold_count <- if (direction == "greater_equal") {
    sum(finite_values >= threshold)
  } else if (direction == "less_equal") {
    sum(finite_values <= threshold)
  } else {
    stop("Unsupported threshold direction: ", direction)
  }
  
  get_quantile <- function(probability) {
    if (length(positive_values) == 0L) NA_real_ else as.numeric(quantile(positive_values, probability, names = FALSE))
  }
  
  tibble(
    metric = metric_name, n_rows = nrow(data), n_hospitalizations = n_distinct(data$hospitalization_id),
    n_missing = sum(is.na(values)), n_nonfinite = sum(!is.na(values) & !is.finite(values)),
    n_nonpositive = sum(!is.na(values) & is.finite(values) & values <= 0),
    n_positive_finite = length(positive_values),
    minimum = if (length(positive_values) > 0L) min(positive_values) else NA_real_,
    p01 = get_quantile(0.01), p05 = get_quantile(0.05), q1 = get_quantile(0.25),
    median = get_quantile(0.50), q3 = get_quantile(0.75), p95 = get_quantile(0.95),
    p99 = get_quantile(0.99), maximum = if (length(positive_values) > 0L) max(positive_values) else NA_real_,
    threshold_definition = if_else(direction == "greater_equal", paste0(">= ", threshold), paste0("<= ", threshold)),
    n_meeting_hypoxemia_threshold = threshold_count
  )
}

oxygenation_metric_distribution_summary <- bind_rows(
  summarise_oxygenation_metric(sf_joined_collapsed, "sf_ratio", "SF ratio", SF_HYPOXEMIA_THRESHOLD, "less_equal"),
  summarise_oxygenation_metric(pf_joined_collapsed, "pf_ratio", "PF ratio", PF_HYPOXEMIA_THRESHOLD, "less_equal"),
  summarise_oxygenation_metric(oi_joined_collapsed, "oi", "OI", OI_HYPOXEMIA_THRESHOLD, "greater_equal"),
  summarise_oxygenation_metric(osi_joined_collapsed, "osi", "OSI", OSI_HYPOXEMIA_THRESHOLD, "greater_equal")
)

oxygenation_metric_validity_summary <- oxygenation_metric_distribution_summary %>%
  select(metric, n_rows, n_hospitalizations, n_missing, n_nonfinite, n_nonpositive, n_positive_finite,
         n_meeting_hypoxemia_threshold)

oxygenation_metric_quantile_summary <- oxygenation_metric_distribution_summary %>%
  select(metric, minimum, p01, p05, q1, median, q3, p95, p99, maximum, threshold_definition)

message("---- Oxygenation-metric validity ----")
print(oxygenation_metric_validity_summary, n = Inf, width = Inf)

message("---- Oxygenation-metric distributions ----")
print(oxygenation_metric_quantile_summary, n = Inf, width = Inf)

rm(summarise_oxygenation_metric)

# =============================================================================
# 3c. Calculate daily pSOFA scores
# =============================================================================

admission_times <- cohort_demographics %>%
  select(hospitalization_id, admission_dttm)

age_lookup <- cohort_demographics %>%
  select(hospitalization_id, age_at_admission)

get_psofa_cv_map_threshold <- function(age_years) {
  case_when(
    age_years < 1 / 12 ~ 46L, age_years < 1 ~ 55L, age_years < 2 ~ 60L, age_years < 5 ~ 62L,
    age_years < 12 ~ 65L, age_years < 18 ~ 67L, TRUE ~ 70L
  )
}

get_psofa_renal_score <- function(creatinine, age_years) {
  case_when(
    age_years < 1 / 12 & creatinine >= 2.6 ~ 4L, age_years < 1 / 12 & creatinine >= 1.7 ~ 3L,
    age_years < 1 / 12 & creatinine >= 1.1 ~ 2L, age_years < 1 / 12 & creatinine >= 0.8 ~ 1L,
    age_years < 1 / 12 ~ 0L,
    age_years < 1 & creatinine >= 2.3 ~ 4L, age_years < 1 & creatinine >= 1.5 ~ 3L,
    age_years < 1 & creatinine >= 1.0 ~ 2L, age_years < 1 & creatinine >= 0.5 ~ 1L, age_years < 1 ~ 0L,
    age_years < 2 & creatinine >= 2.1 ~ 4L, age_years < 2 & creatinine >= 1.4 ~ 3L,
    age_years < 2 & creatinine >= 0.9 ~ 2L, age_years < 2 & creatinine >= 0.4 ~ 1L, age_years < 2 ~ 0L,
    age_years < 5 & creatinine >= 2.0 ~ 4L, age_years < 5 & creatinine >= 1.3 ~ 3L,
    age_years < 5 & creatinine >= 0.8 ~ 2L, age_years < 5 & creatinine >= 0.4 ~ 1L, age_years < 5 ~ 0L,
    age_years < 12 & creatinine >= 2.6 ~ 4L, age_years < 12 & creatinine >= 1.7 ~ 3L,
    age_years < 12 & creatinine >= 1.1 ~ 2L, age_years < 12 & creatinine >= 0.5 ~ 1L,
    age_years < 12 ~ 0L,
    age_years < 18 & creatinine >= 3.5 ~ 4L, age_years < 18 & creatinine >= 2.3 ~ 3L,
    age_years < 18 & creatinine >= 1.4 ~ 2L, age_years < 18 & creatinine >= 0.7 ~ 1L,
    age_years < 18 ~ 0L,
    creatinine >= 5.0 ~ 4L, creatinine >= 3.5 ~ 3L, creatinine >= 2.0 ~ 2L, creatinine >= 1.2 ~ 1L,
    TRUE ~ 0L
  )
}

SOFA_MAX_DAY <- 28L
PSOFA_PREADMISSION_GRACE_HOURS <- 5

add_day <- function(data, datetime_column = "recorded_dttm") {
  data %>%
    inner_join(admission_times, by = "hospitalization_id", relationship = "many-to-one") %>%
    mutate(
      sofa_day = as.integer(floor(as.numeric(difftime(.data[[datetime_column]], admission_dttm, units = "days"))))
    ) %>%
    filter(sofa_day >= 0L, sofa_day <= SOFA_MAX_DAY) %>%
    select(-admission_dttm)
}

labs_daily <- cohort_labs %>%
  filter(
    !is.na(lab_value_numeric),
    lab_category %in% c("creatinine", "bilirubin_total", "platelet_count", "po2_arterial")
  ) %>%
  select(hospitalization_id, lab_result_dttm, lab_category, lab_value_numeric) %>%
  add_day("lab_result_dttm") %>%
  group_by(hospitalization_id, sofa_day, lab_category) %>%
  summarise(
    lab_value_numeric = if_else(
      first(lab_category) %in% c("creatinine", "bilirubin_total"),
      max(lab_value_numeric, na.rm = TRUE),
      min(lab_value_numeric, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = c(hospitalization_id, sofa_day), names_from = lab_category, values_from = lab_value_numeric
  )

vitals_daily <- cohort_vitals %>%
  filter(!is.na(vital_value), vital_category %in% c("map", "spo2")) %>%
  select(hospitalization_id, recorded_dttm, vital_category, vital_value) %>%
  add_day() %>%
  group_by(hospitalization_id, sofa_day, vital_category) %>%
  summarise(vital_value = min(vital_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = c(hospitalization_id, sofa_day), names_from = vital_category, values_from = vital_value)

gcs_daily <- cohort_assessments %>%
  filter(assessment_category == "gcs_total") %>%
  mutate(gcs_total = as.numeric(numeric_value)) %>%
  select(hospitalization_id, recorded_dttm, gcs_total) %>%
  add_day() %>%
  group_by(hospitalization_id, sofa_day) %>%
  summarise(gcs_total = min(gcs_total, na.rm = TRUE), .groups = "drop")

vaso_daily <- cohort_meds %>%
  filter(
    med_category %in% c("norepinephrine", "epinephrine", "dopamine", "dobutamine"),
    med_dose_unit == "MCG/KG/MIN", !is.na(med_dose)
  ) %>%
  select(hospitalization_id, admin_dttm, med_category, med_dose) %>%
  add_day("admin_dttm") %>%
  group_by(hospitalization_id, sofa_day, med_category) %>%
  summarise(dose_mcg_kg_min = max(med_dose, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(hospitalization_id, sofa_day), names_from = med_category, values_from = dose_mcg_kg_min
  ) %>%
  rename_with(~ paste0(.x, "_mcg_kg_min"), .cols = -c(hospitalization_id, sofa_day))

resp_daily_selection <- resp_waterfall %>%
  select(hospitalization_id, recorded_dttm, device_category, fio2_set) %>%
  add_day() %>%
  mutate(device_category = na_if(str_to_title(str_trim(as.character(device_category))), "")) %>%
  group_by(hospitalization_id, sofa_day) %>%
  summarise(
    fio2_set = if (any(!is.na(fio2_set))) max(fio2_set, na.rm = TRUE) else NA_real_,
    n_distinct_device_categories = n_distinct(device_category, na.rm = TRUE),
    first_documented_device = if (all(is.na(device_category))) NA_character_ else first(na.omit(device_category)),
    selected_device_category = {
      available_devices <- sort(unique(na.omit(device_category)))
      
      if ("Imv" %in% available_devices) {
        "Imv"
      } else if ("Nippv" %in% available_devices) {
        "Nippv"
      } else if ("Cpap" %in% available_devices) {
        "Cpap"
      } else if (length(available_devices) == 0L) {
        NA_character_
      } else {
        available_devices[[1]]
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    device_selection_changed = !is.na(first_documented_device) & !is.na(selected_device_category) &
      first_documented_device != selected_device_category
  )

resp_daily_device_selection_summary <- tibble(
  measure = c("Respiratory patient-days", "Patient-days with no device category",
              "Patient-days with multiple device categories",
              "Patient-days where severity selection changed the first device",
              "Patient-days selected as IMV", "Patient-days selected as NIPPV", "Patient-days selected as CPAP"),
  n = c(nrow(resp_daily_selection), sum(is.na(resp_daily_selection$selected_device_category)),
        sum(resp_daily_selection$n_distinct_device_categories > 1),
        sum(resp_daily_selection$device_selection_changed, na.rm = TRUE),
        sum(resp_daily_selection$selected_device_category == "Imv", na.rm = TRUE),
        sum(resp_daily_selection$selected_device_category == "Nippv", na.rm = TRUE),
        sum(resp_daily_selection$selected_device_category == "Cpap", na.rm = TRUE))
)

message("---- Daily pSOFA respiratory-device selection ----")
print(resp_daily_device_selection_summary, n = Inf, width = Inf)

resp_daily <- resp_daily_selection %>%
  transmute(hospitalization_id, sofa_day, fio2_set, device_category = selected_device_category)

psofa_collapsed <- labs_daily %>%
  full_join(vitals_daily, by = c("hospitalization_id", "sofa_day"), relationship = "one-to-one") %>%
  full_join(gcs_daily, by = c("hospitalization_id", "sofa_day"), relationship = "one-to-one") %>%
  full_join(vaso_daily, by = c("hospitalization_id", "sofa_day"), relationship = "one-to-one") %>%
  full_join(resp_daily, by = c("hospitalization_id", "sofa_day"), relationship = "one-to-one") %>%
  left_join(age_lookup, by = "hospitalization_id", relationship = "many-to-one")

for (column_name in c("epinephrine_mcg_kg_min", "norepinephrine_mcg_kg_min", "dopamine_mcg_kg_min",
                      "dobutamine_mcg_kg_min")) {
  if (!column_name %in% names(psofa_collapsed)) {
    psofa_collapsed[[column_name]] <- NA_real_
  }
}

if (anyDuplicated(psofa_collapsed[c("hospitalization_id", "sofa_day")])) {
  stop("Daily pSOFA input contains duplicate hospitalization-day keys.")
}

psofa_daily <- psofa_collapsed %>%
  mutate(
    po2_arterial = if_else(between(po2_arterial, 0, 700), po2_arterial, NA_real_),
    fio2_set = if_else(between(fio2_set, 0.21, 1), fio2_set, NA_real_),
    spo2 = if_else(between(spo2, 50, 100), spo2, NA_real_)
  ) %>%
  .impute_pao2_from_spo2() %>%
  left_join(
    DEVICE_RANK_MAPPING %>%
      mutate(device_category = na_if(str_to_title(str_trim(as.character(device_category))), "")),
    by = "device_category", relationship = "many-to-one"
  ) %>%
  mutate(
    p_f = po2_arterial / fio2_set,
    p_f_imputed = pao2_imputed / fio2_set,
    p_f_for_score = coalesce(p_f, p_f_imputed),
    psofa_resp_oxygenation_source = case_when(
      !is.na(p_f) ~ "Measured PaO2",
      !is.na(p_f_imputed) ~ "SpO2-imputed PaO2",
      TRUE ~ "Unavailable"
    ),
    map_threshold = get_psofa_cv_map_threshold(age_at_admission),
    psofa_cv = case_when(
      epinephrine_mcg_kg_min > 0.1 | norepinephrine_mcg_kg_min > 0.1 ~ 4L,
      epinephrine_mcg_kg_min > 0 | norepinephrine_mcg_kg_min > 0 ~ 3L,
      dopamine_mcg_kg_min > 0 | dobutamine_mcg_kg_min > 0 ~ 2L,
      map < map_threshold ~ 1L,
      map >= map_threshold ~ 0L
    ),
    psofa_coag = case_when(
      platelet_count < 20 ~ 4L, platelet_count < 50 ~ 3L, platelet_count < 100 ~ 2L,
      platelet_count < 150 ~ 1L, platelet_count >= 150 ~ 0L
    ),
    psofa_liver = case_when(
      bilirubin_total >= 12 ~ 4L, bilirubin_total >= 6 ~ 3L, bilirubin_total >= 2 ~ 2L,
      bilirubin_total >= 1.2 ~ 1L, bilirubin_total < 1.2 ~ 0L
    ),
    psofa_resp = case_when(
      p_f_for_score < 100 & device_category %in% c("Imv", "Nippv", "Cpap") ~ 4L,
      p_f_for_score < 200 & device_category %in% c("Imv", "Nippv", "Cpap") ~ 3L,
      p_f_for_score < 300 ~ 2L,
      p_f_for_score < 400 ~ 1L,
      p_f_for_score >= 400 ~ 0L
    ),
    psofa_cns = case_when(
      gcs_total < 6 ~ 4L, gcs_total <= 9 ~ 3L, gcs_total <= 12 ~ 2L, gcs_total <= 14 ~ 1L,
      gcs_total == 15 ~ 0L
    ),
    psofa_renal = get_psofa_renal_score(creatinine, age_at_admission)
  ) %>%
  mutate(
    across(c(psofa_cv, psofa_coag, psofa_liver, psofa_resp, psofa_cns, psofa_renal), ~ coalesce(.x, 0L)),
    psofa_total = psofa_cv + psofa_coag + psofa_liver + psofa_resp + psofa_cns + psofa_renal
  )

# =============================================================================
# QC: daily respiratory pSOFA measured/imputed oxygenation and device linkage
# =============================================================================
# The measured-only score is a diagnostic counterfactual, not a replacement for
# the approved measured-first, imputed-fallback scoring rule.

score_daily_respiratory_psofa <- function(pf_ratio, device_category) {
  case_when(
    pf_ratio < 100 & device_category %in% c("Imv", "Nippv", "Cpap") ~ 4L,
    pf_ratio < 200 & device_category %in% c("Imv", "Nippv", "Cpap") ~ 3L,
    pf_ratio < 300 ~ 2L,
    pf_ratio < 400 ~ 1L,
    pf_ratio >= 400 ~ 0L,
    TRUE ~ NA_integer_
  )
}

psofa_device_mapping_qc <- DEVICE_RANK_MAPPING %>%
  transmute(
    device_category = na_if(str_to_title(str_trim(as.character(device_category))), ""),
    expected_device_rank = device_rank
  )

if (anyNA(psofa_device_mapping_qc$device_category) ||
    anyDuplicated(psofa_device_mapping_qc$device_category)) {
  stop("The normalized pSOFA device mapping contains missing or duplicate labels.")
}

psofa_respiratory_day_qc <- psofa_daily %>%
  transmute(
    hospitalization_id, sofa_day, device_category, device_rank,
    p_f, p_f_imputed, p_f_for_score, psofa_resp_oxygenation_source,
    respiratory_score_measured_only = coalesce(
      score_daily_respiratory_psofa(p_f, device_category), 0L
    ),
    respiratory_score_fallback_expected = coalesce(
      score_daily_respiratory_psofa(p_f_for_score, device_category), 0L
    ),
    respiratory_score_actual = psofa_resp,
    total_score_actual = psofa_total,
    total_score_measured_only = psofa_total - psofa_resp + respiratory_score_measured_only
  ) %>%
  left_join(psofa_device_mapping_qc, by = "device_category", relationship = "many-to-one") %>%
  mutate(
    device_mapping_covered = !is.na(expected_device_rank),
    respiratory_score_change = respiratory_score_actual - respiratory_score_measured_only,
    total_score_change = total_score_actual - total_score_measured_only
  )

psofa_respiratory_source_summary <- psofa_respiratory_day_qc %>%
  group_by(psofa_resp_oxygenation_source) %>%
  summarise(
    n_hospitalization_days = n(),
    n_hospitalizations = n_distinct(hospitalization_id),
    n_days_with_score_change = sum(respiratory_score_change != 0L),
    .groups = "drop"
  )

psofa_respiratory_score_comparison <- psofa_respiratory_day_qc %>%
  count(
    psofa_resp_oxygenation_source,
    respiratory_score_measured_only,
    respiratory_score_actual,
    respiratory_score_change,
    name = "n_hospitalization_days"
  )

psofa_device_linkage_summary <- psofa_respiratory_day_qc %>%
  count(
    device_category, device_rank, expected_device_rank, device_mapping_covered,
    name = "n_hospitalization_days"
  )

psofa_respiratory_derivation_checks <- list(
  hospitalization_day_keys_unique =
    anyDuplicated(psofa_respiratory_day_qc[c("hospitalization_id", "sofa_day")]) == 0L,
  diagnostic_row_count_preserved = nrow(psofa_respiratory_day_qc) == nrow(psofa_daily),
  scoring_ratio_matches_measured_first_fallback = identical(
    psofa_respiratory_day_qc$p_f_for_score,
    coalesce(psofa_respiratory_day_qc$p_f, psofa_respiratory_day_qc$p_f_imputed)
  ),
  oxygenation_source_matches_available_ratio = identical(
    psofa_respiratory_day_qc$psofa_resp_oxygenation_source,
    case_when(
      !is.na(psofa_respiratory_day_qc$p_f) ~ "Measured PaO2",
      !is.na(psofa_respiratory_day_qc$p_f_imputed) ~ "SpO2-imputed PaO2",
      TRUE ~ "Unavailable"
    )
  ),
  nonmissing_scoring_ratios_positive_and_finite = all(
    is.na(psofa_respiratory_day_qc$p_f_for_score) |
      (is.finite(psofa_respiratory_day_qc$p_f_for_score) &
         psofa_respiratory_day_qc$p_f_for_score > 0)
  ),
  actual_respiratory_scores_match_rule = identical(
    psofa_respiratory_day_qc$respiratory_score_actual,
    psofa_respiratory_day_qc$respiratory_score_fallback_expected
  ),
  measured_days_preserve_measured_only_score = all(
    is.na(psofa_respiratory_day_qc$p_f) |
      psofa_respiratory_day_qc$respiratory_score_change == 0L
  ),
  unavailable_days_preserve_zero_fill = all(
    psofa_respiratory_day_qc$psofa_resp_oxygenation_source != "Unavailable" |
      psofa_respiratory_day_qc$respiratory_score_actual == 0L
  ),
  covered_device_labels_have_expected_rank = all(
    !psofa_respiratory_day_qc$device_mapping_covered |
      (!is.na(psofa_respiratory_day_qc$device_rank) &
         psofa_respiratory_day_qc$device_rank == psofa_respiratory_day_qc$expected_device_rank)
  ),
  respiratory_scores_in_valid_range = all(
    psofa_respiratory_day_qc$respiratory_score_actual %in% 0:4
  ),
  total_score_change_equals_respiratory_change = all(
    psofa_respiratory_day_qc$total_score_change ==
      psofa_respiratory_day_qc$respiratory_score_change
  ),
  source_summary_counts_reconcile =
    sum(psofa_respiratory_source_summary$n_hospitalization_days) == nrow(psofa_daily),
  comparison_counts_reconcile =
    sum(psofa_respiratory_score_comparison$n_hospitalization_days) == nrow(psofa_daily)
)

psofa_respiratory_derivation_checks_table <- tibble(
  check = names(psofa_respiratory_derivation_checks),
  passed = vapply(psofa_respiratory_derivation_checks, isTRUE, logical(1))
)

if (!all(psofa_respiratory_derivation_checks_table$passed)) {
  stop(
    "Daily respiratory pSOFA QC failed: ",
    paste(
      psofa_respiratory_derivation_checks_table$check[
        !psofa_respiratory_derivation_checks_table$passed
      ],
      collapse = ", "
    )
  )
}

psofa_respiratory_derivation_qc <- list(
  policy = "Measured PF first; SpO2-imputed PF only when measured PF is unavailable",
  missing_score_policy = "Unavailable respiratory scores retain the existing zero-fill rule",
  timing_note = "Inputs use the existing daily extrema, not newly time-matched SpO2 and FiO2 pairs",
  comparison_note = paste0(
    "Measured-only totals hold all other organ scores fixed. This isolates the respiratory fallback change, ",
    "not the combined effect of earlier laboratory-cleaning changes."
  ),
  source_summary = psofa_respiratory_source_summary,
  score_comparison = psofa_respiratory_score_comparison,
  device_linkage = psofa_device_linkage_summary,
  hospitalization_day_detail = psofa_respiratory_day_qc,
  checks = psofa_respiratory_derivation_checks_table
)

message(
  "Respiratory pSOFA QC passed: ",
  sum(psofa_respiratory_day_qc$psofa_resp_oxygenation_source == "SpO2-imputed PaO2"),
  " imputed-source days; ",
  sum(psofa_respiratory_day_qc$respiratory_score_change != 0L),
  " days changed versus measured-only scoring."
)

rm(score_daily_respiratory_psofa)

sofa_daily <- psofa_daily %>%
  rename(
    sofa_cv = psofa_cv, sofa_coag = psofa_coag, sofa_liver = psofa_liver, sofa_resp = psofa_resp,
    sofa_cns = psofa_cns, sofa_renal = psofa_renal, sofa_total = psofa_total
  )

sofa_scores <- sofa_daily %>%
  filter(sofa_day == 0L) %>%
  select(-sofa_day)

message("Daily pSOFA computed for ", n_distinct(sofa_daily$hospitalization_id), " hospitalizations")

# =============================================================================
# =============================================================================
# 3d. Derive eligible post-intubation analysis timepoints
# =============================================================================

waterfall_context_columns <- c(
  ".resp_source_row_id", "device_category", "mode_category", "fio2_set", "peep_set",
  "plateau_pressure_obs", "resp_rate_set", "peak_inspiratory_pressure_obs",
  "peak_inspiratory_pressure_set", "minute_vent_obs"
)

missing_waterfall_context_columns <- setdiff(waterfall_context_columns, names(resp_waterfall))

if (length(missing_waterfall_context_columns) > 0L) {
  stop(
    "The cleaned respiratory waterfall is missing required analysis-context column(s): ",
    paste(missing_waterfall_context_columns, collapse = ", ")
  )
}

post_intubation_analysis_columns <- c(
  ".resp_source_row_id", "intubation_event_id", "hospitalization_id", "intubation_identification_method",
  "has_31500_in_hospitalization", "procedure_billed_date", "procedure_billed_dttm",
  "intubation_dttm", "intubation_window_end_dttm", "recorded_dttm", "hours_after_intubation",
  "support_resolved", "tidal_volume_source", "tidal_volume_value"
)

post_intubation_context_columns <- c(
  "device_category", "mode_category", "fio2_set", "peep_set", "plateau_pressure_obs", "resp_rate_set",
  "peak_inspiratory_pressure_obs", "peak_inspiratory_pressure_set", "minute_vent_obs"
)

required_post_intubation_analysis_columns <- c(
  post_intubation_analysis_columns,
  post_intubation_context_columns
)

missing_post_intubation_analysis_columns <- setdiff(
  required_post_intubation_analysis_columns,
  names(post_intubation_timepoints)
)

if (length(missing_post_intubation_analysis_columns) > 0L) {
  stop(
    "The post-intubation file is missing required analysis column(s): ",
    paste(missing_post_intubation_analysis_columns, collapse = ", ")
  )
}

waterfall_context <- resp_waterfall %>%
  filter(!is.na(.resp_source_row_id)) %>%
  select(all_of(waterfall_context_columns)) %>%
  rename_with(~ paste0(.x, "_clean"), -".resp_source_row_id") %>%
  mutate(.waterfall_context_matched = TRUE)

eligible_timepoints <- post_intubation_timepoints %>%
  select(all_of(required_post_intubation_analysis_columns)) %>%
  left_join(waterfall_context, by = ".resp_source_row_id", relationship = "one-to-one")

matched_waterfall_context <- eligible_timepoints$.waterfall_context_matched %in% TRUE

for (context_column in post_intubation_context_columns) {
  cleaned_context_column <- paste0(context_column, "_clean")
  
  eligible_timepoints[[context_column]] <- if_else(
    matched_waterfall_context,
    eligible_timepoints[[cleaned_context_column]],
    eligible_timepoints[[context_column]]
  )
}

# Verify that matched rows preserve cleaned values, including deliberately missing values.
matched_context_values_preserved <- vapply(
  post_intubation_context_columns,
  function(context_column) {
    cleaned_context_column <- paste0(context_column, "_clean")
    
    identical(
      eligible_timepoints[[context_column]][matched_waterfall_context],
      eligible_timepoints[[cleaned_context_column]][matched_waterfall_context]
    )
  },
  logical(1)
)

waterfall_context_preservation_checks <- tibble(
  check = paste0("Matched cleaned values preserved: ", post_intubation_context_columns),
  passed = unname(matched_context_values_preserved)
)

stopifnot(isTRUE(all(waterfall_context_preservation_checks$passed)))

eligible_timepoints <- eligible_timepoints %>%
  mutate(
    context_source = if_else(
      matched_waterfall_context,
      "Cleaned respiratory waterfall",
      "Post-intubation source-row fallback"
    )
  ) %>%
  select(-ends_with("_clean"), -.waterfall_context_matched)

rm(
  matched_waterfall_context, context_column, cleaned_context_column,
  matched_context_values_preserved
)

# =============================================================================
# Clean the selected analytic tidal volume using the Script 02 threshold rules
# =============================================================================
# Preserve the original selected value and source.
# Invalid selected values become missing; no post-cleaning source substitution is performed.

selected_tv_threshold_path <- file.path(
  config$thresholds_dir,
  "outlier_thresholds_respiratory_support.csv"
)

if (!file.exists(selected_tv_threshold_path)) {
  stop("The selected tidal-volume threshold file was not found: ", selected_tv_threshold_path)
}

selected_tv_thresholds <- read_csv(selected_tv_threshold_path, show_col_types = FALSE) %>%
  filter(variable_name %in% c("tidal_volume_set", "tidal_volume_obs")) %>%
  transmute(
    tidal_volume_source = if_else(variable_name == "tidal_volume_set", "set", "observed"),
    lower_limit = as.numeric(lower_limit),
    upper_limit = as.numeric(upper_limit)
  )

if (nrow(selected_tv_thresholds) != 2L ||
    anyDuplicated(selected_tv_thresholds$tidal_volume_source) ||
    !setequal(selected_tv_thresholds$tidal_volume_source, c("set", "observed"))) {
  stop("Exactly one threshold rule is required for each selected tidal-volume source.")
}

if (anyNA(selected_tv_thresholds$lower_limit) || anyNA(selected_tv_thresholds$upper_limit) ||
    any(!is.finite(selected_tv_thresholds$lower_limit)) ||
    any(!is.finite(selected_tv_thresholds$upper_limit)) ||
    any(selected_tv_thresholds$lower_limit > selected_tv_thresholds$upper_limit)) {
  stop("Selected tidal-volume thresholds are missing, nonfinite, or incorrectly ordered.")
}

n_timepoints_before_selected_tv_cleaning <- nrow(eligible_timepoints)

eligible_timepoints <- eligible_timepoints %>%
  mutate(tidal_volume_value_before_cleaning = tidal_volume_value) %>%
  left_join(selected_tv_thresholds, by = "tidal_volume_source", relationship = "many-to-one") %>%
  mutate(
    tidal_volume_cleaning_reason = case_when(
      is.na(tidal_volume_value_before_cleaning) ~ "Missing selected value",
      !is.finite(tidal_volume_value_before_cleaning) ~ "Nonfinite selected value",
      tidal_volume_value_before_cleaning <= 0 ~ "Nonpositive selected value",
      tidal_volume_value_before_cleaning < lower_limit ~ "Below lower threshold",
      tidal_volume_value_before_cleaning > upper_limit ~ "Above upper threshold",
      TRUE ~ "Retained"
    ),
    tidal_volume_value = if_else(
      tidal_volume_cleaning_reason == "Retained",
      tidal_volume_value_before_cleaning,
      NA_real_
    )
  ) %>%
  select(-lower_limit, -upper_limit)

selected_tv_cleaning_summary <- eligible_timepoints %>%
  count(tidal_volume_source, tidal_volume_cleaning_reason, name = "n_timepoints")

selected_tv_cleaning_by_hospitalization <- eligible_timepoints %>%
  group_by(hospitalization_id) %>%
  summarise(
    n_selected_timepoints = n(),
    n_usable_before = sum(is.finite(tidal_volume_value_before_cleaning) & tidal_volume_value_before_cleaning > 0),
    n_usable_after = sum(!is.na(tidal_volume_value)),
    lost_all_selected_tidal_volume = n_usable_before > 0L & n_usable_after == 0L,
    .groups = "drop"
  )

selected_tv_cleaning_retention <- selected_tv_cleaning_by_hospitalization %>%
  summarise(
    n_hospitalizations = n(),
    n_with_usable_selected_volume_before = sum(n_usable_before > 0L),
    n_with_usable_selected_volume_after = sum(n_usable_after > 0L),
    n_lost_all_selected_volume = sum(lost_all_selected_tidal_volume)
  )

selected_tv_cleaning_checks <- list(
  timepoint_row_count_preserved = nrow(eligible_timepoints) == n_timepoints_before_selected_tv_cleaning,
  source_rows_unique = anyDuplicated(eligible_timepoints$.resp_source_row_id) == 0L,
  sources_valid = all(eligible_timepoints$tidal_volume_source %in% c("set", "observed")),
  retained_values_positive_and_finite = all(
    is.na(eligible_timepoints$tidal_volume_value) |
      (is.finite(eligible_timepoints$tidal_volume_value) & eligible_timepoints$tidal_volume_value > 0)
  ),
  cleaning_status_matches_missingness = identical(
    eligible_timepoints$tidal_volume_cleaning_reason == "Retained",
    !is.na(eligible_timepoints$tidal_volume_value)
  ),
  retained_values_unchanged = all(
    is.na(eligible_timepoints$tidal_volume_value) |
      eligible_timepoints$tidal_volume_value == eligible_timepoints$tidal_volume_value_before_cleaning
  )
)

selected_tv_cleaning_checks_table <- tibble(
  check = names(selected_tv_cleaning_checks),
  passed = as.logical(unlist(selected_tv_cleaning_checks, use.names = FALSE))
)

message("---- Selected analytic tidal-volume cleaning ----")
print(selected_tv_thresholds, n = Inf, width = Inf)
print(selected_tv_cleaning_summary, n = Inf, width = Inf)
print(selected_tv_cleaning_retention, width = Inf)
print(selected_tv_cleaning_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(selected_tv_cleaning_checks_table$passed)))

eligible_timepoint_context_qc <- eligible_timepoints %>%
  summarise(
    n_post_intubation_rows = n(),
    n_hospitalizations = n_distinct(hospitalization_id),
    n_intubation_events = n_distinct(intubation_event_id),
    n_rows_using_cleaned_waterfall = sum(context_source == "Cleaned respiratory waterfall"),
    n_rows_using_source_fallback = sum(context_source == "Post-intubation source-row fallback"),
    n_missing_identification_method = sum(is.na(intubation_identification_method)),
    n_invalid_identification_method = sum(
      !is.na(intubation_identification_method) &
        !intubation_identification_method %in% allowed_intubation_identification_methods
    ),
    n_rows_with_all_context_values_missing = sum(
      is.na(device_category) & is.na(mode_category) & is.na(fio2_set) & is.na(peep_set) &
        is.na(plateau_pressure_obs) & is.na(resp_rate_set) & is.na(peak_inspiratory_pressure_obs) &
        is.na(peak_inspiratory_pressure_set) & is.na(minute_vent_obs)
    ),
    n_set_rows = sum(tidal_volume_source == "set"),
    n_observed_rows = sum(tidal_volume_source == "observed")
  )

eligible_timepoint_context_source_summary <- eligible_timepoints %>%
  count(context_source, name = "n_timepoints") %>%
  mutate(percent = round(100 * n_timepoints / sum(n_timepoints), 1))

eligible_timepoint_method_summary <- eligible_timepoints %>%
  count(intubation_identification_method, name = "n_rows") %>%
  mutate(percent = round(100 * n_rows / sum(n_rows), 1))

eligible_timepoint_context_checks <- list(
  one_row_per_post_intubation_source_row = nrow(eligible_timepoints) == nrow(post_intubation_timepoints),
  source_rows_unique = anyDuplicated(eligible_timepoints$.resp_source_row_id) == 0L,
  context_source_complete = !any(is.na(eligible_timepoints$context_source)),
  identification_methods_complete =
    eligible_timepoint_context_qc$n_missing_identification_method == 0L,
  identification_methods_valid =
    eligible_timepoint_context_qc$n_invalid_identification_method == 0L,
  identification_method_summary_preserved = identical(
    eligible_timepoint_method_summary,
    post_intubation_method_summary
  ),
  context_source_counts_reconcile =
    eligible_timepoint_context_qc$n_rows_using_cleaned_waterfall +
    eligible_timepoint_context_qc$n_rows_using_source_fallback ==
    eligible_timepoint_context_qc$n_post_intubation_rows,
  all_cohort_hospitalizations_retained = all(cohort_ids %in% eligible_timepoints$hospitalization_id),
  nonmissing_selected_tidal_volume_positive_and_finite = all(
    is.na(eligible_timepoints$tidal_volume_value) |
      (is.finite(eligible_timepoints$tidal_volume_value) & eligible_timepoints$tidal_volume_value > 0)
  ),
  selected_tidal_volume_cleaning_passed = isTRUE(all(selected_tv_cleaning_checks_table$passed)),
  tidal_volume_source_valid = all(eligible_timepoints$tidal_volume_source %in% c("set", "observed")),
  resolved_support_invasive = all(eligible_timepoints$support_resolved == "invasive")
)

eligible_timepoint_context_checks_table <- tibble(
  check = names(eligible_timepoint_context_checks),
  passed = as.logical(unlist(eligible_timepoint_context_checks, use.names = FALSE))
)

message("---- Eligible post-intubation timepoint context ----")
print(eligible_timepoint_context_qc, width = Inf)
print(eligible_timepoint_context_source_summary, n = Inf, width = Inf)
print(eligible_timepoint_method_summary, n = Inf, width = Inf)
print(eligible_timepoint_context_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(eligible_timepoint_context_checks_table$passed)))

# =============================================================================
# Derive respiratory-mechanics variables
# =============================================================================

analysis_data <- eligible_timepoints %>%
  inner_join(
    pbw_pfvc_data %>%
      select(
        hospitalization_id, age_at_admission, sex_category, race_category, sex_numeric, race_numeric,
        height_cm, height_source, height_gap_days, pbw, pbw_nchs, abw, abw_source,
        abw_outside_expected_range, pfvc, deceased
      ),
    by = "hospitalization_id",
    relationship = "many-to-one"
  ) %>%
  mutate(
    pbwpfvc = pbw / pfvc,
    abwpfvc = if_else(!is.na(abw) & abw > 0, abw / pfvc, NA_real_),
    nchspfvc = if_else(!is.na(pbw_nchs) & pbw_nchs > 0, pbw_nchs / pfvc, NA_real_),
    vtpbw = tidal_volume_value / pbw,
    vtabw = if_else(!is.na(abw) & abw > 0, tidal_volume_value / abw, NA_real_),
    vtnchs = if_else(!is.na(pbw_nchs) & pbw_nchs > 0, tidal_volume_value / pbw_nchs, NA_real_),
    vtpfvc = tidal_volume_value / pfvc * 0.1,
    dp_source_raw = case_when(
      !is.na(plateau_pressure_obs) ~ "plateau_obs",
      !is.na(peak_inspiratory_pressure_obs) ~ "pip_obs",
      !is.na(peak_inspiratory_pressure_set) ~ "pip_set",
      TRUE ~ NA_character_
    ),
    dp_raw = case_when(
      is.na(peep_set) ~ NA_real_,
      dp_source_raw == "plateau_obs" ~ plateau_pressure_obs - peep_set,
      dp_source_raw == "pip_obs" ~ peak_inspiratory_pressure_obs - peep_set,
      dp_source_raw == "pip_set" ~ peak_inspiratory_pressure_set - peep_set,
      TRUE ~ NA_real_
    ),
    dp = if_else(!is.na(dp_raw) & dp_raw > 0, dp_raw, NA_real_),
    dp_source = if_else(!is.na(dp), dp_source_raw, NA_character_),
    dp_type = case_when(
      dp_source == "plateau_obs" ~ "plateau-derived static DP",
      dp_source %in% c("pip_obs", "pip_set") ~ "PIP-derived fallback DP",
      TRUE ~ NA_character_
    ),
    dp_invalid_reason = case_when(
      is.na(peep_set) ~ "missing_peep_set",
      is.na(dp_source_raw) ~ "missing_pressure_source",
      !is.na(dp_raw) & dp_raw <= 0 ~ "nonpositive_pressure_minus_peep",
      TRUE ~ NA_character_
    ),
    crs = if_else(!is.na(dp) & tidal_volume_value > 0, tidal_volume_value / dp, NA_real_),
    ers = if_else(!is.na(crs) & crs > 0, 1000 / crs, NA_real_),
    ers_pbw = ers * pbw,
    ers_pfvc = ers * pfvc,
    mp_mode_class = case_when(
      mode_category == "assist control-volume control" ~ "vcv",
      mode_category %in% c("pressure control", "pressure-regulated volume control") ~ "pcv",
      TRUE ~ NA_character_
    ),
    mechanical_power = case_when(
      mp_mode_class == "vcv" & !is.na(resp_rate_set) & !is.na(peak_inspiratory_pressure_obs) & !is.na(dp) & dp > 0 ~
        0.098 * resp_rate_set * tidal_volume_value / 1000 * (peak_inspiratory_pressure_obs - 0.5 * dp),
      mp_mode_class == "pcv" & !is.na(resp_rate_set) & !is.na(peak_inspiratory_pressure_obs) ~
        0.098 * resp_rate_set * tidal_volume_value / 1000 * peak_inspiratory_pressure_obs,
      TRUE ~ NA_real_
    ),
    mp_pbw = if_else(!is.na(pbw) & pbw > 0, mechanical_power / pbw, NA_real_),
    mp_pfvc = if_else(!is.na(pfvc) & pfvc > 0, mechanical_power / pfvc, NA_real_)
  )

message(
  "Eligible post-intubation analysis timepoints: ", nrow(analysis_data), " rows across ",
  n_distinct(analysis_data$hospitalization_id), " hospitalizations"
)

dp_source_qc <- analysis_data %>%
  count(dp_source_raw, dp_source, dp_type, dp_invalid_reason, name = "n", .drop = FALSE) %>%
  mutate(percent = 100 * n / sum(n))

dp_derivation_qc <- analysis_data %>%
  summarise(
    n_timepoints = n(), n_dp_raw = sum(!is.na(dp_raw)), n_dp_valid = sum(!is.na(dp)),
    n_nonpositive_dp_raw = sum(dp_raw <= 0, na.rm = TRUE),
    n_plateau_derived = sum(dp_source == "plateau_obs", na.rm = TRUE),
    n_pip_derived = sum(dp_source %in% c("pip_obs", "pip_set"), na.rm = TRUE),
    n_crs_valid = sum(!is.na(crs)), n_ers_valid = sum(!is.na(ers)),
    minimum_dp_raw = if (any(!is.na(dp_raw))) min(dp_raw, na.rm = TRUE) else NA_real_,
    minimum_valid_dp = if (any(!is.na(dp))) min(dp, na.rm = TRUE) else NA_real_,
    maximum_valid_dp = if (any(!is.na(dp))) max(dp, na.rm = TRUE) else NA_real_
  )

dp_derivation_checks <- list(
  valid_dp_values_positive = all(analysis_data$dp > 0, na.rm = TRUE),
  crs_values_positive = all(analysis_data$crs > 0, na.rm = TRUE),
  ers_values_positive = all(analysis_data$ers > 0, na.rm = TRUE),
  crs_requires_valid_dp = all(is.na(analysis_data$crs) | !is.na(analysis_data$dp)),
  ers_requires_valid_crs = all(is.na(analysis_data$ers) | !is.na(analysis_data$crs)),
  dp_source_present_when_dp_valid = all(is.na(analysis_data$dp) | !is.na(analysis_data$dp_source)),
  dp_type_present_when_dp_valid = all(is.na(analysis_data$dp) | !is.na(analysis_data$dp_type)),
  valid_dp_matches_raw_dp = all(is.na(analysis_data$dp) | analysis_data$dp == analysis_data$dp_raw)
)

print(dp_derivation_qc)
print(dp_source_qc, n = Inf)

stopifnot(
  isTRUE(all(unlist(dp_derivation_checks)))
)

message(
  "Valid driving-pressure source mix: ",
  paste(capture.output(print(table(analysis_data$dp_source, useNA = "ifany"))), collapse = " | ")
)

# =============================================================================
# Norepinephrine-equivalent vasopressor dose
# =============================================================================

catecholamines <- c("norepinephrine", "epinephrine", "dopamine", "phenylephrine", "dobutamine")

ne_equiv_cat <- cohort_meds %>%
  filter(med_category %in% catecholamines, !is.na(med_dose)) %>%
  standardize_pressor_dose(weights = cohort_weights, out_col = "dose_mcg_kg_min", label = "NE-equiv catecholamines") %>%
  filter(!is.na(dose_mcg_kg_min)) %>%
  mutate(
    ne_equiv = case_when(
      med_category == "norepinephrine" ~ dose_mcg_kg_min,
      med_category == "epinephrine" ~ dose_mcg_kg_min,
      med_category == "dopamine" ~ dose_mcg_kg_min / 100,
      med_category == "phenylephrine" ~ dose_mcg_kg_min / 10,
      med_category == "dobutamine" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  select(hospitalization_id, admin_dttm, ne_equiv)

ne_equiv_vaso <- cohort_meds %>%
  filter(med_category == "vasopressin", !is.na(med_dose)) %>%
  mutate(
    med_dose = as.numeric(med_dose),
    unit_clean = str_to_lower(str_replace_all(med_dose_unit, "\\s+", "")),
    units_per_min = case_when(
      unit_clean %in% c("units/min", "u/min") ~ med_dose,
      unit_clean %in% c("units/hr", "u/hr", "units/h", "u/h") ~ med_dose / 60,
      TRUE ~ NA_real_
    ),
    ne_equiv = units_per_min * 2.5
  ) %>%
  filter(!is.na(ne_equiv)) %>%
  select(hospitalization_id, admin_dttm, ne_equiv)

ne_equiv <- bind_rows(ne_equiv_cat, ne_equiv_vaso) %>%
  filter(ne_equiv > 0) %>%
  group_by(hospitalization_id, admin_dttm) %>%
  summarise(ne_equiv_total = sum(ne_equiv, na.rm = TRUE), .groups = "drop")

ne_duplicate_keys <- ne_equiv %>%
  count(hospitalization_id, admin_dttm, name = "n_rows") %>%
  filter(n_rows > 1)

if (nrow(ne_duplicate_keys) > 0L) {
  stop("Norepinephrine-equivalent data contain duplicate hospitalization-time keys.")
}

ne_dt <- as.data.table(ne_equiv)
ne_dt[, join_dttm := as.numeric(admin_dttm)]
setkey(ne_dt, hospitalization_id, join_dttm)

analysis_ne_dt <- analysis_data %>%
  mutate(.analysis_row_id = row_number(), join_dttm = as.numeric(recorded_dttm)) %>%
  as.data.table()

n_analysis_rows_before_ne_join <- nrow(analysis_ne_dt)
setkey(analysis_ne_dt, hospitalization_id, join_dttm)

analysis_ne_joined <- ne_dt[, .(hospitalization_id, join_dttm, ne_equiv_total)][
  analysis_ne_dt, roll = 4 * 3600, on = .(hospitalization_id, join_dttm)
]

if (nrow(analysis_ne_joined) != n_analysis_rows_before_ne_join ||
    anyDuplicated(analysis_ne_joined$.analysis_row_id)) {
  stop("The norepinephrine-equivalent rolling join changed or duplicated IMV timepoints.")
}

analysis_data <- as_tibble(analysis_ne_joined) %>%
  arrange(.analysis_row_id) %>%
  mutate(ne_equiv_total = coalesce(ne_equiv_total, 0)) %>%
  select(-.analysis_row_id, -join_dttm)

message("NE equivalents joined: ", sum(analysis_data$ne_equiv_total > 0), " timepoints with vasopressors")

# =============================================================================
# Ventilatory ratio
# =============================================================================

VR_VE_WINDOW_HOURS <- 4
VR_TIMEPOINT_WINDOW_HOURS <- 1

pco2_data <- cohort_labs %>%
  filter(lab_category == "pco2_arterial", !is.na(lab_value_numeric)) %>%
  transmute(
    hospitalization_id, pco2_row_id = row_number(), paco2 = lab_value_numeric,
    pco2_dttm_num = as.numeric(lab_result_dttm), pco2_measurement_dttm_num = as.numeric(lab_result_dttm)
  )

if (nrow(pco2_data) > 0L) {
  ve_data <- resp_waterfall %>%
    filter(!is.na(minute_vent_obs)) %>%
    transmute(
      hospitalization_id, minute_vent_obs, ve_dttm_num = as.numeric(recorded_dttm),
      ve_source_dttm_num = as.numeric(recorded_dttm)
    )
  
  ve_duplicate_keys <- ve_data %>%
    count(hospitalization_id, ve_dttm_num, name = "n_rows") %>%
    filter(n_rows > 1)
  
  if (nrow(ve_duplicate_keys) > 0L) {
    stop("Minute-ventilation data contain duplicate hospitalization-time keys.")
  }
  
  ve_dt <- as.data.table(ve_data)
  setkey(ve_dt, hospitalization_id, ve_dttm_num)
  
  pco2_dt <- as.data.table(copy(pco2_data))
  setkey(pco2_dt, hospitalization_id, pco2_dttm_num)
  
  vr_at_pco2_raw <- ve_dt[
    pco2_dt, roll = VR_VE_WINDOW_HOURS * 3600,
    on = .(hospitalization_id, ve_dttm_num = pco2_dttm_num)
  ]
  
  if (nrow(vr_at_pco2_raw) != nrow(pco2_data) || anyDuplicated(vr_at_pco2_raw$pco2_row_id)) {
    stop("The minute-ventilation/PaCO2 rolling join changed or duplicated PaCO2 rows.")
  }
  
  vr_at_pco2 <- as_tibble(vr_at_pco2_raw) %>%
    left_join(
      pbw_pfvc_data %>% select(hospitalization_id, pbw),
      by = "hospitalization_id",
      relationship = "many-to-one"
    ) %>%
    mutate(
      ve_lag_hours = if_else(
        !is.na(ve_source_dttm_num) & !is.na(pco2_measurement_dttm_num),
        (pco2_measurement_dttm_num - ve_source_dttm_num) / 3600,
        NA_real_
      ),
      vent_ratio = if_else(
        !is.na(minute_vent_obs) & !is.na(paco2) & !is.na(pbw) & pbw > 0,
        minute_vent_obs * paco2 / (pbw * 0.1 * 37.5),
        NA_real_
      ),
      vr_dttm = pco2_measurement_dttm_num
    ) %>%
    filter(!is.na(vent_ratio)) %>%
    transmute(
      hospitalization_id, pco2_row_id, vr_dttm, vent_ratio, vr_minute_vent = minute_vent_obs,
      vr_paco2 = paco2, ve_lag_hours
    )
  
  vr_duplicate_key_detail <- vr_at_pco2 %>%
    group_by(hospitalization_id, vr_dttm) %>%
    summarise(
      n_rows = n(), n_distinct_minute_vent = n_distinct(vr_minute_vent),
      n_distinct_paco2 = n_distinct(vr_paco2), .groups = "drop"
    ) %>%
    filter(n_rows > 1)
  
  if (nrow(vr_duplicate_key_detail) > 0L && any(vr_duplicate_key_detail$n_distinct_minute_vent > 1)) {
    stop("At least one duplicate VR timestamp contains conflicting minute-ventilation values.")
  }
  
  vr_at_pco2_collapsed <- vr_at_pco2 %>%
    group_by(hospitalization_id, vr_dttm) %>%
    summarise(
      vent_ratio = median(vent_ratio, na.rm = TRUE), vr_minute_vent = first(vr_minute_vent),
      vr_paco2 = median(vr_paco2, na.rm = TRUE), ve_lag_hours = median(ve_lag_hours, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (anyDuplicated(vr_at_pco2_collapsed[c("hospitalization_id", "vr_dttm")])) {
    stop("Duplicate VR hospitalization-time keys remain after aggregation.")
  }
  
  vr_for_join <- as.data.table(vr_at_pco2_collapsed)
  setkey(vr_for_join, hospitalization_id, vr_dttm)
  
  analysis_vr_dt <- analysis_data %>%
    mutate(.analysis_row_id = row_number(), recorded_dttm_num = as.numeric(recorded_dttm)) %>%
    as.data.table()
  
  n_analysis_rows_before_vr_join <- nrow(analysis_vr_dt)
  setkey(analysis_vr_dt, hospitalization_id, recorded_dttm_num)
  
  analysis_vr <- vr_for_join[
    analysis_vr_dt, roll = VR_TIMEPOINT_WINDOW_HOURS * 3600,
    on = .(hospitalization_id, vr_dttm = recorded_dttm_num)
  ]
  
  if (nrow(analysis_vr) != n_analysis_rows_before_vr_join || anyDuplicated(analysis_vr$.analysis_row_id)) {
    stop("The ventilatory-ratio rolling join changed or duplicated IMV timepoints.")
  }
  
  analysis_data <- as_tibble(analysis_vr) %>%
    arrange(.analysis_row_id) %>%
    select(-.analysis_row_id, -vr_dttm)
  
  message("Ventilatory ratio available at ", sum(!is.na(analysis_data$vent_ratio)), " of ",
          nrow(analysis_data), " timepoints")
} else {
  analysis_data <- analysis_data %>%
    mutate(
      vent_ratio = NA_real_, vr_minute_vent = NA_real_, vr_paco2 = NA_real_, ve_lag_hours = NA_real_
    )
  
  message("No PaCO2 data available; ventilatory ratio set to missing.")
}

# =============================================================================
# Join oxygenation metrics to IMV timepoints
# =============================================================================

join_source_to_timepoints <- function(analysis, source_data, source_time, value_columns, window_hours, label) {
  source <- source_data %>%
    select(hospitalization_id, all_of(source_time), all_of(value_columns)) %>%
    rename(.join_dttm = all_of(source_time)) %>%
    as.data.table()
  
  source_key_columns <- c("hospitalization_id", ".join_dttm")
  
  if (anyDuplicated(source[, ..source_key_columns])) {
    stop(label, " source data contain duplicate hospitalization-time keys.")
  }
  
  setkeyv(source, source_key_columns)
  
  target <- analysis %>%
    mutate(.analysis_row_id = row_number(), .join_dttm = as.numeric(recorded_dttm)) %>%
    as.data.table()
  
  n_rows_before <- nrow(target)
  setkey(target, hospitalization_id, .join_dttm)
  
  joined <- source[target, roll = window_hours * 3600, on = .(hospitalization_id, .join_dttm)]
  
  if (nrow(joined) != n_rows_before || anyDuplicated(joined$.analysis_row_id)) {
    stop("The ", label, " rolling join changed or duplicated IMV timepoints.")
  }
  
  as_tibble(joined) %>%
    arrange(.analysis_row_id) %>%
    select(-.analysis_row_id, -.join_dttm)
}

join_metric_to_timepoints <- function(analysis, metric_data, source_time, metric_column, window_hours, label) {
  source <- metric_data %>%
    select(hospitalization_id, all_of(source_time), all_of(metric_column)) %>%
    rename(.source_dttm = all_of(source_time)) %>%
    mutate(.matched_source_dttm = .source_dttm) %>%
    as.data.table()
  
  if (anyDuplicated(source[, .(hospitalization_id, .source_dttm)])) {
    stop(label, " source data contain duplicate hospitalization-time keys.")
  }
  
  setkey(source, hospitalization_id, .source_dttm)
  
  target <- analysis %>%
    mutate(.analysis_row_id = row_number(), .target_dttm = as.numeric(recorded_dttm)) %>%
    as.data.table()
  
  n_rows_before <- nrow(target)
  setkey(target, hospitalization_id, .target_dttm)
  
  prior <- source[
    target,
    on = .(hospitalization_id, .source_dttm = .target_dttm),
    roll = window_hours * 3600
  ]
  
  next_value <- source[
    target,
    on = .(hospitalization_id, .source_dttm = .target_dttm),
    roll = -window_hours * 3600
  ]
  
  if (
    nrow(prior) != n_rows_before || nrow(next_value) != n_rows_before ||
    anyDuplicated(prior$.analysis_row_id) || anyDuplicated(next_value$.analysis_row_id)
  ) {
    stop("The ", label, " bidirectional rolling joins changed or duplicated IMV timepoints.")
  }
  
  metric_prior <- paste0(metric_column, "_prior")
  metric_next <- paste0(metric_column, "_next")
  source_prior <- paste0(metric_column, "_source_dttm_prior")
  source_next <- paste0(metric_column, "_source_dttm_next")
  matched_source <- paste0(metric_column, "_matched_dttm_num")
  matched_offset <- paste0(metric_column, "_match_offset_hours")
  matched_direction <- paste0(metric_column, "_match_direction")
  
  prior_selected <- as_tibble(prior) %>%
    transmute(
      .analysis_row_id,
      !!metric_prior := .data[[metric_column]],
      !!source_prior := .matched_source_dttm
    )
  
  next_selected <- as_tibble(next_value) %>%
    transmute(
      .analysis_row_id,
      !!metric_next := .data[[metric_column]],
      !!source_next := .matched_source_dttm
    )
  
  result <- target %>%
    as_tibble() %>%
    left_join(prior_selected, by = ".analysis_row_id", relationship = "one-to-one") %>%
    left_join(next_selected, by = ".analysis_row_id", relationship = "one-to-one") %>%
    mutate(
      .prior_offset_hours = (.data[[source_prior]] - .target_dttm) / 3600,
      .next_offset_hours = (.data[[source_next]] - .target_dttm) / 3600,
      .use_prior = case_when(
        !is.na(.data[[metric_prior]]) & is.na(.data[[metric_next]]) ~ TRUE,
        is.na(.data[[metric_prior]]) & !is.na(.data[[metric_next]]) ~ FALSE,
        !is.na(.data[[metric_prior]]) & !is.na(.data[[metric_next]]) ~
          abs(.prior_offset_hours) <= abs(.next_offset_hours),
        TRUE ~ NA
      ),
      !!metric_column := case_when(
        .use_prior %in% TRUE ~ as.numeric(.data[[metric_prior]]),
        .use_prior %in% FALSE ~ as.numeric(.data[[metric_next]]),
        TRUE ~ NA_real_
      ),
      !!matched_source := case_when(
        .use_prior %in% TRUE ~ as.numeric(.data[[source_prior]]),
        .use_prior %in% FALSE ~ as.numeric(.data[[source_next]]),
        TRUE ~ NA_real_
      ),
      !!matched_offset := case_when(
        .use_prior %in% TRUE ~ .prior_offset_hours,
        .use_prior %in% FALSE ~ .next_offset_hours,
        TRUE ~ NA_real_
      ),
      !!matched_direction := case_when(
        .data[[matched_offset]] < 0 ~ "prior",
        .data[[matched_offset]] == 0 ~ "exact",
        .data[[matched_offset]] > 0 ~ "future",
        TRUE ~ NA_character_
      ),
      .nearest_selection_valid =
        is.na(.data[[matched_offset]]) |
        (
          (is.na(.prior_offset_hours) |
             abs(.data[[matched_offset]]) <= abs(.prior_offset_hours) + 1e-12) &
            (is.na(.next_offset_hours) |
               abs(.data[[matched_offset]]) <= abs(.next_offset_hours) + 1e-12)
        ),
      .prior_tie_preference_valid =
        is.na(.prior_offset_hours) |
        is.na(.next_offset_hours) |
        abs(abs(.prior_offset_hours) - abs(.next_offset_hours)) > 1e-12 |
        .use_prior %in% TRUE
    )
  
  if (!isTRUE(all(result$.nearest_selection_valid))) {
    stop("The ", label, " nearest-value selection failed.")
  }
  
  if (!isTRUE(all(result$.prior_tie_preference_valid))) {
    stop("The ", label, " exact-distance tie did not prefer the prior value.")
  }
  
  result %>%
    arrange(.analysis_row_id) %>%
    select(
      -any_of(
        c(
          ".analysis_row_id", ".target_dttm", metric_prior, metric_next, source_prior, source_next,
          ".prior_offset_hours", ".next_offset_hours", ".use_prior", ".nearest_selection_valid",
          ".prior_tie_preference_valid", ".source_dttm", ".matched_source_dttm"
        )
      )
    )
}

analysis_with_sf <- join_metric_to_timepoints(
  analysis_data, sf_joined_collapsed, "sf_dttm", "sf_ratio", OXYGENATION_TIMEPOINT_MATCH_WINDOW_HOURS, "sf"
)

analysis_with_sf <- join_metric_to_timepoints(
  analysis_with_sf, pf_joined_collapsed, "pf_dttm", "pf_ratio", OXYGENATION_TIMEPOINT_MATCH_WINDOW_HOURS, "pf"
)

analysis_with_sf <- join_metric_to_timepoints(
  analysis_with_sf, oi_joined_collapsed, "oi_dttm", "oi", OXYGENATION_TIMEPOINT_MATCH_WINDOW_HOURS, "oi"
)

analysis_with_sf <- join_metric_to_timepoints(
  analysis_with_sf, osi_joined_collapsed, "osi_dttm", "osi", OXYGENATION_TIMEPOINT_MATCH_WINDOW_HOURS, "osi"
)

oxygenation_join_integrity_summary <- tibble(
  metric = c("SF ratio", "PF ratio", "OI", "OSI"),
  n_eligible_timepoints = nrow(analysis_with_sf),
  n_available = c(
    sum(!is.na(analysis_with_sf$sf_ratio)), sum(!is.na(analysis_with_sf$pf_ratio)),
    sum(!is.na(analysis_with_sf$oi)), sum(!is.na(analysis_with_sf$osi))
  ),
  net_row_count_change = nrow(analysis_with_sf) - nrow(analysis_data)
)

message("---- Oxygenation joins to IMV timepoints ----")
print(oxygenation_join_integrity_summary, n = Inf, width = Inf)

oxygenation_timepoint_match_summary <- bind_rows(
  analysis_with_sf %>%
    filter(!is.na(sf_ratio)) %>%
    transmute(
      metric = "SF", match_direction = sf_ratio_match_direction,
      match_offset_hours = sf_ratio_match_offset_hours
    ),
  analysis_with_sf %>%
    filter(!is.na(pf_ratio)) %>%
    transmute(
      metric = "PF", match_direction = pf_ratio_match_direction,
      match_offset_hours = pf_ratio_match_offset_hours
    ),
  analysis_with_sf %>%
    filter(!is.na(oi)) %>%
    transmute(
      metric = "OI", match_direction = oi_match_direction,
      match_offset_hours = oi_match_offset_hours
    ),
  analysis_with_sf %>%
    filter(!is.na(osi)) %>%
    transmute(
      metric = "OSI", match_direction = osi_match_direction,
      match_offset_hours = osi_match_offset_hours
    )
) %>%
  group_by(metric, match_direction) %>%
  summarise(
    n_matches = n(),
    minimum_offset_hours = min(match_offset_hours),
    median_offset_hours = median(match_offset_hours),
    maximum_offset_hours = max(match_offset_hours),
    .groups = "drop"
  )

oxygenation_timepoint_match_checks <- list(
  all_sf_offsets_within_window = all(
    abs(analysis_with_sf$sf_ratio_match_offset_hours) <= OXYGENATION_TIMEPOINT_MATCH_WINDOW_HOURS,
    na.rm = TRUE
  ),
  all_pf_offsets_within_window = all(
    abs(analysis_with_sf$pf_ratio_match_offset_hours) <= OXYGENATION_TIMEPOINT_MATCH_WINDOW_HOURS,
    na.rm = TRUE
  ),
  all_oi_offsets_within_window = all(
    abs(analysis_with_sf$oi_match_offset_hours) <= OXYGENATION_TIMEPOINT_MATCH_WINDOW_HOURS,
    na.rm = TRUE
  ),
  all_osi_offsets_within_window = all(
    abs(analysis_with_sf$osi_match_offset_hours) <= OXYGENATION_TIMEPOINT_MATCH_WINDOW_HOURS,
    na.rm = TRUE
  ),
  sf_value_and_offset_missingness_agree = identical(
    is.na(analysis_with_sf$sf_ratio),
    is.na(analysis_with_sf$sf_ratio_match_offset_hours)
  ),
  pf_value_and_offset_missingness_agree = identical(
    is.na(analysis_with_sf$pf_ratio),
    is.na(analysis_with_sf$pf_ratio_match_offset_hours)
  ),
  oi_value_and_offset_missingness_agree = identical(
    is.na(analysis_with_sf$oi),
    is.na(analysis_with_sf$oi_match_offset_hours)
  ),
  osi_value_and_offset_missingness_agree = identical(
    is.na(analysis_with_sf$osi),
    is.na(analysis_with_sf$osi_match_offset_hours)
  ),
  sf_value_and_timestamp_missingness_agree = identical(
    is.na(analysis_with_sf$sf_ratio),
    is.na(analysis_with_sf$sf_ratio_matched_dttm_num)
  ),
  pf_value_and_timestamp_missingness_agree = identical(
    is.na(analysis_with_sf$pf_ratio),
    is.na(analysis_with_sf$pf_ratio_matched_dttm_num)
  ),
  oi_value_and_timestamp_missingness_agree = identical(
    is.na(analysis_with_sf$oi),
    is.na(analysis_with_sf$oi_matched_dttm_num)
  ),
  osi_value_and_timestamp_missingness_agree = identical(
    is.na(analysis_with_sf$osi),
    is.na(analysis_with_sf$osi_matched_dttm_num)
  ),
  sf_offsets_reconstruct_source_time = all(
    is.na(analysis_with_sf$sf_ratio) |
      abs(
        analysis_with_sf$sf_ratio_matched_dttm_num -
          (
            as.numeric(analysis_with_sf$recorded_dttm) +
              analysis_with_sf$sf_ratio_match_offset_hours * 3600
          )
      ) < 1e-8
  ),
  pf_offsets_reconstruct_source_time = all(
    is.na(analysis_with_sf$pf_ratio) |
      abs(
        analysis_with_sf$pf_ratio_matched_dttm_num -
          (
            as.numeric(analysis_with_sf$recorded_dttm) +
              analysis_with_sf$pf_ratio_match_offset_hours * 3600
          )
      ) < 1e-8
  ),
  oi_offsets_reconstruct_source_time = all(
    is.na(analysis_with_sf$oi) |
      abs(
        analysis_with_sf$oi_matched_dttm_num -
          (
            as.numeric(analysis_with_sf$recorded_dttm) +
              analysis_with_sf$oi_match_offset_hours * 3600
          )
      ) < 1e-8
  ),
  osi_offsets_reconstruct_source_time = all(
    is.na(analysis_with_sf$osi) |
      abs(
        analysis_with_sf$osi_matched_dttm_num -
          (
            as.numeric(analysis_with_sf$recorded_dttm) +
              analysis_with_sf$osi_match_offset_hours * 3600
          )
      ) < 1e-8
  )
)

# Nonexact matches are sample characteristics, not requirements for valid matching.
# Attach this information to the saved matching summary without adding pass/fail checks.
oxygenation_nonexact_match_information <- tibble(
  metric = c("SF", "PF"),
  has_nonexact_matches = c(
    any(
      abs(analysis_with_sf$sf_ratio_match_offset_hours) > 1e-8,
      na.rm = TRUE
    ),
    any(
      abs(analysis_with_sf$pf_ratio_match_offset_hours) > 1e-8,
      na.rm = TRUE
    )
  )
)

attr(
  oxygenation_timepoint_match_summary,
  "nonexact_match_information"
) <- oxygenation_nonexact_match_information

oxygenation_timepoint_match_checks_table <- tibble(
  check = names(oxygenation_timepoint_match_checks),
  passed = as.logical(unlist(oxygenation_timepoint_match_checks, use.names = FALSE))
)

message("---- Oxygenation-to-timepoint nearest-match timing ----")
print(oxygenation_timepoint_match_summary, n = Inf, width = Inf)
print(oxygenation_timepoint_match_checks_table, n = Inf, width = Inf)

message("---- Informational: presence of nonexact SF and PF matches ----")
print(oxygenation_nonexact_match_information, n = Inf, width = Inf)

stopifnot(isTRUE(all(oxygenation_timepoint_match_checks_table$passed)))

# =============================================================================
# Join creatinine, platelets, pSOFA, and survival variables
# =============================================================================

join_source_to_timepoints <- function(analysis, source_data, source_time, value_columns, window_hours, label) {
  source <- source_data %>%
    select(hospitalization_id, all_of(source_time), all_of(value_columns)) %>%
    rename(.source_dttm = all_of(source_time)) %>%
    as.data.table()
  
  if (anyDuplicated(source[, .(hospitalization_id, .source_dttm)])) {
    stop(label, " source data contain duplicate hospitalization-time keys.")
  }
  
  setkey(source, hospitalization_id, .source_dttm)
  
  target <- analysis %>%
    mutate(.analysis_row_id = row_number(), .target_dttm = as.numeric(recorded_dttm)) %>%
    as.data.table()
  
  n_rows_before <- nrow(target)
  setkey(target, hospitalization_id, .target_dttm)
  
  joined <- source[
    target, on = .(hospitalization_id, .source_dttm = .target_dttm), roll = window_hours * 3600
  ]
  
  if (nrow(joined) != n_rows_before || anyDuplicated(joined$.analysis_row_id)) {
    stop("The ", label, " rolling join changed or duplicated IMV timepoints.")
  }
  
  as_tibble(joined) %>%
    arrange(.analysis_row_id) %>%
    select(-any_of(c(".analysis_row_id", ".source_dttm", ".target_dttm")))
}

# =============================================================================
# Creatinine
# =============================================================================

creat_labs <- cohort_labs %>%
  filter(lab_category == "creatinine", !is.na(lab_value_numeric)) %>%
  transmute(hospitalization_id, creat_dttm = as.numeric(lab_result_dttm), creatinine = lab_value_numeric)

creat_labs_collapsed <- creat_labs %>%
  group_by(hospitalization_id, creat_dttm) %>%
  summarise(creatinine = median(creatinine, na.rm = TRUE), .groups = "drop")

analysis_with_sf <- join_source_to_timepoints(
  analysis_with_sf, creat_labs_collapsed, "creat_dttm", "creatinine", 4, "creatinine"
) %>%
  group_by(hospitalization_id) %>%
  mutate(
    creatinine_baseline = {
      ordered_values <- creatinine[order(recorded_dttm)]
      available_values <- ordered_values[!is.na(ordered_values)]
      if (length(available_values) > 0L) available_values[[1]] else NA_real_
    },
    delta_creatinine = creatinine - creatinine_baseline
  ) %>%
  ungroup()

# =============================================================================
# Platelet count
# =============================================================================

platelet_labs <- cohort_labs %>%
  filter(lab_category == "platelet_count", !is.na(lab_value_numeric)) %>%
  transmute(hospitalization_id, platelet_dttm = as.numeric(lab_result_dttm), platelet_count = lab_value_numeric)

platelet_labs_collapsed <- platelet_labs %>%
  group_by(hospitalization_id, platelet_dttm) %>%
  summarise(platelet_count = median(platelet_count, na.rm = TRUE), .groups = "drop")

analysis_with_sf <- join_source_to_timepoints(
  analysis_with_sf, platelet_labs_collapsed, "platelet_dttm", "platelet_count", 4, "platelet"
) %>%
  group_by(hospitalization_id) %>%
  mutate(
    platelet_baseline = {
      ordered_values <- platelet_count[order(recorded_dttm)]
      available_values <- ordered_values[!is.na(ordered_values)]
      if (length(available_values) > 0L) available_values[[1]] else NA_real_
    },
    delta_platelet = platelet_count - platelet_baseline
  ) %>%
  ungroup()

rm(join_source_to_timepoints)

# =============================================================================
# Daily pSOFA and survival variables
# =============================================================================

analysis_with_sf <- analysis_with_sf %>%
  left_join(
    cohort_demographics %>% select(hospitalization_id, admission_dttm, discharge_dttm, death_dttm),
    by = "hospitalization_id",
    relationship = "many-to-one"
  ) %>%
  mutate(
    hours_from_admission = as.numeric(difftime(recorded_dttm, admission_dttm, units = "hours")),
    sofa_day_raw = as.integer(floor(hours_from_admission / 24)),
    in_psofa_preadmission_grace =
      !is.na(hours_from_admission) &
      hours_from_admission >= -PSOFA_PREADMISSION_GRACE_HOURS &
      hours_from_admission < 0,
    sofa_day = if_else(in_psofa_preadmission_grace, 0L, sofa_day_raw)
  ) %>%
  left_join(
    sofa_daily %>%
      select(
        hospitalization_id, sofa_day, sofa_total, sofa_cv, sofa_coag, sofa_liver, sofa_resp, sofa_cns,
        sofa_renal
      ),
    by = c("hospitalization_id", "sofa_day"),
    relationship = "many-to-one"
  ) %>%
  mutate(
    psofa_match_source = case_when(
      in_psofa_preadmission_grace & !is.na(sofa_total) ~ "Day 0 via pre-admission grace",
      !is.na(sofa_total) ~ "Same admission-relative day",
      TRUE ~ NA_character_
    ),
    death_day = as.numeric(difftime(death_dttm, admission_dttm, units = "days")),
    mortality_event_60 = if_else(!is.na(death_day) & death_day >= 0 & death_day <= 60, 1L, 0L),
    surv_time = if_else(mortality_event_60 == 1L, death_day, 60)
  )

if (nrow(analysis_with_sf) != nrow(analysis_data)) {
  stop("Joining pSOFA and survival variables unexpectedly changed the number of IMV timepoints.")
}

psofa_timepoint_match_summary <- analysis_with_sf %>%
  count(psofa_match_source, name = "n_timepoints", .drop = FALSE) %>%
  mutate(percent = round(100 * n_timepoints / sum(n_timepoints), 1))

psofa_grace_summary <- analysis_with_sf %>%
  summarise(
    preadmission_grace_hours = PSOFA_PREADMISSION_GRACE_HOURS,
    maximum_sofa_day = SOFA_MAX_DAY,
    n_timepoints = n(),
    n_preadmission_timepoints = sum(hours_from_admission < 0, na.rm = TRUE),
    n_timepoints_in_grace = sum(in_psofa_preadmission_grace),
    n_grace_timepoints_with_psofa = sum(
      in_psofa_preadmission_grace & !is.na(sofa_total)
    ),
    n_timepoints_before_grace = sum(
      !is.na(hours_from_admission) & hours_from_admission < -PSOFA_PREADMISSION_GRACE_HOURS
    ),
    n_timepoints_after_maximum_day = sum(sofa_day_raw > SOFA_MAX_DAY, na.rm = TRUE)
  )

psofa_grace_checks <- list(
  grace_rows_mapped_to_day_zero =
    all(analysis_with_sf$sofa_day[analysis_with_sf$in_psofa_preadmission_grace] == 0L),
  non_grace_rows_preserve_raw_day = all(
    analysis_with_sf$in_psofa_preadmission_grace |
      analysis_with_sf$sofa_day == analysis_with_sf$sofa_day_raw
  ),
  no_psofa_before_grace = all(
    is.na(analysis_with_sf$sofa_total) |
      is.na(analysis_with_sf$hours_from_admission) |
      analysis_with_sf$hours_from_admission >= -PSOFA_PREADMISSION_GRACE_HOURS
  ),
  no_psofa_after_maximum_day = all(
    is.na(analysis_with_sf$sofa_total) |
      is.na(analysis_with_sf$sofa_day_raw) |
      analysis_with_sf$sofa_day_raw <= SOFA_MAX_DAY
  ),
  psofa_match_source_complete_when_score_present = all(
    is.na(analysis_with_sf$sofa_total) | !is.na(analysis_with_sf$psofa_match_source)
  )
)

psofa_grace_checks_table <- tibble(
  check = names(psofa_grace_checks),
  passed = as.logical(unlist(psofa_grace_checks, use.names = FALSE))
)

message("---- pSOFA matching with pre-admission grace ----")
print(psofa_timepoint_match_summary, n = Inf, width = Inf)
print(psofa_grace_summary, width = Inf)
print(psofa_grace_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(psofa_grace_checks_table$passed)))

# =============================================================================
# 3e. Select the cross-sectional cohort
# =============================================================================
# A hospitalization is eligible if at least one post-intubation timepoint:
#   1. has VT/PBW, VT/PFVC, same-day or grace-linked pSOFA, and at least one oxygenation metric; and
#   2. meets any hypoxemia criterion: OI >= 4, OSI >= 5, PF <= 300, or SF <= 250.
#
# VT/PBW remains an analytic variable but is not restricted to 4-8 mL/kg.
# Eligible pre-admission timepoints within five hours of admission use day-0 pSOFA.
# Completed oxygenation metrics are matched independently using the closest value within ±5 hours.

analysis_with_completeness <- analysis_with_sf %>%
  mutate(
    is_hypoxemic = case_when(
      is.na(oi) & is.na(osi) & is.na(pf_ratio) & is.na(sf_ratio) ~ NA,
      !is.na(oi) & oi >= OI_HYPOXEMIA_THRESHOLD ~ TRUE,
      !is.na(osi) & osi >= OSI_HYPOXEMIA_THRESHOLD ~ TRUE,
      !is.na(pf_ratio) & pf_ratio <= PF_HYPOXEMIA_THRESHOLD ~ TRUE,
      !is.na(sf_ratio) & sf_ratio <= SF_HYPOXEMIA_THRESHOLD ~ TRUE,
      TRUE ~ FALSE
    ),
    has_oxygenation = !is.na(oi) | !is.na(osi) | !is.na(pf_ratio) | !is.na(sf_ratio),
    has_all_data = !is.na(vtpbw) & !is.na(vtpfvc) & has_oxygenation & !is.na(sofa_total)
  )

completeness_timepoint_summary <- tibble(
  measure = c(
    "Eligible post-intubation timepoints entering completeness gate", "Timepoints with VT/PBW",
    "Timepoints with VT/PFVC", "Timepoints with at least one oxygenation metric",
    "Timepoints with eligible pSOFA", "Timepoints satisfying all completeness requirements"
  ),
  n_timepoints = c(
    nrow(analysis_with_completeness), sum(!is.na(analysis_with_completeness$vtpbw)),
    sum(!is.na(analysis_with_completeness$vtpfvc)), sum(analysis_with_completeness$has_oxygenation),
    sum(!is.na(analysis_with_completeness$sofa_total)), sum(analysis_with_completeness$has_all_data)
  )
) %>%
  mutate(percent_of_eligible_timepoints = round(100 * n_timepoints / first(n_timepoints), 1))

completeness_by_hospitalization <- analysis_with_completeness %>%
  group_by(hospitalization_id) %>%
  summarise(
    n_eligible_timepoints = n(),
    n_complete_timepoints = sum(has_all_data),
    ever_has_vtpbw = any(!is.na(vtpbw)),
    ever_has_vtpfvc = any(!is.na(vtpfvc)),
    ever_has_oxygenation = any(has_oxygenation),
    ever_has_sofa = any(!is.na(sofa_total)),
    has_complete_timepoint = any(has_all_data),
    .groups = "drop"
  ) %>%
  mutate(
    never_has_vtpbw = !ever_has_vtpbw,
    never_has_vtpfvc = !ever_has_vtpfvc,
    never_has_oxygenation = !ever_has_oxygenation,
    never_has_sofa = !ever_has_sofa,
    components_never_coincide = !has_complete_timepoint & ever_has_vtpbw & ever_has_vtpfvc &
      ever_has_oxygenation & ever_has_sofa
  )

completeness_retention_summary <- completeness_by_hospitalization %>%
  count(has_complete_timepoint, name = "n_hospitalizations") %>%
  mutate(
    completeness_status = if_else(
      has_complete_timepoint,
      "Retained; at least one complete eligible invasive timepoint",
      "Excluded; no complete eligible invasive timepoint"
    ),
    percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1)
  ) %>%
  select(completeness_status, n_hospitalizations, percent)

completeness_failure_summary <- tibble(
  failure_reason = c(
    "Never has VT/PBW", "Never has VT/PFVC", "Never has an oxygenation metric",
    "Never has pSOFA within retained timing bounds", "All components occur, but never at the same timepoint"
  ),
  n_hospitalizations = c(
    sum(completeness_by_hospitalization$never_has_vtpbw),
    sum(completeness_by_hospitalization$never_has_vtpfvc),
    sum(completeness_by_hospitalization$never_has_oxygenation),
    sum(completeness_by_hospitalization$never_has_sofa),
    sum(completeness_by_hospitalization$components_never_coincide)
  )
)

message("---- Complete-data gate: timepoint availability ----")
print(completeness_timepoint_summary, n = Inf, width = Inf)

message("---- Complete-data gate: hospitalization retention ----")
print(completeness_retention_summary, n = Inf, width = Inf)

message("---- Complete-data gate: nonexclusive failure reasons ----")
print(completeness_failure_summary, n = Inf, width = Inf)

write_parquet(analysis_with_completeness, file.path(output_dir, "analysis_all_eligible_timepoints.parquet"))

# =============================================================================
# QC: hypoxemia gate
# =============================================================================

hypoxemia_gate_timepoints <- analysis_with_completeness %>%
  filter(has_all_data) %>%
  mutate(
    meets_oi = !is.na(oi) & oi >= OI_HYPOXEMIA_THRESHOLD,
    meets_osi = !is.na(osi) & osi >= OSI_HYPOXEMIA_THRESHOLD,
    meets_pf = !is.na(pf_ratio) & pf_ratio <= PF_HYPOXEMIA_THRESHOLD,
    meets_sf = !is.na(sf_ratio) & sf_ratio <= SF_HYPOXEMIA_THRESHOLD,
    n_hypoxemia_criteria_met = as.integer(meets_oi) + as.integer(meets_osi) + as.integer(meets_pf) +
      as.integer(meets_sf)
  )

hypoxemia_gate_timepoint_summary <- tibble(
  measure = c(
    "Complete eligible invasive timepoints entering hypoxemia gate",
    "Timepoints meeting at least one hypoxemia criterion",
    "Timepoints meeting no hypoxemia criterion",
    "Timepoints meeting more than one hypoxemia criterion"
  ),
  n_timepoints = c(
    nrow(hypoxemia_gate_timepoints),
    sum(hypoxemia_gate_timepoints$is_hypoxemic == TRUE, na.rm = TRUE),
    sum(hypoxemia_gate_timepoints$is_hypoxemic == FALSE, na.rm = TRUE),
    sum(hypoxemia_gate_timepoints$n_hypoxemia_criteria_met > 1, na.rm = TRUE)
  )
) %>%
  mutate(percent_of_entering_timepoints = round(100 * n_timepoints / first(n_timepoints), 1))

hypoxemia_criterion_summary <- tibble(
  criterion = c("OI >= 4", "OSI >= 5", "PF ratio <= 300", "SF ratio <= 250", "Any hypoxemia criterion"),
  n_timepoints = c(
    sum(hypoxemia_gate_timepoints$meets_oi), sum(hypoxemia_gate_timepoints$meets_osi),
    sum(hypoxemia_gate_timepoints$meets_pf), sum(hypoxemia_gate_timepoints$meets_sf),
    sum(hypoxemia_gate_timepoints$is_hypoxemic == TRUE, na.rm = TRUE)
  ),
  n_hospitalizations = c(
    n_distinct(hypoxemia_gate_timepoints$hospitalization_id[hypoxemia_gate_timepoints$meets_oi]),
    n_distinct(hypoxemia_gate_timepoints$hospitalization_id[hypoxemia_gate_timepoints$meets_osi]),
    n_distinct(hypoxemia_gate_timepoints$hospitalization_id[hypoxemia_gate_timepoints$meets_pf]),
    n_distinct(hypoxemia_gate_timepoints$hospitalization_id[hypoxemia_gate_timepoints$meets_sf]),
    n_distinct(hypoxemia_gate_timepoints$hospitalization_id[hypoxemia_gate_timepoints$is_hypoxemic == TRUE])
  )
)

hypoxemia_gate_by_hospitalization <- analysis_with_completeness %>%
  group_by(hospitalization_id) %>%
  summarise(
    n_complete_timepoints = sum(has_all_data),
    n_qualifying_hypoxemic_timepoints = sum(has_all_data & is_hypoxemic == TRUE, na.rm = TRUE),
    has_complete_timepoint = any(has_all_data),
    has_qualifying_hypoxemic_timepoint = any(has_all_data & is_hypoxemic == TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hypoxemia_gate_status = case_when(
      !has_complete_timepoint ~ "Did not enter; no complete eligible invasive timepoint",
      has_qualifying_hypoxemic_timepoint ~
        "Retained; at least one complete eligible invasive timepoint met a hypoxemia criterion",
      TRUE ~ "Excluded; no complete eligible invasive timepoint met a hypoxemia criterion"
    )
  )

hypoxemia_gate_retention_summary <- hypoxemia_gate_by_hospitalization %>%
  filter(has_complete_timepoint) %>%
  count(hypoxemia_gate_status, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1)) %>%
  arrange(desc(n_hospitalizations))

message("---- Hypoxemia gate: timepoint summary ----")
print(hypoxemia_gate_timepoint_summary, n = Inf, width = Inf)

message("---- Hypoxemia gate: overlapping criterion counts ----")
print(hypoxemia_criterion_summary, n = Inf, width = Inf)

message("---- Hypoxemia gate: hospitalization retention ----")
print(hypoxemia_gate_retention_summary, n = Inf, width = Inf)

# =============================================================================
# Select one index timepoint per eligible hospitalization
# =============================================================================

qualifying_timepoints <- analysis_with_completeness %>%
  filter(has_all_data, is_hypoxemic == TRUE)

INDEX_WINDOW_HOURS <- 6

eligible_start <- analysis_with_completeness %>%
  group_by(hospitalization_id) %>%
  summarise(eligible_start_dttm = min(recorded_dttm), .groups = "drop")

index_tier1 <- qualifying_timepoints %>%
  left_join(eligible_start, by = "hospitalization_id", relationship = "many-to-one") %>%
  filter(
    !is.na(dp), recorded_dttm >= eligible_start_dttm,
    recorded_dttm <= eligible_start_dttm + hours(INDEX_WINDOW_HOURS)
  ) %>%
  group_by(hospitalization_id) %>%
  slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-eligible_start_dttm)

index_tier2 <- qualifying_timepoints %>%
  filter(!hospitalization_id %in% index_tier1$hospitalization_id) %>%
  group_by(hospitalization_id) %>%
  slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>%
  ungroup()

selected_index_timepoints <- bind_rows(
  index_tier1 %>% mutate(index_selection_tier = "Tier 1: pressure-complete within 6 hours"),
  index_tier2 %>% mutate(index_selection_tier = "Tier 2: first qualifying timepoint")
)

first_qualifying_timepoint <- qualifying_timepoints %>%
  group_by(hospitalization_id) %>%
  summarise(first_qualifying_dttm = min(recorded_dttm), .groups = "drop")

first_tier1_candidate <- qualifying_timepoints %>%
  left_join(eligible_start, by = "hospitalization_id", relationship = "many-to-one") %>%
  filter(
    !is.na(dp), recorded_dttm >= eligible_start_dttm,
    recorded_dttm <= eligible_start_dttm + hours(INDEX_WINDOW_HOURS)
  ) %>%
  group_by(hospitalization_id) %>%
  summarise(first_tier1_candidate_dttm = min(recorded_dttm), .groups = "drop")

index_selection_by_hospitalization <- selected_index_timepoints %>%
  left_join(eligible_start, by = "hospitalization_id", relationship = "many-to-one") %>%
  left_join(first_qualifying_timepoint, by = "hospitalization_id", relationship = "many-to-one") %>%
  left_join(first_tier1_candidate, by = "hospitalization_id", relationship = "many-to-one") %>%
  mutate(
    hours_from_eligible_start = as.numeric(difftime(recorded_dttm, eligible_start_dttm, units = "hours")),
    selected_row_is_eligible = has_all_data & is_hypoxemic == TRUE,
    tier1_requirements_met = index_selection_tier == "Tier 1: pressure-complete within 6 hours" &
      !is.na(dp) & !is.na(hours_from_eligible_start) & hours_from_eligible_start >= 0 &
      hours_from_eligible_start <= INDEX_WINDOW_HOURS & recorded_dttm == first_tier1_candidate_dttm,
    tier2_requirements_met = index_selection_tier == "Tier 2: first qualifying timepoint" &
      is.na(first_tier1_candidate_dttm) & recorded_dttm == first_qualifying_dttm,
    index_rule_satisfied = case_when(
      index_selection_tier == "Tier 1: pressure-complete within 6 hours" ~ tier1_requirements_met,
      index_selection_tier == "Tier 2: first qualifying timepoint" ~ tier2_requirements_met,
      TRUE ~ FALSE
    )
  )

qualifying_hospitalization_ids <- qualifying_timepoints %>%
  distinct(hospitalization_id)

selected_hospitalization_ids <- selected_index_timepoints %>%
  distinct(hospitalization_id)

index_selection_integrity_summary <- tibble(
  measure = c("Hospitalizations with at least one qualifying timepoint", "Selected index rows",
              "Unique hospitalizations with a selected index row", "Duplicate selected hospitalization IDs",
              "Qualifying hospitalizations missing an index row",
              "Selected hospitalizations not in qualifying cohort",
              "Selected rows failing eligibility requirements", "Selected rows failing the two-tier index rule"),
  n = c(nrow(qualifying_hospitalization_ids), nrow(selected_index_timepoints),
        n_distinct(selected_index_timepoints$hospitalization_id),
        sum(duplicated(selected_index_timepoints$hospitalization_id)),
        nrow(anti_join(qualifying_hospitalization_ids, selected_hospitalization_ids, by = "hospitalization_id")),
        nrow(anti_join(selected_hospitalization_ids, qualifying_hospitalization_ids, by = "hospitalization_id")),
        sum(!index_selection_by_hospitalization$selected_row_is_eligible, na.rm = TRUE),
        sum(!index_selection_by_hospitalization$index_rule_satisfied, na.rm = TRUE))
)

index_selection_tier_summary <- index_selection_by_hospitalization %>%
  group_by(index_selection_tier) %>%
  summarise(
    n_hospitalizations = n_distinct(hospitalization_id),
    n_with_driving_pressure = sum(!is.na(dp)),
    median_hours_from_eligible_start = median(hours_from_eligible_start, na.rm = TRUE),
    maximum_hours_from_eligible_start = max(hours_from_eligible_start, na.rm = TRUE),
    n_failing_index_rule = sum(!index_rule_satisfied, na.rm = TRUE),
    .groups = "drop"
  )

message("---- Index-timepoint selection integrity ----")
print(index_selection_integrity_summary, n = Inf, width = Inf)

message("---- Index-timepoint selection by tier ----")
print(index_selection_tier_summary, n = Inf, width = Inf)

index_error_measures <- c("Duplicate selected hospitalization IDs",
                          "Qualifying hospitalizations missing an index row",
                          "Selected hospitalizations not in qualifying cohort",
                          "Selected rows failing eligibility requirements",
                          "Selected rows failing the two-tier index rule")

if (any(index_selection_integrity_summary$n[index_selection_integrity_summary$measure %in% index_error_measures] > 0)) {
  stop("The two-tier index-timepoint selection failed an integrity check.")
}

if (anyDuplicated(cohort_weights$hospitalization_id)) {
  stop("cohort_weights contains duplicate hospitalization IDs.")
}

cross_sectional <- selected_index_timepoints %>%
  left_join(cohort_weights, by = "hospitalization_id", relationship = "one-to-one") %>%
  mutate(
    bmi = if_else(!is.na(weight_kg) & height_cm > 0, weight_kg / (height_cm / 100)^2, NA_real_)
  )

abw_weight_join_matches <- (
  is.na(cross_sectional$abw) & is.na(cross_sectional$weight_kg)
) | (
  !is.na(cross_sectional$abw) & !is.na(cross_sectional$weight_kg) &
    abs(cross_sectional$abw - cross_sectional$weight_kg) <= 1e-12
)

if (!all(abw_weight_join_matches)) {
  stop("The ABW carried from pbw_pfvc_data does not match weight_kg after the cross-sectional weight join.")
}

missing_cross_sectional_provenance_columns <- setdiff(required_height_provenance_columns, names(cross_sectional))

if (length(missing_cross_sectional_provenance_columns) > 0L) {
  stop(
    "cross_sectional is missing required height-provenance column(s): ",
    paste(missing_cross_sectional_provenance_columns, collapse = ", ")
  )
}

cross_sectional_height_provenance_summary <- cross_sectional %>%
  count(height_source, name = "n_hospitalizations", .drop = FALSE) %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

cross_sectional_height_provenance_checks <- list(
  required_columns_present = TRUE,
  no_missing_height_source = !any(is.na(cross_sectional$height_source)),
  donor_rows_have_gap = all(
    cross_sectional$height_source != "Another encounter within window" | !is.na(cross_sectional$height_gap_days)
  ),
  non_donor_rows_have_zero_gap = all(
    cross_sectional$height_source == "Another encounter within window" |
      (!is.na(cross_sectional$height_gap_days) & cross_sectional$height_gap_days == 0)
  )
)

cross_sectional_height_provenance_checks_table <- tibble(
  check = names(cross_sectional_height_provenance_checks),
  passed = as.logical(unlist(cross_sectional_height_provenance_checks, use.names = FALSE))
)

message("---- Cross-sectional height provenance ----")
print(cross_sectional_height_provenance_summary, n = Inf, width = Inf)
print(cross_sectional_height_provenance_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(cross_sectional_height_provenance_checks_table$passed)))

eligible_patients <- cross_sectional$hospitalization_id

message(
  "Index selection: ", nrow(index_tier1), " Tier 1 hospitalizations and ", nrow(index_tier2),
  " Tier 2 hospitalizations"
)

message(
  "Index timepoints with driving pressure: ", sum(!is.na(cross_sectional$dp)), " of ", nrow(cross_sectional)
)

analysis_all <- analysis_with_sf %>%
  filter(hospitalization_id %in% eligible_patients)

# =============================================================================
# 3f. Derive 28-day ventilator-free days
# =============================================================================
# VFD-28 requires successful final liberation from IMV before day 28.
#
# Rules:
#   * death within 28 days -> 0 VFDs;
#   * ongoing or unresolved IMV at the observation endpoint -> 0 VFDs;
#   * successful final liberation -> days from final liberation through day 28.
#
# Respiratory-support rows represent state intervals rather than fixed one-hour observations. Each state is carried until
# the next waterfall timestamp, with a maximum carry-forward interval of one hour.

vfd_demographic_times <- cohort_demographics %>%
  transmute(
    hospitalization_id, admission_dttm = as.POSIXct(admission_dttm, tz = "UTC"),
    discharge_dttm = as.POSIXct(discharge_dttm, tz = "UTC"), death_dttm = as.POSIXct(death_dttm, tz = "UTC")
  )

vfd_duplicate_timestamp_summary <- resp_waterfall %>%
  filter(hospitalization_id %in% eligible_patients) %>%
  count(hospitalization_id, recorded_dttm, name = "n_rows") %>%
  filter(n_rows > 1)

if (nrow(vfd_duplicate_timestamp_summary) > 0L) {
  stop("Duplicate respiratory-waterfall timestamps prevent unambiguous VFD interval construction.")
}

vfd_support_intervals <- resp_waterfall %>%
  filter(hospitalization_id %in% eligible_patients) %>%
  inner_join(vfd_demographic_times, by = "hospitalization_id", relationship = "many-to-one") %>%
  mutate(
    day28_dttm = admission_dttm + days(28),
    observation_end_dttm = if_else(
      !is.na(discharge_dttm) & discharge_dttm >= admission_dttm & discharge_dttm < day28_dttm,
      discharge_dttm,
      day28_dttm
    ),
    is_imv = coalesce(
      str_to_lower(str_trim(as.character(device_category))) == "imv" |
        str_to_lower(str_trim(as.character(support))) == "invasive",
      FALSE
    )
  ) %>%
  filter(!is.na(recorded_dttm), recorded_dttm >= admission_dttm, recorded_dttm < observation_end_dttm) %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  group_by(hospitalization_id) %>%
  mutate(
    next_waterfall_dttm = lead(recorded_dttm),
    interval_end_dttm = pmin(
      coalesce(next_waterfall_dttm, observation_end_dttm), recorded_dttm + hours(1), observation_end_dttm
    ),
    interval_hours = as.numeric(difftime(interval_end_dttm, recorded_dttm, units = "hours"))
  ) %>%
  ungroup()

if (n_distinct(vfd_support_intervals$hospitalization_id) != n_distinct(eligible_patients)) {
  stop("At least one final-cohort hospitalization has no respiratory-support interval available for VFD derivation.")
}

vfd_interval_integrity_summary <- tibble(
  measure = c("Support intervals", "Intervals with missing duration", "Intervals with nonpositive duration",
              "Intervals longer than one hour"),
  n = c(nrow(vfd_support_intervals), sum(is.na(vfd_support_intervals$interval_hours)),
        sum(vfd_support_intervals$interval_hours <= 0, na.rm = TRUE),
        sum(vfd_support_intervals$interval_hours > 1, na.rm = TRUE))
)

vfd_interval_error_measures <- c("Intervals with missing duration", "Intervals with nonpositive duration",
                                 "Intervals longer than one hour")

if (any(vfd_interval_integrity_summary$n[
  vfd_interval_integrity_summary$measure %in% vfd_interval_error_measures
] > 0)) {
  stop("The VFD respiratory-support intervals failed an integrity check.")
}

# Identify contiguous IMV episodes.
vfd_imv_episodes <- vfd_support_intervals %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  group_by(hospitalization_id) %>%
  mutate(
    starts_imv_episode = is_imv & !lag(is_imv, default = FALSE),
    imv_episode_number = cumsum(starts_imv_episode)
  ) %>%
  ungroup() %>%
  filter(is_imv) %>%
  group_by(hospitalization_id, imv_episode_number) %>%
  summarise(
    episode_start_dttm = min(recorded_dttm), episode_end_dttm = max(interval_end_dttm), .groups = "drop"
  )

# Determine final respiratory status and the end of the final IMV episode.
vfd_final_support_status <- vfd_support_intervals %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  group_by(hospitalization_id) %>%
  summarise(
    observation_end_dttm = first(observation_end_dttm),
    observation_ended_at_discharge = !is.na(first(discharge_dttm)) &
      first(observation_end_dttm) == first(discharge_dttm),
    observation_ended_at_day28 = first(observation_end_dttm) == first(day28_dttm),
    final_documented_is_imv = last(is_imv),
    .groups = "drop"
  ) %>%
  left_join(
    vfd_imv_episodes %>%
      group_by(hospitalization_id) %>%
      summarise(final_imv_end_dttm = max(episode_end_dttm), .groups = "drop"),
    by = "hospitalization_id",
    relationship = "one-to-one"
  )

if (nrow(vfd_final_support_status) != n_distinct(eligible_patients) ||
    anyDuplicated(vfd_final_support_status$hospitalization_id)) {
  stop("Final respiratory-support status did not produce exactly one row per eligible hospitalization.")
}

vfd_data_detailed <- vfd_final_support_status %>%
  left_join(vfd_demographic_times, by = "hospitalization_id", relationship = "one-to-one") %>%
  mutate(
    death_day = as.numeric(difftime(death_dttm, admission_dttm, units = "days")),
    discharge_day = as.numeric(difftime(discharge_dttm, admission_dttm, units = "days")),
    final_liberation_day = as.numeric(difftime(final_imv_end_dttm, admission_dttm, units = "days")),
    died_within_28 = !is.na(death_day) & death_day >= 0 & death_day <= 28,
    successful_final_liberation = !died_within_28 & !final_documented_is_imv &
      !is.na(final_liberation_day) & final_liberation_day >= 0 & final_liberation_day < 28,
    vfd_28 = case_when(
      died_within_28 ~ 0,
      successful_final_liberation ~ pmax(28 - final_liberation_day, 0),
      TRUE ~ 0
    ),
    # 0 = no documented successful liberation; 1 = successful final liberation; 2 = death within 28 days.
    vfd_status = case_when(
      died_within_28 ~ 2L,
      successful_final_liberation ~ 1L,
      TRUE ~ 0L
    ),
    vfd_time = case_when(
      vfd_status == 2L ~ death_day,
      vfd_status == 1L ~ final_liberation_day,
      observation_ended_at_discharge & !is.na(discharge_day) ~ pmin(discharge_day, 28),
      TRUE ~ 28
    ),
    vfd_time = pmax(vfd_time, 1 / 24)
  )

vfd_data <- vfd_data_detailed %>%
  select(hospitalization_id, vfd_28, vfd_time, vfd_status)

cross_sectional <- cross_sectional %>%
  left_join(vfd_data, by = "hospitalization_id", relationship = "one-to-one")

vfd_final_integrity_summary <- tibble(
  measure = c("Cross-sectional rows after VFD join", "Duplicate hospitalization IDs after VFD join",
              "Rows with missing VFD-28", "Rows with missing VFD time", "Rows with missing VFD status",
              "VFD-28 values outside 0-28", "VFD times outside 0-28", "Deaths with nonzero VFD-28",
              "No-liberation cases with nonzero VFD-28", "Liberation events inconsistent with VFD-28",
              "Unexpected VFD status codes"),
  n = c(nrow(cross_sectional), sum(duplicated(cross_sectional$hospitalization_id)),
        sum(is.na(cross_sectional$vfd_28)), sum(is.na(cross_sectional$vfd_time)),
        sum(is.na(cross_sectional$vfd_status)),
        sum(cross_sectional$vfd_28 < 0 | cross_sectional$vfd_28 > 28, na.rm = TRUE),
        sum(cross_sectional$vfd_time < 0 | cross_sectional$vfd_time > 28, na.rm = TRUE),
        sum(cross_sectional$vfd_status == 2L & cross_sectional$vfd_28 != 0, na.rm = TRUE),
        sum(cross_sectional$vfd_status == 0L & cross_sectional$vfd_28 != 0, na.rm = TRUE),
        sum(cross_sectional$vfd_status == 1L &
              abs(cross_sectional$vfd_28 - (28 - cross_sectional$vfd_time)) > 1e-8, na.rm = TRUE),
        sum(!cross_sectional$vfd_status %in% c(0L, 1L, 2L)))
)

vfd_error_measures <- c("Duplicate hospitalization IDs after VFD join", "Rows with missing VFD-28",
                        "Rows with missing VFD time", "Rows with missing VFD status", "VFD-28 values outside 0-28",
                        "VFD times outside 0-28", "Deaths with nonzero VFD-28",
                        "No-liberation cases with nonzero VFD-28",
                        "Liberation events inconsistent with VFD-28", "Unexpected VFD status codes")

if (any(vfd_final_integrity_summary$n[
  vfd_final_integrity_summary$measure %in% vfd_error_measures
] > 0)) {
  stop("The final VFD-28 derivation failed an integrity check.")
}

vfd_outcome_status_summary <- cross_sectional %>%
  mutate(
    vfd_status_label = case_when(
      vfd_status == 0L ~ "No documented successful liberation",
      vfd_status == 1L ~ "Successful final liberation",
      vfd_status == 2L ~ "Death within 28 days",
      TRUE ~ "Unexpected"
    )
  ) %>%
  count(vfd_status, vfd_status_label, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

message("---- VFD interval integrity ----")
print(vfd_interval_integrity_summary, n = Inf, width = Inf)

message("---- Final VFD-28 integrity ----")
print(vfd_final_integrity_summary, n = Inf, width = Inf)

message("---- Final VFD-28 outcome status ----")
print(vfd_outcome_status_summary, n = Inf, width = Inf)

message("Final-liberation-based median VFD-28: ", round(median(cross_sectional$vfd_28, na.rm = TRUE), 1))

# =============================================================================
# 3g. VFR terciles
# =============================================================================

cross_sectional <- cross_sectional %>%
  mutate(
    vfr_tercile = ntile(vtpfvc, 3),
    vfr_tercile = factor(vfr_tercile, labels = c("T1 (Low)", "T2 (Mid)", "T3 (High)"))
  )

# =============================================================================
# QC: final-output integrity
# =============================================================================

final_qualifying_ids <- qualifying_timepoints %>%
  distinct(hospitalization_id)

final_cross_sectional_ids <- cross_sectional %>%
  distinct(hospitalization_id)

final_analysis_all_ids <- analysis_all %>%
  distinct(hospitalization_id)

final_output_integrity_summary <- tibble(
  measure = c("Cross-sectional rows", "Unique cross-sectional hospitalization IDs",
              "Duplicate cross-sectional hospitalization IDs", "Missing cross-sectional hospitalization IDs",
              "Missing cross-sectional index timestamps",
              "Qualifying hospitalizations missing from cross-sectional output",
              "Cross-sectional hospitalizations outside qualifying cohort", "Selected rows without complete data",
              "Selected rows with invalid tidal-volume source", "Selected rows with missing or nonpositive tidal-volume value",
              "Selected rows with missing intubation-identification method", "Selected rows with invalid intubation-identification method",
              "Selected rows not classified as hypoxemic", "Cross-sectional rows with missing VFD-28",
              "Cross-sectional rows with missing VFR tercile", "Unique hospitalizations in all-timepoints output",
              "All-timepoints rows with missing intubation-identification method",
              "All-timepoints rows with invalid intubation-identification method",
              "Final hospitalizations missing from all-timepoints output",
              "All-timepoints hospitalizations outside final cohort"),
  n = c(nrow(cross_sectional), n_distinct(cross_sectional$hospitalization_id),
        sum(duplicated(cross_sectional$hospitalization_id)),
        sum(is.na(cross_sectional$hospitalization_id) |
              str_trim(as.character(cross_sectional$hospitalization_id)) == ""),
        sum(is.na(cross_sectional$recorded_dttm)),
        nrow(anti_join(final_qualifying_ids, final_cross_sectional_ids, by = "hospitalization_id")),
        nrow(anti_join(final_cross_sectional_ids, final_qualifying_ids, by = "hospitalization_id")),
        sum(!cross_sectional$has_all_data, na.rm = TRUE),
        sum(is.na(cross_sectional$tidal_volume_source) | !cross_sectional$tidal_volume_source %in% c("set", "observed")),
        sum(is.na(cross_sectional$tidal_volume_value) | cross_sectional$tidal_volume_value <= 0),
        sum(is.na(cross_sectional$intubation_identification_method)),
        sum(!is.na(cross_sectional$intubation_identification_method) & !cross_sectional$intubation_identification_method %in%
            allowed_intubation_identification_methods),
        sum(is.na(cross_sectional$is_hypoxemic) | cross_sectional$is_hypoxemic != TRUE),
        sum(is.na(cross_sectional$vfd_28)), sum(is.na(cross_sectional$vfr_tercile)),
        n_distinct(analysis_all$hospitalization_id),
        sum(is.na(analysis_all$intubation_identification_method)),
        sum(!is.na(analysis_all$intubation_identification_method) & !analysis_all$intubation_identification_method %in%
            allowed_intubation_identification_methods),
        nrow(anti_join(final_cross_sectional_ids, final_analysis_all_ids, by = "hospitalization_id")),
        nrow(anti_join(final_analysis_all_ids, final_cross_sectional_ids, by = "hospitalization_id")))
)

final_output_error_measures <- c("Duplicate cross-sectional hospitalization IDs",
                                 "Missing cross-sectional hospitalization IDs",
                                 "Missing cross-sectional index timestamps",
                                 "Qualifying hospitalizations missing from cross-sectional output",
                                 "Cross-sectional hospitalizations outside qualifying cohort",
                                 "Selected rows without complete data", "Selected rows with invalid tidal-volume source",
                                 "Selected rows with missing or nonpositive tidal-volume value",
                                 "Selected rows with missing intubation-identification method",
                                 "Selected rows with invalid intubation-identification method",
                                 "Selected rows not classified as hypoxemic",
                                 "Cross-sectional rows with missing VFD-28",
                                 "Cross-sectional rows with missing VFR tercile",
                                 "All-timepoints rows with missing intubation-identification method",
                                 "All-timepoints rows with invalid intubation-identification method",
                                 "Final hospitalizations missing from all-timepoints output",
                                 "All-timepoints hospitalizations outside final cohort")

if (any(final_output_integrity_summary$n[
  final_output_integrity_summary$measure %in% final_output_error_measures
] > 0)) {
  stop("The final analysis output failed at least one integrity check.")
}

final_output_variable_availability <- tibble(
  variable = c(
    "PBW", "PFVC", "PBW/PFVC", "VT/PBW", "VT/PFVC", "SOFA total", "Tidal-volume value", "Tidal-volume source", "Intubation-identification method", 
    "At least one oxygenation metric", "Driving pressure", "Height source", "Height gap days", "VFD-28", "VFR tercile"
  ),
  n_available = c(
    sum(!is.na(cross_sectional$pbw)), sum(!is.na(cross_sectional$pfvc)),
    sum(!is.na(cross_sectional$pbwpfvc)), sum(!is.na(cross_sectional$vtpbw)),
    sum(!is.na(cross_sectional$vtpfvc)), sum(!is.na(cross_sectional$sofa_total)),
    sum(!is.na(cross_sectional$tidal_volume_value)), sum(!is.na(cross_sectional$tidal_volume_source)),
    sum(!is.na(cross_sectional$intubation_identification_method)),
    sum(!is.na(cross_sectional$oi) | !is.na(cross_sectional$osi) | !is.na(cross_sectional$pf_ratio) | !is.na(cross_sectional$sf_ratio)),
    sum(!is.na(cross_sectional$dp)), sum(!is.na(cross_sectional$height_source)),
    sum(!is.na(cross_sectional$height_gap_days)), sum(!is.na(cross_sectional$vfd_28)),
    sum(!is.na(cross_sectional$vfr_tercile))
  )
) %>%
  mutate(
    n_missing = nrow(cross_sectional) - n_available,
    percent_available = round(100 * n_available / nrow(cross_sectional), 1)
  )

final_vfr_tercile_summary <- cross_sectional %>%
  count(vfr_tercile, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

final_intubation_method_summary <- cross_sectional %>%
  count(intubation_identification_method, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

message("---- Final cross-sectional output integrity ----")
print(final_output_integrity_summary, n = Inf, width = Inf)

message("---- Final cross-sectional variable availability ----")
print(final_output_variable_availability, n = Inf, width = Inf)

message("---- Final VFR tercile distribution ----")
print(final_vfr_tercile_summary, n = Inf, width = Inf)

message("---- Final intubation-identification method distribution ----")
print(final_intubation_method_summary, n = Inf, width = Inf)

# =============================================================================
# Reconcile the complete Script 03 filtering funnel
# =============================================================================

script03_filter_funnel <- tibble(
  filter_step = c(
    "PBW and positive PFVC available",
    "At least one complete eligible invasive timepoint",
    "At least one complete eligible invasive timepoint meeting any hypoxemia criterion"
  ),
  n_hospitalizations = c(
    n_distinct(pbw_pfvc_data$hospitalization_id),
    sum(completeness_by_hospitalization$has_complete_timepoint),
    n_distinct(cross_sectional$hospitalization_id)
  )
) %>%
  mutate(
    n_removed_at_step = lag(n_hospitalizations) - n_hospitalizations,
    percent_retained_from_previous_step = round(100 * n_hospitalizations / lag(n_hospitalizations), 1),
    percent_retained_from_pbw_pfvc = round(100 * n_hospitalizations / first(n_hospitalizations), 1)
  )

script03_filter_reconciliation <- tibble(
  checkpoint = c("PBW/PFVC cohort", "Complete eligible timepoint cohort", "Final hypoxemic cohort"),
  funnel_count = script03_filter_funnel$n_hospitalizations,
  independently_calculated_count = c(
    nrow(pbw_pfvc_data),
    n_distinct(analysis_with_completeness$hospitalization_id[analysis_with_completeness$has_all_data]),
    n_distinct(qualifying_timepoints$hospitalization_id)
  )
) %>%
  mutate(
    counts_agree = funnel_count == independently_calculated_count,
    difference = funnel_count - independently_calculated_count
  )

message("---- Complete Script 03 filtering funnel ----")
print(script03_filter_funnel, n = Inf, width = Inf)

message("---- Script 03 funnel reconciliation ----")
print(script03_filter_reconciliation, n = Inf, width = Inf)

if (any(!script03_filter_reconciliation$counts_agree)) {
  stop("At least one Script 03 filtering count failed reconciliation.")
}

# =============================================================================
# 3h. Save analysis outputs
# =============================================================================

write_parquet(analysis_all, file.path(output_dir, "analysis_all_timepoints.parquet"))
write_parquet(cross_sectional, file.path(output_dir, "analysis_cross_sectional.parquet"))
write_parquet(sofa_scores, file.path(output_dir, "sofa_scores.parquet"))
write_parquet(sofa_daily, file.path(output_dir, "sofa_daily.parquet"))

required_analysis_outputs <- c("analysis_broad_pfvc.parquet", "analysis_all_eligible_timepoints.parquet",
                               "analysis_all_timepoints.parquet", "analysis_cross_sectional.parquet",
                               "sofa_scores.parquet", "sofa_daily.parquet")

missing_analysis_outputs <- required_analysis_outputs[
  !file.exists(file.path(output_dir, required_analysis_outputs))
]

if (length(missing_analysis_outputs) > 0L) {
  stop("Script 03 failed to create required output(s): ", paste(missing_analysis_outputs, collapse = ", "))
}

# =============================================================================
# 3i. Complete the attrition log
# =============================================================================

n_step4 <- n_distinct(pbw_pfvc_data$hospitalization_id)
n_step5 <- n_distinct(
  analysis_with_completeness$hospitalization_id[analysis_with_completeness$has_all_data]
)
n_step6 <- n_distinct(eligible_patients)

partial_path <- file.path(output_dir, "attrition_log_partial.csv")

if (!file.exists(partial_path)) {
  stop("Missing attrition_log_partial.csv from Script 01: ", partial_path)
}

attrition <- read_csv(partial_path, show_col_types = FALSE) %>%
  attrition_add(
    ATTRITION_STEPS[4], n_step4,
    exclusion_reason = "Not selected as the latest eligible hospitalization per patient, height outside 80-200 cm, or PBW/PFVC unavailable"
  ) %>%
  attrition_add(
    ATTRITION_STEPS[5], n_step5,
    exclusion_reason = "No complete eligible invasive timepoint with required index data"
  ) %>%
  attrition_add(
    ATTRITION_STEPS[6], n_step6,
    exclusion_reason = "No complete eligible invasive timepoint meeting a hypoxemia criterion"
  ) %>%
  mutate(site = site_name, .before = 1)

expected_attrition_exclusions <- head(attrition$n_remaining, -1L) - tail(attrition$n_remaining, -1L)

attrition_checks <- list(
  attrition_has_six_steps = nrow(attrition) == length(ATTRITION_STEPS) && length(ATTRITION_STEPS) == 6L,
  attrition_step_order_valid = identical(
    as.integer(attrition$step_order),
    seq_along(ATTRITION_STEPS)
  ),
  attrition_labels_match_protocol = identical(
    as.character(attrition$step_label),
    as.character(ATTRITION_STEPS)
  ),
  attrition_counts_complete = !any(is.na(attrition$n_remaining)),
  attrition_counts_are_whole_numbers = all(
    attrition$n_remaining == as.integer(attrition$n_remaining)
  ),
  attrition_counts_nonnegative = all(attrition$n_remaining >= 0),
  attrition_counts_nonincreasing = all(diff(attrition$n_remaining) <= 0),
  attrition_first_exclusion_is_missing = is.na(attrition$n_excluded[[1]]),
  attrition_remaining_exclusions_complete = !any(is.na(attrition$n_excluded[-1])),
  attrition_exclusions_are_whole_numbers = all(
    attrition$n_excluded[-1] == as.integer(attrition$n_excluded[-1])
  ),
  attrition_exclusions_nonnegative = all(attrition$n_excluded[-1] >= 0),
  attrition_exclusions_reconcile = isTRUE(all.equal(
    as.numeric(attrition$n_excluded[-1]),
    as.numeric(expected_attrition_exclusions),
    tolerance = 0,
    check.attributes = FALSE
  )),
  attrition_final_count_matches_cross_sectional = last(attrition$n_remaining) == nrow(cross_sectional),
  attrition_final_count_matches_unique_hospitalizations = last(attrition$n_remaining) ==
    n_distinct(cross_sectional$hospitalization_id)
)

attrition_check_qc <- tibble(
  check = names(attrition_checks),
  passed = as.logical(unlist(attrition_checks, use.names = FALSE))
)

print(attrition, n = Inf, width = Inf)
print(attrition_check_qc, n = Inf, width = Inf)

stopifnot(
  isTRUE(all(attrition_check_qc$passed))
)

# Shared CSV deliverables belong in the configured final/tables directory.
table_dir <- config$table_dir

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(table_dir)) {
  stop("Script 03 could not create the configured final table directory: ", table_dir)
}

attrition_output_path <- file.path(table_dir, paste0("attrition_log_", site_name, ".csv"))

write_csv(attrition, attrition_output_path)

attrition_readback <- read_csv(attrition_output_path, show_col_types = FALSE)

attrition_readback_matches <- identical(names(attrition_readback), names(attrition)) &&
  nrow(attrition_readback) == nrow(attrition) &&
  identical(as.character(attrition_readback$site), as.character(attrition$site)) &&
  identical(as.integer(attrition_readback$step_order), as.integer(attrition$step_order)) &&
  identical(as.character(attrition_readback$step_label), as.character(attrition$step_label)) &&
  identical(as.integer(attrition_readback$n_remaining), as.integer(attrition$n_remaining)) &&
  identical(as.integer(attrition_readback$n_excluded), as.integer(attrition$n_excluded)) &&
  identical(as.character(attrition_readback$exclusion_reason), as.character(attrition$exclusion_reason))

stopifnot(
  file.exists(attrition_output_path),
  attrition_readback_matches
)

message(
  "Hospitalization-level attrition log written: ",
  paste(attrition$n_remaining, collapse = " -> ")
)

# =============================================================================
# 3j. Federated PBW:PFVC distribution exports
# =============================================================================

dist_groups <- bind_rows(
  cross_sectional %>%
    transmute(group_type = "sex", group_value = as.character(sex_category), value = pbwpfvc),
  cross_sectional %>%
    transmute(group_type = "race", group_value = as.character(race_category), value = pbwpfvc),
  cross_sectional %>%
    mutate(group_value = as.character(cut(age_at_admission, AGE_BIN_EDGES, right = FALSE))) %>%
    transmute(group_type = "age_bin", group_value, value = pbwpfvc),
  cross_sectional %>%
    mutate(group_value = as.character(cut(height_cm, HEIGHT_BIN_EDGES, right = FALSE))) %>%
    transmute(group_type = "height_bin", group_value, value = pbwpfvc)
) %>%
  filter(!is.na(value), !is.na(group_value), group_value != "NA")

dist_keep <- dist_groups %>%
  count(group_type, group_value) %>%
  filter(n >= DIST_MIN_CELL)

dist_groups <- dist_groups %>%
  semi_join(dist_keep, by = c("group_type", "group_value"))

bin_minimum <- min(PBWPFVC_BIN_EDGES)
bin_maximum <- max(PBWPFVC_BIN_EDGES)

dist_histograms <- dist_groups %>%
  mutate(
    value = pmin(pmax(value, bin_minimum), bin_maximum - 1e-9),
    bin_index = findInterval(value, PBWPFVC_BIN_EDGES, rightmost.closed = TRUE)
  ) %>%
  count(group_type, group_value, bin_index, name = "count") %>%
  mutate(bin_left = PBWPFVC_BIN_EDGES[bin_index], bin_right = PBWPFVC_BIN_EDGES[bin_index + 1]) %>%
  transmute(site = site_name, group_type, group_value, bin_left, bin_right, count) %>%
  arrange(group_type, group_value, bin_left)

dist_quantiles <- dist_groups %>%
  group_by(group_type, group_value) %>%
  summarise(
    n = n(), median = median(value), q1 = quantile(value, 0.25), q3 = quantile(value, 0.75),
    p10 = quantile(value, 0.10), p20 = quantile(value, 0.20), p30 = quantile(value, 0.30),
    p40 = quantile(value, 0.40), p50 = quantile(value, 0.50), p60 = quantile(value, 0.60),
    p70 = quantile(value, 0.70), p80 = quantile(value, 0.80), p90 = quantile(value, 0.90),
    .groups = "drop"
  ) %>%
  mutate(site = site_name, .before = 1)

write_csv(dist_histograms, file.path(table_dir, paste0("dist_histograms_", site_name, ".csv")))
write_csv(dist_quantiles, file.path(table_dir, paste0("dist_quantiles_", site_name, ".csv")))

message("Federated distribution exports written for ", nrow(dist_keep), " demographic groups")

# =============================================================================
# Append concise Script 03 QC to the combined diagnostic RDS
# =============================================================================
# Use the configured cohort diagnostic path established at startup.
# Reload the latest saved object below before appending Script 03 diagnostics.

if (!file.exists(diagnostic_rds_path)) {
  stop("The combined Script 01/02 diagnostic RDS was not found: ", diagnostic_rds_path)
}

oxygenation_qc_passed <- all(
  oxygenation_metric_validity_summary$n_nonfinite == 0 &
    oxygenation_metric_validity_summary$n_nonpositive == 0
) &&
  all(oxygenation_join_integrity_summary$net_row_count_change == 0) &&
  all(oxygenation_timepoint_match_checks_table$passed)

psofa_grace_qc_passed <- all(psofa_grace_checks_table$passed)

index_qc_passed <- all(
  index_selection_integrity_summary$n[
    index_selection_integrity_summary$measure %in% index_error_measures
  ] == 0
)

vfd_interval_qc_passed <- all(
  vfd_interval_integrity_summary$n[
    vfd_interval_integrity_summary$measure %in% vfd_interval_error_measures
  ] == 0
)

vfd_qc_passed <- vfd_interval_qc_passed && all(
  vfd_final_integrity_summary$n[vfd_final_integrity_summary$measure %in% vfd_error_measures] == 0
)

final_output_qc_passed <- all(
  final_output_integrity_summary$n[
    final_output_integrity_summary$measure %in% final_output_error_measures
  ] == 0
)

funnel_qc_passed <- all(script03_filter_reconciliation$counts_agree)
dp_qc_passed <- isTRUE(all(unlist(dp_derivation_checks)))
height_provenance_qc_passed <- all(
  height_provenance_input_checks_table$passed,
  pbw_pfvc_height_provenance_checks_table$passed,
  cross_sectional_height_provenance_checks_table$passed
)

overall_script03_qc_passed <- all(
  oxygenation_qc_passed, psofa_grace_qc_passed, index_qc_passed, dp_qc_passed,
  height_provenance_qc_passed, vfd_qc_passed, final_output_qc_passed, funnel_qc_passed,
  isTRUE(all(selected_tv_cleaning_checks_table$passed))
)

script03_qc_summary <- tibble(
  qc_area = c(
    "Cohort funnel", "Post-intubation timepoints", "Complete data", "Hypoxemia",
    "Oxygenation calculations", "pSOFA respiratory device", "Index selection", "Driving pressure",
    "Height provenance", "VFD-28", "Final output", "Overall"
  ),
  result = c(
    paste0(
      paste(script03_filter_funnel$n_hospitalizations, collapse = " -> "),
      " hospitalizations across the Script 03 funnel."
    ),
    paste0(
      nrow(analysis_data), " eligible post-intubation rows across ",
      n_distinct(analysis_data$hospitalization_id), " hospitalizations; ",
      sum(analysis_data$tidal_volume_source == "set"), " set and ",
      sum(analysis_data$tidal_volume_source == "observed"),
      " observed. Intubation-identification provenance was complete and valid."
    ),
    paste0(
      sum(completeness_by_hospitalization$has_complete_timepoint),
      " hospitalizations had at least one complete eligible invasive timepoint."
    ),
    paste0(
      n_distinct(qualifying_timepoints$hospitalization_id),
      " hospitalizations had a complete eligible invasive timepoint meeting a hypoxemia criterion."
    ),
    "SF, PF, OI, and OSI were finite and positive; all joins preserved eligible timepoint row counts.",
    paste0(
      sum(resp_daily_selection$device_selection_changed, na.rm = TRUE),
      " patient-days required severity-based device selection; ",
      psofa_grace_summary$n_grace_timepoints_with_psofa,
      " eligible timepoints received day-0 pSOFA through the ",
      PSOFA_PREADMISSION_GRACE_HOURS, "-hour pre-admission grace."
    ),
    paste0(nrow(index_tier1), " Tier 1 and ", nrow(index_tier2), " Tier 2 selections; no rule violations."),
    paste0(
      dp_derivation_qc$n_dp_valid, " valid driving-pressure estimates among ",
      dp_derivation_qc$n_timepoints, " eligible timepoints; ",
      dp_derivation_qc$n_plateau_derived, " plateau-derived and ",
      dp_derivation_qc$n_pip_derived, " PIP-derived. ",
      dp_derivation_qc$n_nonpositive_dp_raw, " nonpositive candidate values were invalidated. ",
      sum(!is.na(cross_sectional$dp)), " of ", nrow(cross_sectional),
      " selected index rows had valid driving pressure."
    ),
    paste0(
      nrow(cross_sectional), " final hospitalizations retained height source and gap provenance; ",
      sum(cross_sectional$height_source == "Another encounter within window"),
      " used an in-window donor height."
    ),
    paste0(
      sum(cross_sectional$vfd_status == 1L), " successful final liberations; ",
      sum(cross_sectional$vfd_status == 0L), " without documented liberation; ",
      sum(cross_sectional$vfd_status == 2L), " deaths. Median VFD-28 = ",
      round(median(cross_sectional$vfd_28, na.rm = TRUE), 1), " days."
    ),
    paste0(
      nrow(cross_sectional), " rows represented ",
      n_distinct(cross_sectional$hospitalization_id), " unique final hospitalizations."
    ),
    if_else(
      overall_script03_qc_passed,
      "All retained Script 03 QC checks passed.",
      "At least one retained Script 03 QC check requires review."
    )
  ),
  status = c(
    "INFORMATIONAL", "INFORMATIONAL", "INFORMATIONAL", "INFORMATIONAL",
    if_else(oxygenation_qc_passed, "PASS", "REVIEW"),
    if_else(psofa_grace_qc_passed, "PASS", "REVIEW"),
    if_else(index_qc_passed, "PASS", "REVIEW"),
    if_else(dp_qc_passed, "PASS", "REVIEW"),
    if_else(height_provenance_qc_passed, "PASS", "REVIEW"),
    if_else(vfd_qc_passed, "PASS", "REVIEW"),
    if_else(final_output_qc_passed, "PASS", "REVIEW"),
    if_else(overall_script03_qc_passed, "PASS", "REVIEW")
  )
)

# =============================================================================
# QC: height-bin coverage before distribution-export group filtering
# =============================================================================

height_bin_coverage_detail <- cross_sectional %>%
  transmute(
    hospitalization_id,
    height_cm,
    height_bin = as.character(cut(height_cm, breaks = HEIGHT_BIN_EDGES, right = FALSE)),
    analytic_height_valid = is.finite(height_cm) & height_cm >= 80 & height_cm <= 200,
    height_bin_assigned = !is.na(height_bin)
  )

height_bin_coverage_summary <- height_bin_coverage_detail %>%
  summarise(
    n_hospitalizations = n(),
    n_valid_analytic_heights = sum(analytic_height_valid),
    n_assigned_to_height_bin = sum(height_bin_assigned),
    n_unassigned = sum(!height_bin_assigned),
    n_in_lower_extension = sum(height_cm >= 80 & height_cm < 92, na.rm = TRUE),
    n_in_upper_extension = sum(height_cm >= 182 & height_cm <= 200, na.rm = TRUE)
  )

height_bin_counts_before_filtering <- height_bin_coverage_detail %>%
  count(height_bin, name = "n_hospitalizations", .drop = FALSE)

height_bin_coverage_checks <- list(
  one_row_per_final_hospitalization =
    nrow(height_bin_coverage_detail) == nrow(cross_sectional) &&
    anyDuplicated(height_bin_coverage_detail$hospitalization_id) == 0L,
  all_final_heights_within_analytic_bounds =
    all(height_bin_coverage_detail$analytic_height_valid),
  all_final_heights_assigned_to_bins =
    all(height_bin_coverage_detail$height_bin_assigned),
  height_bin_counts_reconcile =
    sum(height_bin_counts_before_filtering$n_hospitalizations) == nrow(cross_sectional),
  accepted_endpoints_receive_bins =
    !anyNA(cut(c(80, 200), breaks = HEIGHT_BIN_EDGES, right = FALSE))
)

height_bin_coverage_checks_table <- tibble(
  check = names(height_bin_coverage_checks),
  passed = vapply(height_bin_coverage_checks, isTRUE, logical(1))
)

if (!all(height_bin_coverage_checks_table$passed)) {
  stop(
    "Height-bin coverage QC failed: ",
    paste(
      height_bin_coverage_checks_table$check[!height_bin_coverage_checks_table$passed],
      collapse = ", "
    )
  )
}

height_bin_coverage_qc <- list(
  bin_edges_cm = HEIGHT_BIN_EDGES,
  interval_rule = "Lower boundary inclusive; upper boundary exclusive",
  population = "Final cross-sectional cohort before distribution-export group filtering",
  interpretation = paste0(
    "Coverage checks validate height-group assignment only. They do not establish that every group is exported ",
    "or that the distribution exports satisfy an approved disclosure policy."
  ),
  summary = height_bin_coverage_summary,
  bin_counts_before_filtering = height_bin_counts_before_filtering,
  hospitalization_level = height_bin_coverage_detail,
  checks = height_bin_coverage_checks_table
)

message(
  "Height-bin coverage QC passed: ",
  height_bin_coverage_summary$n_assigned_to_height_bin,
  " of ", height_bin_coverage_summary$n_hospitalizations,
  " final hospitalizations assigned."
)

script03_compact_qc <- list(
  metadata = tibble(
    script = "Script 03: Variable Derivation", site = site_name, generated_at = as.character(Sys.time()),
    overall_status = if_else(overall_script03_qc_passed, "PASS", "REVIEW"),
    note = "Report-ready QC with hospitalization-level selected tidal-volume cleaning diagnostics."
  ),
  summary = script03_qc_summary,
  cohort_funnel = script03_filter_funnel,
  funnel_reconciliation = script03_filter_reconciliation,
  post_intubation_timepoints = list(
    allowed_identification_methods = allowed_intubation_identification_methods,
    input_summary = post_intubation_input_summary,
    method_summary = post_intubation_method_summary,
    method_event_consistency = post_intubation_method_event_consistency,
    input_checks = post_intubation_input_checks_table,
    context_summary = eligible_timepoint_context_qc,
    context_source_summary = eligible_timepoint_context_source_summary,
    context_method_summary = eligible_timepoint_method_summary,
    context_checks = eligible_timepoint_context_checks_table,
    context_preservation_checks = waterfall_context_preservation_checks,
    final_method_summary = final_intubation_method_summary,
    selected_tidal_volume_cleaning = list(
      policy = "Clean the originally selected value; preserve rows and source; no source substitution",
      thresholds = selected_tv_thresholds,
      summary = selected_tv_cleaning_summary,
      hospitalization_level = selected_tv_cleaning_by_hospitalization,
      retention = selected_tv_cleaning_retention,
      checks = selected_tv_cleaning_checks_table
    )
  ),
  oxygenation_timepoint_matching = list(
    component_match_window_hours = OXYGENATION_COMPONENT_MATCH_WINDOW_HOURS,
    timepoint_match_window_hours = OXYGENATION_TIMEPOINT_MATCH_WINDOW_HOURS,
    summary = oxygenation_timepoint_match_summary,
    checks = oxygenation_timepoint_match_checks_table
  ),
  psofa_timepoint_matching = list(
    maximum_day = SOFA_MAX_DAY,
    preadmission_grace_hours = PSOFA_PREADMISSION_GRACE_HOURS,
    source_summary = psofa_timepoint_match_summary,
    grace_summary = psofa_grace_summary,
    checks = psofa_grace_checks_table
  ),
  attrition_log = list(
    unit = "hospitalization",
    data = attrition,
    checks = tibble(
      check = names(attrition_checks),
      passed = as.logical(unlist(attrition_checks, use.names = FALSE))
    )
  ),
  respiratory_mechanics = list(
    driving_pressure_summary = dp_derivation_qc,
    driving_pressure_sources = dp_source_qc,
    driving_pressure_checks = tibble(
      check = names(dp_derivation_checks),
      passed = as.logical(unlist(dp_derivation_checks, use.names = FALSE))
    )
  ),
  height_provenance = list(
    required_columns = required_height_provenance_columns,
    input_summary = height_provenance_input_summary,
    input_checks = height_provenance_input_checks_table,
    pbw_pfvc_summary = pbw_pfvc_height_provenance_summary,
    pbw_pfvc_checks = pbw_pfvc_height_provenance_checks_table,
    cross_sectional_summary = cross_sectional_height_provenance_summary,
    cross_sectional_checks = cross_sectional_height_provenance_checks_table
  ),
  vfd_outcome_status = vfd_outcome_status_summary,
  final_variable_availability = final_output_variable_availability
)

script03_compact_qc$psofa_respiratory_derivation <- psofa_respiratory_derivation_qc

script03_compact_qc$height_bin_coverage <- height_bin_coverage_qc

cohort_diagnostics <- readRDS(diagnostic_rds_path)

if (!is.list(cohort_diagnostics)) {
  stop("The combined diagnostic RDS does not contain a list.")
}

required_existing_sections <- c("quality_checks_script02")
missing_existing_sections <- setdiff(required_existing_sections, names(cohort_diagnostics))

if (length(missing_existing_sections) > 0L) {
  stop(
    "The combined diagnostic RDS is missing required prior section(s): ",
    paste(missing_existing_sections, collapse = ", "),
    ". Complete Script 02 before executing Script 03."
  )
}

cohort_diagnostics$quality_checks_script03 <- script03_compact_qc
saveRDS(cohort_diagnostics, diagnostic_rds_path)

cohort_diagnostics_check <- readRDS(diagnostic_rds_path)

if (!identical(cohort_diagnostics_check, cohort_diagnostics)) {
  stop("The combined diagnostic RDS failed its complete Script 03 read-back comparison.")
}

script03_qc_readback_checks <- list(
  height_bin_coverage_matches = identical(
    cohort_diagnostics_check$quality_checks_script03$height_bin_coverage,
    script03_compact_qc$height_bin_coverage
  ),
  height_bin_coverage_checks_passed = isTRUE(all(
    cohort_diagnostics_check$quality_checks_script03$height_bin_coverage$checks$passed
  )),
  psofa_respiratory_derivation_matches = identical(
    cohort_diagnostics_check$quality_checks_script03$psofa_respiratory_derivation,
    script03_compact_qc$psofa_respiratory_derivation
  ),
  psofa_respiratory_derivation_checks_passed = isTRUE(all(
    cohort_diagnostics_check$quality_checks_script03$psofa_respiratory_derivation$checks$passed
  )),
  cohort_funnel_matches = identical(
    cohort_diagnostics_check$quality_checks_script03$cohort_funnel,
    script03_compact_qc$cohort_funnel
  ),
  driving_pressure_summary_matches = identical(
    cohort_diagnostics_check$quality_checks_script03$respiratory_mechanics$driving_pressure_summary,
    script03_compact_qc$respiratory_mechanics$driving_pressure_summary
  ),
  driving_pressure_sources_match = identical(
    cohort_diagnostics_check$quality_checks_script03$respiratory_mechanics$driving_pressure_sources,
    script03_compact_qc$respiratory_mechanics$driving_pressure_sources
  ),
  driving_pressure_checks_match = identical(
    cohort_diagnostics_check$quality_checks_script03$respiratory_mechanics$driving_pressure_checks,
    script03_compact_qc$respiratory_mechanics$driving_pressure_checks
  ),
  height_provenance_matches = identical(
    cohort_diagnostics_check$quality_checks_script03$height_provenance,
    script03_compact_qc$height_provenance
  ),
  post_intubation_timepoints_match = identical(
    cohort_diagnostics_check$quality_checks_script03$post_intubation_timepoints,
    script03_compact_qc$post_intubation_timepoints
  ),
  oxygenation_timepoint_matching_matches = identical(
    cohort_diagnostics_check$quality_checks_script03$oxygenation_timepoint_matching,
    script03_compact_qc$oxygenation_timepoint_matching
  ),
  psofa_timepoint_matching_matches = identical(
    cohort_diagnostics_check$quality_checks_script03$psofa_timepoint_matching,
    script03_compact_qc$psofa_timepoint_matching
  )
)

if (!isTRUE(all(unlist(script03_qc_readback_checks)))) {
  stop("The Script 03 QC summary failed its RDS read-back check.")
}

message("Concise Script 03 QC summary added to: ", diagnostic_rds_path)
message("Script 03 QC status: ", cohort_diagnostics_check$quality_checks_script03$metadata$overall_status)

rm(cohort_diagnostics_check)

# =============================================================================
# Script 03 completion
# =============================================================================
# Diagnostic report rendering is a separate workflow step.
# The report must consume the saved diagnostic RDS without refitting models.
# Script 03 does not render, copy, or delete report files.

message("Script 03 complete.")
message("  Cross-sectional cohort: ", nrow(cross_sectional), " hospitalizations")
message("  All retained timepoints: ", nrow(analysis_all))
message("  Analysis intermediates: ", output_dir)
message("  CSV deliverables: ", table_dir)
message("  Combined diagnostic RDS: ", diagnostic_rds_path)
message("  Diagnostic PDF rendering is deferred to the separate reporting workflow.")