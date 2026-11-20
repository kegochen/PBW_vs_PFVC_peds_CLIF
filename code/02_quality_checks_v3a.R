# =============================================================================
# Script 02: Quality Checks
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================

rm(list = ls())

library(tidyverse)
library(arrow)
library(lubridate)

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
thresholds <- config$thresholds_dir
diagnostic_dir <- config$cohort_diagnostic_dir
diagnostic_rds_path <- file.path(diagnostic_dir, "cohort_diagnostics.rds")

# Script 02 consumes Script 01 outputs. Do not create empty prerequisite directories
# that could conceal an incorrect output location or a missing Script 01 run.

required_script02_directories <- c(output_dir, diagnostic_dir, thresholds)
missing_script02_directories <- required_script02_directories[
  !dir.exists(required_script02_directories)
]

if (length(missing_script02_directories) > 0L) {
  stop(
    "Script 02 prerequisite directory(s) not found: ",
    paste(missing_script02_directories, collapse = "; "),
    ". Confirm configuration and complete Script 01 before executing Script 02."
  )
}

message("Script 02 starting | site: ", site_name)
  
# =============================================================================
# Utility functions
# =============================================================================
  
  replace_outliers_with_na_long <- function(df, threshold_table, category_variable, numeric_variable) {
    df %>%
      left_join(threshold_table, by = category_variable, relationship = "many-to-one") %>%
      mutate(
        !!sym(numeric_variable) := if_else(
          !is.na(lower_limit) & !is.na(upper_limit) &
            (.data[[numeric_variable]] < lower_limit | .data[[numeric_variable]] > upper_limit),
          NA_real_,
          as.numeric(.data[[numeric_variable]])
        )
      ) %>%
      select(-lower_limit, -upper_limit)
  }
  
  generate_summary_stats <- function(data, category_variable, numeric_variable) {
    data %>%
      group_by({{ category_variable }}) %>%
      summarise(
        N = sum(!is.na({{ numeric_variable }})),
        Min = if_else(N > 0, min({{ numeric_variable }}, na.rm = TRUE), NA_real_),
        Max = if_else(N > 0, max({{ numeric_variable }}, na.rm = TRUE), NA_real_),
        Mean = if_else(N > 0, mean({{ numeric_variable }}, na.rm = TRUE), NA_real_),
        Median = if_else(N > 0, median({{ numeric_variable }}, na.rm = TRUE), NA_real_),
        Q1 = if_else(N > 0, as.numeric(quantile({{ numeric_variable }}, 0.25, na.rm = TRUE)), NA_real_),
        Q3 = if_else(N > 0, as.numeric(quantile({{ numeric_variable }}, 0.75, na.rm = TRUE)), NA_real_),
        .groups = "drop"
      ) %>%
      arrange({{ category_variable }})
  }
  
# =============================================================================
# Load intermediate data from Script 01
# =============================================================================
  
  required_input_files <- c(
    "cohort_hospitalization_ids.rds", "resp_support_waterfall.parquet", "post_intubation_tidal_volume_rows.parquet",
    "cohort_vitals.parquet", "cohort_labs.parquet", "cohort_heights.parquet", "cohort_weights.parquet"
  )
  
  missing_input_files <- required_input_files[!file.exists(file.path(output_dir, required_input_files))]
  
  if (length(missing_input_files) > 0L) {
    stop("Script 02 is missing required input file(s): ", paste(missing_input_files, collapse = ", "))
  }
  
  if (!file.exists(diagnostic_rds_path)) {
    stop("The Script 01 diagnostic RDS was not found: ", diagnostic_rds_path)
  }
  
  message("Script 02 | loading and validating Script 01 inputs")
  cohort_ids <- readRDS(file.path(output_dir, "cohort_hospitalization_ids.rds"))
  resp_waterfall <- read_parquet(file.path(output_dir, "resp_support_waterfall.parquet"))
  post_intubation_tidal_volume_rows <- read_parquet(file.path(output_dir, "post_intubation_tidal_volume_rows.parquet"))
  cohort_vitals <- read_parquet(file.path(output_dir, "cohort_vitals.parquet"))
  cohort_labs <- read_parquet(file.path(output_dir, "cohort_labs.parquet"))
  cohort_heights <- read_parquet(file.path(output_dir, "cohort_heights.parquet"))
  cohort_weights <- read_parquet(file.path(output_dir, "cohort_weights.parquet"))
  
  required_post_intubation_columns <- c(
    ".resp_source_row_id", "intubation_event_id", "hospitalization_id", "intubation_dttm",
    "intubation_window_end_dttm", "recorded_dttm", "support_resolved",
    "intubation_identification_method", "tidal_volume_source", "tidal_volume_value"
  )
  
  missing_post_intubation_columns <- setdiff(required_post_intubation_columns, names(post_intubation_tidal_volume_rows))
  
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
  
  post_intubation_input_qc <- post_intubation_tidal_volume_rows %>%
    summarise(
      n_rows = n(),
      n_hospitalizations = n_distinct(hospitalization_id),
      n_intubation_events = n_distinct(intubation_event_id),
      n_set_rows = sum(tidal_volume_source == "set", na.rm = TRUE),
      n_observed_rows = sum(tidal_volume_source == "observed", na.rm = TRUE),
      n_rows_outside_cohort = sum(!hospitalization_id %in% cohort_ids),
      n_missing_hospitalization_id = sum(is.na(hospitalization_id)),
      n_missing_event_id = sum(is.na(intubation_event_id)),
      n_missing_intubation_dttm = sum(is.na(intubation_dttm)),
      n_missing_window_end_dttm = sum(is.na(intubation_window_end_dttm)),
      n_missing_recorded_dttm = sum(is.na(recorded_dttm)),
      n_missing_identification_method = sum(is.na(intubation_identification_method)),
      n_invalid_identification_method = sum(
        !is.na(intubation_identification_method) &
          !intubation_identification_method %in% allowed_intubation_identification_methods
      ),
      n_noninvasive_rows = sum(is.na(support_resolved) | support_resolved != "invasive"),
      n_invalid_sources = sum(is.na(tidal_volume_source) | !tidal_volume_source %in% c("set", "observed")),
      n_missing_or_nonpositive_values = sum(is.na(tidal_volume_value) | tidal_volume_value <= 0),
      n_duplicate_source_rows = sum(duplicated(.resp_source_row_id)),
      n_rows_before_intubation = sum(recorded_dttm < intubation_dttm, na.rm = TRUE),
      n_rows_after_window = sum(recorded_dttm > intubation_window_end_dttm, na.rm = TRUE),
      n_window_end_mismatches = sum(intubation_window_end_dttm != intubation_dttm + hours(1), na.rm = TRUE)
    )
  
  post_intubation_event_mapping_qc <- post_intubation_tidal_volume_rows %>%
    group_by(intubation_event_id) %>%
    summarise(
      n_hospitalizations = n_distinct(hospitalization_id, na.rm = FALSE),
      n_anchor_timestamps = n_distinct(intubation_dttm, na.rm = FALSE),
      n_identification_methods = n_distinct(intubation_identification_method, na.rm = FALSE),
      .groups = "drop"
    )
  
  post_intubation_input_checks <- list(
    required_columns_present = length(missing_post_intubation_columns) == 0L,
    input_nonempty = nrow(post_intubation_tidal_volume_rows) > 0L,
    all_cohort_hospitalizations_represented =
      all(cohort_ids %in% post_intubation_tidal_volume_rows$hospitalization_id),
    no_rows_outside_cohort = post_intubation_input_qc$n_rows_outside_cohort == 0L,
    identifiers_and_datetimes_complete = all(
      post_intubation_input_qc$n_missing_hospitalization_id == 0L,
      post_intubation_input_qc$n_missing_event_id == 0L,
      post_intubation_input_qc$n_missing_intubation_dttm == 0L,
      post_intubation_input_qc$n_missing_window_end_dttm == 0L,
      post_intubation_input_qc$n_missing_recorded_dttm == 0L
    ),
    identification_methods_complete =
      post_intubation_input_qc$n_missing_identification_method == 0L,
    identification_methods_valid =
      post_intubation_input_qc$n_invalid_identification_method == 0L,
    event_ids_map_to_one_hospitalization_anchor_and_method = all(
      post_intubation_event_mapping_qc$n_hospitalizations == 1L &
        post_intubation_event_mapping_qc$n_anchor_timestamps == 1L &
        post_intubation_event_mapping_qc$n_identification_methods == 1L
    ),
    all_rows_confirmed_invasive = post_intubation_input_qc$n_noninvasive_rows == 0L,
    sources_valid = post_intubation_input_qc$n_invalid_sources == 0L,
    values_positive = post_intubation_input_qc$n_missing_or_nonpositive_values == 0L,
    respiratory_source_rows_unique = post_intubation_input_qc$n_duplicate_source_rows == 0L,
    rows_within_inclusive_one_hour_window =
      post_intubation_input_qc$n_rows_before_intubation == 0L &&
      post_intubation_input_qc$n_rows_after_window == 0L,
    window_end_matches_one_hour = post_intubation_input_qc$n_window_end_mismatches == 0L,
    source_counts_reconcile =
      post_intubation_input_qc$n_set_rows + post_intubation_input_qc$n_observed_rows ==
      post_intubation_input_qc$n_rows
  )
  
  post_intubation_input_checks_table <- tibble(
    check = names(post_intubation_input_checks),
    passed = as.logical(unlist(post_intubation_input_checks, use.names = FALSE))
  )
  
  post_intubation_method_summary <- post_intubation_tidal_volume_rows %>%
    count(intubation_identification_method, name = "n_rows") %>%
    mutate(percent = round(100 * n_rows / sum(n_rows), 1))
  
  if (!isTRUE(all(post_intubation_input_checks_table$passed))) {
    print(post_intubation_input_qc, width = Inf)
    print(post_intubation_method_summary, n = Inf, width = Inf)
    print(post_intubation_input_checks_table, n = Inf, width = Inf)
  }

  stopifnot(isTRUE(all(post_intubation_input_checks_table$passed)))
  
  required_height_provenance_columns <- c("height_source", "height_gap_days")
  missing_height_provenance_columns <- setdiff(required_height_provenance_columns, names(cohort_heights))
  
  height_provenance_input_qc <- tibble(
    diagnostic = c(
      "Required provenance columns", "Required provenance columns present", "Required provenance columns missing",
      "Height rows with missing source", "Donor-source rows with missing gap", "Non-donor rows with nonzero gap"
    ),
    n = c(
      length(required_height_provenance_columns),
      length(required_height_provenance_columns) - length(missing_height_provenance_columns),
      length(missing_height_provenance_columns),
      if ("height_source" %in% names(cohort_heights)) sum(is.na(cohort_heights$height_source)) else NA_integer_,
      if (all(required_height_provenance_columns %in% names(cohort_heights))) {
        sum(cohort_heights$height_source == "Another encounter within window" & is.na(cohort_heights$height_gap_days), na.rm = TRUE)
      } else {
        NA_integer_
      },
      if (all(required_height_provenance_columns %in% names(cohort_heights))) {
        sum(
          cohort_heights$height_source != "Another encounter within window" &
            (is.na(cohort_heights$height_gap_days) | cohort_heights$height_gap_days != 0),
          na.rm = TRUE
        )
      } else {
        NA_integer_
      }
    )
  )
  
  if (length(missing_height_provenance_columns) > 0L) {
    print(height_provenance_input_qc, n = Inf, width = Inf)
  }
  
  if (length(missing_height_provenance_columns) > 0L) {
    stop(
      "cohort_heights.parquet is missing required provenance column(s): ",
      paste(missing_height_provenance_columns, collapse = ", ")
    )
  }
  
  if (anyNA(cohort_heights$height_source)) {
    print(height_provenance_input_qc, n = Inf, width = Inf)
    stop("cohort_heights.parquet contains missing height_source values.")
  }
  
  if (any(
    cohort_heights$height_source == "Another encounter within window" & is.na(cohort_heights$height_gap_days),
    na.rm = TRUE
  )) {
    print(height_provenance_input_qc, n = Inf, width = Inf)
    stop("At least one donor-derived height is missing height_gap_days.")
  }
  
  if (any(
    cohort_heights$height_source != "Another encounter within window" &
    (is.na(cohort_heights$height_gap_days) | cohort_heights$height_gap_days != 0),
    na.rm = TRUE
  )) {
    print(height_provenance_input_qc, n = Inf, width = Inf)
    stop("At least one non-donor height has a missing or nonzero height_gap_days value.")
  }
  
  
# =============================================================================
# Index-hospitalization ID integrity
# =============================================================================
  
  cohort_id_qc <- tibble(
    diagnostic = c("Total entries in cohort ID vector", "Distinct index hospitalization IDs",
                   "Duplicate entries in cohort ID vector", "Missing entries in cohort ID vector"),
    n = c(length(cohort_ids), n_distinct(cohort_ids, na.rm = TRUE), sum(duplicated(cohort_ids)), sum(is.na(cohort_ids)))
  )
  
  if (anyDuplicated(cohort_ids) || anyNA(cohort_ids)) {
    print(cohort_id_qc, n = Inf)
  }
  
  if (anyDuplicated(cohort_ids)) {
    stop("cohort_ids contains duplicate hospitalization IDs. Resolve this before running Script 02.")
  }
  
  if (anyNA(cohort_ids)) {
    stop("cohort_ids contains missing hospitalization IDs. Resolve this before running Script 02.")
  }
  
# =============================================================================
# Laboratory cleaning
# =============================================================================
  
  message("Script 02 | cleaning clinical variables and building diagnostics")
  labs_threshold_path <- file.path(thresholds, "outlier_thresholds_labs.csv")
  
  if (!file.exists(labs_threshold_path)) {
    stop("The laboratory threshold file was not found: ", labs_threshold_path)
  }
  
  labs_thresholds <- read_csv(labs_threshold_path, show_col_types = FALSE) %>%
    mutate(across(c(lower_limit, upper_limit), as.numeric))
  
  lab_categories <- sort(unique(cohort_labs$lab_category))
  
  labs_thresholds_analysis <- labs_thresholds %>%
    filter(lab_category %in% lab_categories)
  
  relevant_duplicate_thresholds <- labs_thresholds_analysis %>%
    count(lab_category, name = "n_rules") %>%
    filter(n_rules > 1)
  
  if (nrow(relevant_duplicate_thresholds) > 0L) {
    print(relevant_duplicate_thresholds, n = Inf)
    stop("Duplicate threshold rules exist for a lab category present in cohort_labs.")
  }
  
  cohort_labs_clean <- replace_outliers_with_na_long(
    cohort_labs, labs_thresholds_analysis, "lab_category", "lab_value_numeric"
  )
  
  if (nrow(cohort_labs_clean) != nrow(cohort_labs)) {
    stop("Laboratory threshold cleaning unexpectedly changed the number of rows.")
  }
  
  lab_value_diagnostic <- cohort_labs %>%
    transmute(
      hospitalization_id, lab_category, value_before = lab_value_numeric,
      value_after = cohort_labs_clean$lab_value_numeric
    ) %>%
    left_join(
      labs_thresholds_analysis %>% select(lab_category, lower_limit, upper_limit),
      by = "lab_category",
      relationship = "many-to-one"
    ) %>%
    mutate(
      removed_by_cleaning = !is.na(value_before) & is.na(value_after),
      removal_rule = case_when(
        !removed_by_cleaning ~ NA_character_,
        !is.na(lower_limit) & value_before < lower_limit ~ paste0("Below lower threshold (", lower_limit, ")"),
        !is.na(upper_limit) & value_before > upper_limit ~ paste0("Above upper threshold (", upper_limit, ")"),
        TRUE ~ "Removed for unidentified reason"
      )
    )
  
  lab_dropout_by_hospitalization <- crossing(hospitalization_id = cohort_ids, lab_category = lab_categories) %>%
    left_join(
      lab_value_diagnostic %>%
        filter(hospitalization_id %in% cohort_ids) %>%
        group_by(hospitalization_id, lab_category) %>%
        summarise(
          n_usable_before = sum(!is.na(value_before)),
          n_usable_after = sum(!is.na(value_after)),
          removal_rule = if_else(
            n_usable_before > 0 & n_usable_after == 0,
            paste(sort(unique(na.omit(removal_rule))), collapse = "; "),
            NA_character_
          ),
          .groups = "drop"
        ),
      by = c("hospitalization_id", "lab_category"),
      relationship = "one-to-one"
    ) %>%
    mutate(
      across(c(n_usable_before, n_usable_after), ~ replace_na(.x, 0L)),
      had_usable_data_before = n_usable_before > 0,
      has_usable_data_after = n_usable_after > 0,
      lost_all_usable_data = had_usable_data_before & !has_usable_data_after
    )
  
  lab_dropout_summary <- lab_dropout_by_hospitalization %>%
    group_by(lab_category) %>%
    summarise(
      n_index_hospitalizations = n(),
      n_with_usable_data_before = sum(had_usable_data_before),
      n_with_usable_data_after = sum(has_usable_data_after),
      n_lost_all_usable_data = sum(lost_all_usable_data),
      percent_lost_among_those_with_data_before = if_else(
        n_with_usable_data_before > 0,
        round(100 * n_lost_all_usable_data / n_with_usable_data_before, 1),
        NA_real_
      ),
      .groups = "drop"
    )
  
  lab_dropout_rule_summary <- lab_dropout_by_hospitalization %>%
    filter(lost_all_usable_data) %>%
    count(lab_category, removal_rule, name = "n_hospitalizations") %>%
    arrange(lab_category, desc(n_hospitalizations))
  
  lab_integrity_summary <- tibble(
    diagnostic = c("Duplicate entries in cohort ID vector", "Lab rows outside index cohort",
                   "Index hospitalizations with no lab rows", "Exact duplicate lab rows",
                   "Row-count change after cleaning"),
    n = c(sum(duplicated(cohort_ids)), sum(!cohort_labs$hospitalization_id %in% cohort_ids),
          sum(!cohort_ids %in% cohort_labs$hospitalization_id), sum(duplicated(as.data.frame(cohort_labs))),
          nrow(cohort_labs_clean) - nrow(cohort_labs))
  )
  
  
  
  
  lab_summary <- generate_summary_stats(cohort_labs_clean, lab_category, lab_value_numeric)
  
# =============================================================================
# Vital-sign cleaning
# =============================================================================
  
  cohort_vitals_before_deduplication <- cohort_vitals
  cohort_vitals <- distinct(cohort_vitals)
  
  vital_deduplication_summary <- tibble(
    diagnostic = c("Vital rows before exact deduplication", "Vital rows after exact deduplication",
                   "Exact duplicate rows removed"),
    n = c(nrow(cohort_vitals_before_deduplication), nrow(cohort_vitals),
          nrow(cohort_vitals_before_deduplication) - nrow(cohort_vitals))
  )
  
  
  vitals_threshold_path <- file.path(thresholds, "outlier_thresholds_adults_vitals.csv")
  
  if (!file.exists(vitals_threshold_path)) {
    stop("The vital-sign threshold file was not found: ", vitals_threshold_path)
  }
  
  vitals_thresholds <- read_csv(vitals_threshold_path, show_col_types = FALSE) %>%
    mutate(across(c(lower_limit, upper_limit), as.numeric))
  
  vitals_thresholds_analysis <- vitals_thresholds %>%
    filter(vital_category %in% unique(cohort_vitals$vital_category))
  
  duplicate_vital_thresholds <- vitals_thresholds_analysis %>%
    count(vital_category, name = "n_rules") %>%
    filter(n_rules > 1)
  
  if (nrow(duplicate_vital_thresholds) > 0L) {
    print(duplicate_vital_thresholds, n = Inf)
    stop("Duplicate threshold rules exist for a vital category present in cohort_vitals.")
  }
  
  cohort_vitals_qc <- cohort_vitals %>%
    mutate(.qc_row_id = row_number())
  
  cohort_vitals_clean_qc <- replace_outliers_with_na_long(
    cohort_vitals_qc, vitals_thresholds_analysis, "vital_category", "vital_value"
  )
  
  if (nrow(cohort_vitals_clean_qc) != nrow(cohort_vitals_qc) || anyDuplicated(cohort_vitals_clean_qc$.qc_row_id)) {
    stop("Vital-sign threshold cleaning changed row count or duplicated source rows.")
  }
  
  cohort_vitals_clean <- cohort_vitals_clean_qc %>%
    select(-.qc_row_id)
  
  vital_value_diagnostic <- cohort_vitals_qc %>%
    select(.qc_row_id, hospitalization_id, vital_category, value_before = vital_value) %>%
    left_join(
      cohort_vitals_clean_qc %>% select(.qc_row_id, value_after = vital_value),
      by = ".qc_row_id",
      relationship = "one-to-one"
    ) %>%
    left_join(
      vitals_thresholds_analysis %>% select(vital_category, lower_limit, upper_limit),
      by = "vital_category",
      relationship = "many-to-one"
    ) %>%
    mutate(
      removed_by_cleaning = !is.na(value_before) & is.na(value_after),
      removal_rule = case_when(
        removed_by_cleaning & !is.na(lower_limit) & value_before < lower_limit ~
          paste0("Below lower threshold (", lower_limit, ")"),
        removed_by_cleaning & !is.na(upper_limit) & value_before > upper_limit ~
          paste0("Above upper threshold (", upper_limit, ")"),
        removed_by_cleaning ~ "Removed for unidentified reason",
        TRUE ~ NA_character_
      )
    )
  
  vital_categories <- sort(unique(cohort_vitals$vital_category))
  
  vital_dropout_by_hospitalization <- crossing(hospitalization_id = cohort_ids, vital_category = vital_categories) %>%
    left_join(
      vital_value_diagnostic %>%
        filter(hospitalization_id %in% cohort_ids) %>%
        group_by(hospitalization_id, vital_category) %>%
        summarise(
          n_usable_before = sum(!is.na(value_before)),
          n_usable_after = sum(!is.na(value_after)),
          rules_causing_removal = paste(sort(unique(na.omit(removal_rule))), collapse = "; "),
          .groups = "drop"
        ),
      by = c("hospitalization_id", "vital_category"),
      relationship = "one-to-one"
    ) %>%
    mutate(
      across(c(n_usable_before, n_usable_after), ~ replace_na(.x, 0L)),
      had_usable_data_before = n_usable_before > 0,
      has_usable_data_after = n_usable_after > 0,
      lost_all_usable_data = had_usable_data_before & !has_usable_data_after,
      removal_rule = if_else(lost_all_usable_data, rules_causing_removal, NA_character_)
    ) %>%
    select(-rules_causing_removal)
  
  vital_dropout_summary <- vital_dropout_by_hospitalization %>%
    group_by(vital_category) %>%
    summarise(
      n_index_hospitalizations = n(),
      n_with_usable_data_before = sum(had_usable_data_before),
      n_with_usable_data_after = sum(has_usable_data_after),
      n_lost_all_usable_data = sum(lost_all_usable_data),
      percent_lost_among_those_with_data_before = if_else(
        n_with_usable_data_before > 0,
        round(100 * n_lost_all_usable_data / n_with_usable_data_before, 1),
        NA_real_
      ),
      .groups = "drop"
    )
  
  vital_dropout_rule_summary <- vital_dropout_by_hospitalization %>%
    filter(lost_all_usable_data) %>%
    count(vital_category, removal_rule, name = "n_hospitalizations")
  
  vital_integrity_summary <- tibble(
    diagnostic = c("Duplicate entries in cohort ID vector", "Vital rows outside index cohort",
                   "Index hospitalizations with no vital rows", "Exact duplicate vital rows detected",
                   "Exact duplicate vital rows removed", "Row-count change during threshold cleaning"),
    n = c(sum(duplicated(cohort_ids)),
          sum(!cohort_vitals_before_deduplication$hospitalization_id %in% cohort_ids),
          sum(!cohort_ids %in% cohort_vitals_before_deduplication$hospitalization_id),
          sum(duplicated(as.data.frame(cohort_vitals_before_deduplication))),
          nrow(cohort_vitals_before_deduplication) - nrow(cohort_vitals),
          nrow(cohort_vitals_clean) - nrow(cohort_vitals))
  )
  
  vital_duplicate_groups <- cohort_vitals_before_deduplication %>%
    group_by(across(everything())) %>%
    summarise(n_copies = n(), .groups = "drop") %>%
    filter(n_copies > 1)
  
  vital_duplicate_summary <- vital_duplicate_groups %>%
    summarise(
      n_duplicate_groups = n(),
      n_rows_in_duplicate_groups = sum(n_copies),
      n_extra_duplicate_rows = sum(n_copies - 1L),
      n_hospitalizations_affected = n_distinct(hospitalization_id)
    )
  
  vital_duplicate_by_category <- vital_duplicate_groups %>%
    group_by(vital_category) %>%
    summarise(
      n_duplicate_groups = n(),
      n_rows_in_duplicate_groups = sum(n_copies),
      n_extra_duplicate_rows = sum(n_copies - 1L),
      n_hospitalizations_affected = n_distinct(hospitalization_id),
      maximum_number_of_copies = max(n_copies),
      .groups = "drop"
    )
  
  
  
  
  
  
  vital_summary <- generate_summary_stats(cohort_vitals_clean, vital_category, vital_value)
  
# =============================================================================
# Height cleaning
# =============================================================================
  
  cohort_heights_clean <- cohort_heights %>%
    mutate(height_cm = if_else(height_cm < 80 | height_cm > 200, NA_real_, height_cm))
  
  if (nrow(cohort_heights_clean) != nrow(cohort_heights)) {
    stop("Height cleaning unexpectedly changed the number of rows.")
  }
  
  missing_height_provenance_columns_after_cleaning <- setdiff(
    required_height_provenance_columns, names(cohort_heights_clean)
  )
  
  height_provenance_output_qc <- tibble(
    diagnostic = c(
      "Required provenance columns retained", "Required provenance columns missing after cleaning",
      "Height source values changed during cleaning", "Height gap values changed during cleaning"
    ),
    n = c(
      length(required_height_provenance_columns) - length(missing_height_provenance_columns_after_cleaning),
      length(missing_height_provenance_columns_after_cleaning),
      if ("height_source" %in% names(cohort_heights_clean)) {
        sum(cohort_heights_clean$height_source != cohort_heights$height_source, na.rm = TRUE)
      } else {
        NA_integer_
      },
      if ("height_gap_days" %in% names(cohort_heights_clean)) {
        sum(
          !is.na(cohort_heights_clean$height_gap_days) != !is.na(cohort_heights$height_gap_days) |
            (!is.na(cohort_heights_clean$height_gap_days) & !is.na(cohort_heights$height_gap_days) &
               cohort_heights_clean$height_gap_days != cohort_heights$height_gap_days)
        )
      } else {
        NA_integer_
      }
    )
  )
  
  if (
    length(missing_height_provenance_columns_after_cleaning) > 0L ||
    !identical(cohort_heights_clean$height_source, cohort_heights$height_source) ||
    !identical(cohort_heights_clean$height_gap_days, cohort_heights$height_gap_days)
  ) {
    print(height_provenance_output_qc, n = Inf, width = Inf)
  }
  
  if (length(missing_height_provenance_columns_after_cleaning) > 0L) {
    stop(
      "Height cleaning removed required provenance column(s): ",
      paste(missing_height_provenance_columns_after_cleaning, collapse = ", ")
    )
  }
  
  if (!identical(cohort_heights_clean$height_source, cohort_heights$height_source)) {
    stop("Height cleaning unexpectedly changed height_source.")
  }
  
  if (!identical(cohort_heights_clean$height_gap_days, cohort_heights$height_gap_days)) {
    stop("Height cleaning unexpectedly changed height_gap_days.")
  }
  
  height_value_diagnostic <- cohort_heights %>%
    transmute(
      hospitalization_id, height_before = height_cm, height_after = cohort_heights_clean$height_cm,
      removed_by_cleaning = !is.na(height_before) & is.na(height_after),
      removal_rule = case_when(
        removed_by_cleaning & height_before < 80 ~ "Below lower threshold (80 cm)",
        removed_by_cleaning & height_before > 200 ~ "Above upper threshold (200 cm)",
        removed_by_cleaning ~ "Removed for unidentified reason",
        TRUE ~ NA_character_
      )
    )
  
  height_dropout_by_hospitalization <- tibble(hospitalization_id = cohort_ids) %>%
    left_join(
      height_value_diagnostic %>%
        filter(hospitalization_id %in% cohort_ids) %>%
        group_by(hospitalization_id) %>%
        summarise(
          n_height_records = n(),
          n_usable_before = sum(!is.na(height_before)),
          n_usable_after = sum(!is.na(height_after)),
          rules_causing_removal = paste(sort(unique(na.omit(removal_rule))), collapse = "; "),
          .groups = "drop"
        ),
      by = "hospitalization_id",
      relationship = "one-to-one"
    ) %>%
    mutate(
      across(c(n_height_records, n_usable_before, n_usable_after), ~ replace_na(.x, 0L)),
      had_usable_data_before = n_usable_before > 0,
      has_usable_data_after = n_usable_after > 0,
      lost_all_usable_data = had_usable_data_before & !has_usable_data_after,
      removal_rule = case_when(
        lost_all_usable_data & rules_causing_removal != "" ~ rules_causing_removal,
        lost_all_usable_data ~ "Unexpected or unidentified",
        TRUE ~ NA_character_
      )
    ) %>%
    select(-rules_causing_removal)
  
  height_dropout_summary <- height_dropout_by_hospitalization %>%
    summarise(
      n_index_hospitalizations = n(),
      n_with_usable_data_before = sum(had_usable_data_before),
      n_with_usable_data_after = sum(has_usable_data_after),
      n_lost_all_usable_data = sum(lost_all_usable_data),
      n_without_usable_data_before = sum(!had_usable_data_before),
      percent_lost_among_those_with_data_before = if_else(
        n_with_usable_data_before > 0,
        round(100 * n_lost_all_usable_data / n_with_usable_data_before, 1),
        NA_real_
      )
    )
  
  height_dropout_rule_summary <- height_dropout_by_hospitalization %>%
    filter(lost_all_usable_data) %>%
    count(removal_rule, name = "n_hospitalizations") %>%
    arrange(desc(n_hospitalizations))
  
  height_original_rule_summary <- cohort_heights %>%
    summarise(
      n_with_available_height = sum(!is.na(height_cm)),
      n_retained_by_original_rule = sum(!is.na(height_cm) & height_cm >= 90 & height_cm <= 182),
      n_removed_by_original_rule = sum(!is.na(height_cm) & (height_cm < 90 | height_cm > 182)),
      n_removed_below_original_rule = sum(!is.na(height_cm) & height_cm < 90),
      n_removed_above_original_rule = sum(!is.na(height_cm) & height_cm > 182)
    )
  
  height_selected_rule_summary <- cohort_heights %>%
    summarise(
      n_with_available_height = sum(!is.na(height_cm)),
      n_retained_by_selected_rule = sum(!is.na(height_cm) & height_cm >= 80 & height_cm <= 200),
      n_removed_by_selected_rule = sum(!is.na(height_cm) & (height_cm < 80 | height_cm > 200)),
      n_removed_below_selected_rule = sum(!is.na(height_cm) & height_cm < 80),
      n_removed_above_selected_rule = sum(!is.na(height_cm) & height_cm > 200)
    )
  
  height_threshold_decision <- tibble(
    variable = "height_cm", original_rule = "90-182 cm", selected_rule = "80-200 cm",
    rationale = paste0(
      "In the current cohort, the original range would remove ",
      height_original_rule_summary$n_removed_by_original_rule, " of ",
      height_original_rule_summary$n_with_available_height, " hospitalizations with available height. ",
      "The selected broad plausibility range removed ", height_selected_rule_summary$n_removed_by_selected_rule, " of ",
      height_selected_rule_summary$n_with_available_height,
      " and retained plausible heights among shorter children and taller adolescents."
    ),
    literature_status = "Citation needed"
  )
  
  height_threshold_checks <- list(
    original_rule_counts_reconcile =
      height_original_rule_summary$n_retained_by_original_rule +
      height_original_rule_summary$n_removed_by_original_rule ==
      height_original_rule_summary$n_with_available_height,
    selected_rule_counts_reconcile =
      height_selected_rule_summary$n_retained_by_selected_rule +
      height_selected_rule_summary$n_removed_by_selected_rule ==
      height_selected_rule_summary$n_with_available_height,
    selected_rule_matches_cleaned_output =
      height_selected_rule_summary$n_retained_by_selected_rule ==
      sum(!is.na(cohort_heights_clean$height_cm)),
    selected_rule_removals_match_dropout =
      height_selected_rule_summary$n_removed_by_selected_rule ==
      height_dropout_summary$n_lost_all_usable_data
  )
  
  height_threshold_checks_table <- tibble(
    check = names(height_threshold_checks),
    passed = as.logical(unlist(height_threshold_checks, use.names = FALSE))
  )
  
  if (!isTRUE(all(height_threshold_checks_table$passed))) {
    print(height_original_rule_summary, width = Inf)
    print(height_selected_rule_summary, width = Inf)
    print(height_threshold_checks_table, n = Inf, width = Inf)
  }

  stopifnot(isTRUE(all(height_threshold_checks_table$passed)))
  
  height_integrity_summary <- tibble(
    diagnostic = c("Duplicate entries in cohort ID vector", "Height rows outside index cohort",
                   "Index hospitalizations with no height row", "Hospitalization IDs with multiple height rows",
                   "Exact duplicate height rows", "Row-count change after cleaning"),
    n = c(sum(duplicated(cohort_ids)), sum(!cohort_heights$hospitalization_id %in% cohort_ids),
          sum(!cohort_ids %in% cohort_heights$hospitalization_id),
          cohort_heights %>% count(hospitalization_id) %>% filter(n > 1) %>% nrow(),
          sum(duplicated(as.data.frame(cohort_heights))), nrow(cohort_heights_clean) - nrow(cohort_heights))
  )
  
  
  
  
  
# =============================================================================
# Respiratory-support cleaning
# =============================================================================
  
  respiratory_threshold_path <- file.path(thresholds, "outlier_thresholds_respiratory_support.csv")
  
  if (!file.exists(respiratory_threshold_path)) {
    stop("The respiratory-support threshold file was not found: ", respiratory_threshold_path)
  }
  
  resp_thresholds <- read_csv(respiratory_threshold_path, show_col_types = FALSE) %>%
    mutate(across(c(lower_limit, upper_limit), as.numeric))
  
  resp_thresholds_analysis <- resp_thresholds %>%
    filter(variable_name %in% names(resp_waterfall)) %>%
    mutate(
      lower_limit = case_when(
        variable_name %in% c("flow_rate_set", "peep_set", "pressure_control_set", "pressure_support_set") ~ 0,
        TRUE ~ lower_limit
      ),
      upper_limit = if_else(variable_name %in% c("resp_rate_set", "resp_rate_obs"), 80, upper_limit)
    )
  
  duplicate_resp_thresholds <- resp_thresholds_analysis %>%
    count(variable_name, name = "n_rules") %>%
    filter(n_rules > 1)
  
  if (nrow(duplicate_resp_thresholds) > 0L) {
    print(duplicate_resp_thresholds, n = Inf)
    stop("Duplicate threshold rules exist for respiratory variables present in resp_waterfall.")
  }
  
  resp_waterfall_before_cleaning <- resp_waterfall
  resp_waterfall_clean <- resp_waterfall
  
  for (i in seq_len(nrow(resp_thresholds_analysis))) {
    var_name <- resp_thresholds_analysis$variable_name[[i]]
    lower <- resp_thresholds_analysis$lower_limit[[i]]
    upper <- resp_thresholds_analysis$upper_limit[[i]]
    
    if (var_name %in% names(resp_waterfall_clean) && !is.na(lower) && !is.na(upper)) {
      n_before <- sum(!is.na(resp_waterfall_clean[[var_name]]))
      
      resp_waterfall_clean[[var_name]] <- if_else(
        resp_waterfall_clean[[var_name]] < lower | resp_waterfall_clean[[var_name]] > upper,
        NA_real_,
        as.numeric(resp_waterfall_clean[[var_name]])
      )
      
      n_after <- sum(!is.na(resp_waterfall_clean[[var_name]]))
    }
  }
  
  if (nrow(resp_waterfall_clean) != nrow(resp_waterfall_before_cleaning)) {
    stop("Respiratory-support cleaning unexpectedly changed the number of rows.")
  }
  
  resp_variables <- resp_thresholds_analysis %>%
    filter(!is.na(lower_limit), !is.na(upper_limit)) %>%
    pull(variable_name) %>%
    unique()
  
  resp_value_diagnostic <- map_dfr(resp_variables, function(var_name) {
    threshold_row <- resp_thresholds_analysis %>%
      filter(variable_name == var_name)
    
    lower <- threshold_row$lower_limit[[1]]
    upper <- threshold_row$upper_limit[[1]]
    
    tibble(
      hospitalization_id = resp_waterfall_before_cleaning$hospitalization_id,
      variable_name = var_name,
      value_before = resp_waterfall_before_cleaning[[var_name]],
      value_after = resp_waterfall_clean[[var_name]],
      lower_limit = lower,
      upper_limit = upper
    ) %>%
      mutate(
        removed_by_cleaning = !is.na(value_before) & is.na(value_after),
        removed_by_lower_threshold = removed_by_cleaning & value_before < lower_limit,
        removed_by_upper_threshold = removed_by_cleaning & value_before > upper_limit
      )
  })
  
  resp_dropout_by_hospitalization <- crossing(hospitalization_id = cohort_ids, variable_name = resp_variables) %>%
    left_join(
      resp_value_diagnostic %>%
        filter(hospitalization_id %in% cohort_ids) %>%
        group_by(hospitalization_id, variable_name) %>%
        summarise(
          n_usable_before = sum(!is.na(value_before)),
          n_usable_after = sum(!is.na(value_after)),
          n_removed_lower = sum(removed_by_lower_threshold, na.rm = TRUE),
          n_removed_upper = sum(removed_by_upper_threshold, na.rm = TRUE),
          lower_limit = first(lower_limit),
          upper_limit = first(upper_limit),
          .groups = "drop"
        ),
      by = c("hospitalization_id", "variable_name"),
      relationship = "one-to-one"
    ) %>%
    mutate(
      across(c(n_usable_before, n_usable_after, n_removed_lower, n_removed_upper), ~ replace_na(.x, 0L)),
      had_usable_data_before = n_usable_before > 0,
      has_usable_data_after = n_usable_after > 0,
      lost_all_usable_data = had_usable_data_before & !has_usable_data_after,
      threshold_rule_responsible = case_when(
        !lost_all_usable_data ~ NA_character_,
        n_removed_lower > 0 & n_removed_upper > 0 ~
          paste0("Both lower (", lower_limit, ") and upper (", upper_limit, ") thresholds"),
        n_removed_lower > 0 ~ paste0("Lower threshold (", lower_limit, ")"),
        n_removed_upper > 0 ~ paste0("Upper threshold (", upper_limit, ")"),
        TRUE ~ "Unexpected or unidentified"
      )
    )
  
  resp_dropout_summary <- resp_dropout_by_hospitalization %>%
    group_by(variable_name) %>%
    summarise(
      n_index_hospitalizations = n(),
      n_with_usable_data_before = sum(had_usable_data_before),
      n_with_usable_data_after = sum(has_usable_data_after),
      n_lost_all_usable_data = sum(lost_all_usable_data),
      n_without_usable_data_before = sum(!had_usable_data_before),
      percent_lost_among_those_with_data_before = if_else(
        n_with_usable_data_before > 0,
        round(100 * n_lost_all_usable_data / n_with_usable_data_before, 1),
        NA_real_
      ),
      .groups = "drop"
    )
  
  resp_dropout_rule_summary <- resp_dropout_by_hospitalization %>%
    filter(lost_all_usable_data) %>%
    count(variable_name, threshold_rule_responsible, name = "n_hospitalizations") %>%
    arrange(variable_name, desc(n_hospitalizations))
  
  resp_integrity_summary <- tibble(
    diagnostic = c("Duplicate entries in cohort ID vector", "Waterfall rows outside index cohort",
                   "Index hospitalizations with no waterfall rows", "Exact duplicate waterfall rows",
                   "Row-count change after cleaning", "Threshold variables evaluated"),
    n = c(sum(duplicated(cohort_ids)), sum(!resp_waterfall$hospitalization_id %in% cohort_ids),
          sum(!cohort_ids %in% resp_waterfall$hospitalization_id), sum(duplicated(as.data.frame(resp_waterfall))),
          nrow(resp_waterfall_clean) - nrow(resp_waterfall), length(resp_variables))
  )
  
  
  
  
# =============================================================================
# Values responsible for complete respiratory-variable loss
# =============================================================================
  
  resp_complete_loss_keys <- resp_dropout_by_hospitalization %>%
    filter(lost_all_usable_data) %>%
    select(hospitalization_id, variable_name, threshold_rule_responsible)
  
  resp_complete_loss_detail <- map_dfr(resp_variables, function(var_name) {
    tibble(
      hospitalization_id = resp_waterfall_before_cleaning$hospitalization_id,
      device_category = resp_waterfall_before_cleaning$device_category,
      is_scaffold = resp_waterfall_before_cleaning$is_scaffold,
      variable_name = var_name,
      value_before = resp_waterfall_before_cleaning[[var_name]],
      value_after = resp_waterfall_clean[[var_name]]
    )
  }) %>%
    inner_join(
      resp_complete_loss_keys,
      by = c("hospitalization_id", "variable_name"),
      relationship = "many-to-one"
    ) %>%
    filter(!is.na(value_before), is.na(value_after))
  
  resp_complete_loss_value_summary <- resp_complete_loss_detail %>%
    group_by(variable_name, threshold_rule_responsible) %>%
    summarise(
      n_hospitalizations = n_distinct(hospitalization_id),
      n_removed_measurements = n(),
      minimum_removed_value = min(value_before, na.rm = TRUE),
      q1_removed_value = as.numeric(quantile(value_before, 0.25, na.rm = TRUE)),
      median_removed_value = median(value_before, na.rm = TRUE),
      q3_removed_value = as.numeric(quantile(value_before, 0.75, na.rm = TRUE)),
      maximum_removed_value = max(value_before, na.rm = TRUE),
      .groups = "drop"
    )
  
  
# =============================================================================
# Targeted set tidal-volume unit review
# =============================================================================
  
  if (anyDuplicated(cohort_weights$hospitalization_id)) {
    stop("cohort_weights contains multiple rows for at least one hospitalization ID.")
  }
  
  tidal_volume_complete_loss_ids <- resp_complete_loss_keys %>%
    filter(variable_name == "tidal_volume_set") %>%
    pull(hospitalization_id)
  
  has_tidal_volume_inference_flag <- "tidal_volume_set_inferred" %in% names(resp_waterfall_before_cleaning)
  
  tidal_volume_complete_loss_detail <- resp_waterfall_before_cleaning %>%
    filter(hospitalization_id %in% tidal_volume_complete_loss_ids, !is.na(tidal_volume_set)) %>%
    select(
      hospitalization_id,
      any_of(c("recorded_dttm", "device_category", "support", "mode_name", "mode_category", "is_scaffold",
               "tidal_volume_set_inferred")),
      tidal_volume_set,
      tidal_volume_obs
    )
  
  if (!has_tidal_volume_inference_flag) {
    tidal_volume_complete_loss_detail <- tidal_volume_complete_loss_detail %>%
      mutate(tidal_volume_set_inferred = NA)
  }
  
  tidal_volume_complete_loss_detail <- tidal_volume_complete_loss_detail %>%
    left_join(
      cohort_weights %>% select(hospitalization_id, weight_kg),
      by = "hospitalization_id",
      relationship = "many-to-one"
    ) %>%
    mutate(
      row_type = case_when(
        is.na(is_scaffold) ~ "Missing scaffold indicator",
        is_scaffold ~ "Scaffold row",
        TRUE ~ "Real row"
      ),
      observed_to_weight_scaled_set_ratio = if_else(
        !is.na(tidal_volume_obs) & !is.na(weight_kg) & weight_kg > 0 & tidal_volume_set > 0,
        tidal_volume_obs / (tidal_volume_set * weight_kg),
        NA_real_
      )
    )
  
  tidal_volume_complete_loss_summary <- tidal_volume_complete_loss_detail %>%
    filter(row_type == "Real row") %>%
    group_by(hospitalization_id) %>%
    summarise(
      weight_kg = first(weight_kg),
      n_real_rows = n(),
      minimum_set_value = min(tidal_volume_set, na.rm = TRUE),
      median_set_value = median(tidal_volume_set, na.rm = TRUE),
      maximum_set_value = max(tidal_volume_set, na.rm = TRUE),
      n_rows_with_observed_tv = sum(!is.na(tidal_volume_obs)),
      median_observed_tv = if_else(
        any(!is.na(tidal_volume_obs)),
        median(tidal_volume_obs, na.rm = TRUE),
        NA_real_
      ),
      median_observed_to_weight_scaled_set_ratio = if_else(
        any(!is.na(observed_to_weight_scaled_set_ratio)),
        median(observed_to_weight_scaled_set_ratio, na.rm = TRUE),
        NA_real_
      ),
      n_inferred_rows = sum(tidal_volume_set_inferred %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(median_set_value)
  
  tidal_volume_unit_pattern_summary <- tidal_volume_complete_loss_summary %>%
    summarise(
      n_hospitalizations = n(),
      n_with_observed_tv = sum(!is.na(median_observed_tv)),
      n_consistent_with_ml_per_kg_mapping = sum(
        !is.na(median_observed_to_weight_scaled_set_ratio) &
          median_observed_to_weight_scaled_set_ratio >= 0.5 &
          median_observed_to_weight_scaled_set_ratio <= 2
      ),
      n_not_consistent_with_ml_per_kg_mapping = sum(
        !is.na(median_observed_to_weight_scaled_set_ratio) &
          (median_observed_to_weight_scaled_set_ratio < 0.5 |
             median_observed_to_weight_scaled_set_ratio > 2)
      ),
      n_without_enough_information = sum(is.na(median_observed_to_weight_scaled_set_ratio))
    )
  
  
# =============================================================================
# Respiratory-threshold decisions
# =============================================================================
  
  respiratory_threshold_decisions <- tibble(
    variable_name = c("flow_rate_set", "peep_set", "pressure_control_set", "pressure_support_set", "resp_rate_set",
                      "resp_rate_obs", "tidal_volume_set"),
    original_bound = c("Lower = 10", "Lower = 5", "Lower = 5", "Lower = 5", "Upper = 35", "Upper = 35",
                       "Lower = 25"),
    selected_bound = c("Lower = 0", "Lower = 0", "Lower = 0", "Lower = 0", "Upper = 80", "Upper = 80",
                       "Lower = 25"),
    decision = c("Changed", "Changed", "Changed", "Changed", "Changed", "Changed", "Retained"),
    rationale = c(
      "Values below 10 were common for conventional nasal cannula; the original threshold caused substantial loss.",
      "PEEP of 4 is clinically possible; the original threshold caused complete loss for one hospitalization.",
      "Nonnegative settings of 0-4 were retained; negative values remain excluded.",
      "Zero may indicate pressure support turned off rather than an erroneous measurement.",
      "Rates above 35 can occur in pediatric patients; values up to 80 are retained provisionally.",
      "Rates above 35 can occur in pediatric patients; values up to 80 are retained provisionally.",
      paste0(
        "Values below 25 show mixed units and frequently resemble mL/kg rather than absolute mL; ",
        "lowering the threshold could combine incompatible measurements."
      )
    ),
    literature_status = c("Citation needed", "Citation needed", "Citation needed", "Citation needed",
                          "Citation needed", "Citation needed", "ETL/data-dictionary confirmation needed")
  )
  
  
# =============================================================================
# Combined Script 02 summaries
# =============================================================================
  
  script02_dropout_summary <- bind_rows(
    lab_dropout_summary %>%
      transmute(
        data_domain = "Laboratory", data_element = lab_category, n_index_hospitalizations,
        n_with_usable_data_before, n_with_usable_data_after, n_lost_all_usable_data,
        percent_lost_among_those_with_data_before
      ),
    vital_dropout_summary %>%
      transmute(
        data_domain = "Vital sign", data_element = vital_category, n_index_hospitalizations,
        n_with_usable_data_before, n_with_usable_data_after, n_lost_all_usable_data,
        percent_lost_among_those_with_data_before
      ),
    height_dropout_summary %>%
      transmute(
        data_domain = "Anthropometric", data_element = "height_cm", n_index_hospitalizations,
        n_with_usable_data_before, n_with_usable_data_after, n_lost_all_usable_data,
        percent_lost_among_those_with_data_before
      ),
    resp_dropout_summary %>%
      transmute(
        data_domain = "Respiratory support", data_element = variable_name, n_index_hospitalizations,
        n_with_usable_data_before, n_with_usable_data_after, n_lost_all_usable_data,
        percent_lost_among_those_with_data_before
      )
  ) %>%
    arrange(data_domain, data_element)
  
  script02_complete_loss_rules <- bind_rows(
    lab_dropout_rule_summary %>%
      transmute(
        data_domain = "Laboratory", data_element = lab_category, threshold_rule = removal_rule,
        n_hospitalizations
      ),
    vital_dropout_rule_summary %>%
      transmute(
        data_domain = "Vital sign", data_element = vital_category, threshold_rule = removal_rule,
        n_hospitalizations
      ),
    height_dropout_rule_summary %>%
      transmute(
        data_domain = "Anthropometric", data_element = "height_cm", threshold_rule = removal_rule,
        n_hospitalizations
      ),
    resp_dropout_rule_summary %>%
      transmute(
        data_domain = "Respiratory support", data_element = variable_name,
        threshold_rule = threshold_rule_responsible, n_hospitalizations
      )
  ) %>%
    arrange(data_domain, data_element, threshold_rule)
  
  
  
# =============================================================================
# Build the complete Script 02 diagnostic section required by the QMD
# =============================================================================
  
  script02_diagnostics <- list(
    cohort_id_qc = cohort_id_qc,
    post_intubation_tidal_volume = list(
      input_summary = post_intubation_input_qc,
      method_summary = post_intubation_method_summary,
      event_mapping_qc = post_intubation_event_mapping_qc,
      checks = post_intubation_input_checks_table
    ),
    combined_dropout_summary = script02_dropout_summary,
    combined_complete_loss_rules = script02_complete_loss_rules,
    labs = list(
      dropout_summary = lab_dropout_summary,
      dropout_rule_summary = lab_dropout_rule_summary,
      integrity_summary = lab_integrity_summary
    ),
    vitals = list(
      dropout_summary = vital_dropout_summary,
      dropout_rule_summary = vital_dropout_rule_summary,
      integrity_summary = vital_integrity_summary,
      deduplication_summary = vital_deduplication_summary,
      duplicate_summary = vital_duplicate_summary,
      duplicate_by_category = vital_duplicate_by_category
    ),
    heights = list(
      dropout_summary = height_dropout_summary,
      dropout_rule_summary = height_dropout_rule_summary,
      integrity_summary = height_integrity_summary,
      original_rule_summary = height_original_rule_summary,
      selected_rule_summary = height_selected_rule_summary,
      threshold_decision = height_threshold_decision,
      threshold_checks = height_threshold_checks_table,
      provenance = list(
        required_columns = required_height_provenance_columns,
        missing_input_columns = missing_height_provenance_columns,
        missing_output_columns = missing_height_provenance_columns_after_cleaning,
        input_qc = height_provenance_input_qc,
        output_qc = height_provenance_output_qc,
        passed = length(missing_height_provenance_columns) == 0L &&
          length(missing_height_provenance_columns_after_cleaning) == 0L &&
          identical(cohort_heights_clean$height_source, cohort_heights$height_source) &&
          identical(cohort_heights_clean$height_gap_days, cohort_heights$height_gap_days)
      )
    ),
    respiratory_support = list(
      dropout_summary = resp_dropout_summary,
      dropout_rule_summary = resp_dropout_rule_summary,
      integrity_summary = resp_integrity_summary,
      complete_loss_value_summary = resp_complete_loss_value_summary,
      threshold_decisions = respiratory_threshold_decisions,
      tidal_volume_complete_loss = list(
        summary = tidal_volume_unit_pattern_summary,
        hospitalization_level = tidal_volume_complete_loss_summary
      )
    )
  )
  
# =============================================================================
# Save cleaned outputs
# =============================================================================
  
  message("Script 02 | saving cleaned outputs and validating diagnostics")
  summary_stats_dir <- file.path(output_dir, "summary_stats")
  dir.create(summary_stats_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_parquet(cohort_labs_clean, file.path(output_dir, "cohort_labs_clean.parquet"))
  write_parquet(cohort_vitals_clean, file.path(output_dir, "cohort_vitals_clean.parquet"))
  write_parquet(cohort_heights_clean, file.path(output_dir, "cohort_heights_clean.parquet"))
  write_parquet(resp_waterfall_clean, file.path(output_dir, "resp_support_waterfall_clean.parquet"))
  write_csv(lab_summary, file.path(summary_stats_dir, paste0("lab_summary_", site_name, ".csv")))
  write_csv(vital_summary, file.path(summary_stats_dir, paste0("vital_summary_", site_name, ".csv")))
  
  required_output_files <- c("cohort_labs_clean.parquet", "cohort_vitals_clean.parquet",
                             "cohort_heights_clean.parquet", "resp_support_waterfall_clean.parquet")
  
  missing_output_files <- required_output_files[!file.exists(file.path(output_dir, required_output_files))]
  
  if (length(missing_output_files) > 0L) {
    stop("Script 02 failed to create required output(s): ", paste(missing_output_files, collapse = ", "))
  }
  
# =============================================================================
# Append Script 02 diagnostics to the existing combined RDS
# =============================================================================
  
  cohort_diagnostics <- readRDS(diagnostic_rds_path)
  
  if (!is.list(cohort_diagnostics)) {
    stop("The Script 01 diagnostic RDS does not contain a list.")
  }
  
  cohort_diagnostics$quality_checks_script02 <- script02_diagnostics
  saveRDS(cohort_diagnostics, diagnostic_rds_path)
  
# =============================================================================
# Validate the saved RDS in a fresh object
# =============================================================================
  
  diagnostic_readback <- readRDS(diagnostic_rds_path)
  
  if (!identical(diagnostic_readback, cohort_diagnostics)) {
    stop("The combined diagnostic RDS failed its complete Script 02 read-back comparison.")
  }
  
  required_script02_components <- c(
    "cohort_id_qc", "post_intubation_tidal_volume", "combined_dropout_summary", "combined_complete_loss_rules",
    "labs", "vitals", "heights", "respiratory_support"
  )
  
  required_post_intubation_components <- c("input_summary", "method_summary", "event_mapping_qc", "checks")
  required_lab_components <- c("dropout_summary", "dropout_rule_summary", "integrity_summary")
  
  required_vital_components <- c(
    "dropout_summary", "dropout_rule_summary", "integrity_summary", "deduplication_summary", "duplicate_summary",
    "duplicate_by_category"
  )
  
  required_height_components <- c(
    "dropout_summary", "dropout_rule_summary", "integrity_summary", "original_rule_summary",
    "selected_rule_summary", "threshold_decision", "threshold_checks", "provenance"
  )
  
  required_height_provenance_components <- c(
    "required_columns", "missing_input_columns", "missing_output_columns", "input_qc", "output_qc", "passed"
  )
  
  required_resp_components <- c(
    "dropout_summary", "dropout_rule_summary", "integrity_summary", "complete_loss_value_summary",
    "threshold_decisions", "tidal_volume_complete_loss"
  )
  
  required_tidal_components <- c("summary", "hospitalization_level")
  
  find_missing_components <- function(required_components, diagnostic_object, prefix = "") {
    missing_components <- setdiff(required_components, names(diagnostic_object))
    
    if (length(missing_components) == 0L) {
      return(character(0))
    }
    
    paste0(prefix, missing_components)
  }
  
  missing_report_components <- c(
    find_missing_components(required_script02_components, diagnostic_readback$quality_checks_script02),
    find_missing_components(
      required_post_intubation_components,
      diagnostic_readback$quality_checks_script02$post_intubation_tidal_volume,
      prefix = "post_intubation_tidal_volume$"
    ),
    find_missing_components(required_lab_components, diagnostic_readback$quality_checks_script02$labs, prefix = "labs$"),
    find_missing_components(
      required_vital_components, diagnostic_readback$quality_checks_script02$vitals, prefix = "vitals$"
    ),
    find_missing_components(
      required_height_components, diagnostic_readback$quality_checks_script02$heights, prefix = "heights$"
    ),
    find_missing_components(
      required_height_provenance_components,
      diagnostic_readback$quality_checks_script02$heights$provenance,
      prefix = "heights$provenance$"
    ),
    find_missing_components(
      required_resp_components, diagnostic_readback$quality_checks_script02$respiratory_support,
      prefix = "respiratory_support$"
    ),
    find_missing_components(
      required_tidal_components,
      diagnostic_readback$quality_checks_script02$respiratory_support$tidal_volume_complete_loss,
      prefix = "respiratory_support$tidal_volume_complete_loss$"
    )
  )
  
  if (length(missing_report_components) > 0L) {
    stop("The saved Script 02 diagnostic section is missing: ", paste(missing_report_components, collapse = ", "))
  }
  
  if (!identical(
    diagnostic_readback$quality_checks_script02$combined_dropout_summary,
    script02_diagnostics$combined_dropout_summary
  )) {
    stop("The Script 02 diagnostic RDS failed its read-back comparison.")
  }
  
  if (!identical(
    diagnostic_readback$quality_checks_script02$heights$threshold_checks,
    script02_diagnostics$heights$threshold_checks
  )) {
    stop("The Script 02 height-threshold diagnostics failed their read-back comparison.")
  }
  
  if (!isTRUE(all(diagnostic_readback$quality_checks_script02$heights$threshold_checks$passed))) {
    stop("At least one saved Script 02 height-threshold integrity check failed.")
  }
  
  if (!isTRUE(diagnostic_readback$quality_checks_script02$heights$provenance$passed)) {
    stop("The saved Script 02 height-provenance integrity check did not pass.")
  }
  
  if (!identical(
    diagnostic_readback$quality_checks_script02$heights$provenance,
    script02_diagnostics$heights$provenance
  )) {
    stop("The Script 02 height-provenance diagnostics failed their read-back comparison.")
  }
  
  if (!identical(
    diagnostic_readback$quality_checks_script02$post_intubation_tidal_volume,
    script02_diagnostics$post_intubation_tidal_volume
  )) {
    stop("The Script 02 post-intubation tidal-volume diagnostics failed their read-back comparison.")
  }
  
  if (!isTRUE(all(diagnostic_readback$quality_checks_script02$post_intubation_tidal_volume$checks$passed))) {
    stop("At least one saved Script 02 post-intubation tidal-volume integrity check failed.")
  }
  
  message("Script 02 complete.")
