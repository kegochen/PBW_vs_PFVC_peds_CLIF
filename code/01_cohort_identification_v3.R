# =============================================================================
# Script 01: Cohort Identification
# =============================================================================

rm(list = ls())

library(tidyverse)
library(arrow)
library(lubridate)
library(data.table)

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

attrition_utility_path <- file.path(config$utils_dir, "attrition_log.R")

if (!file.exists(attrition_utility_path)) {
  stop("Attrition utility not found: ", attrition_utility_path)
}

source(attrition_utility_path)

# Source-table standardization will be completed separately; paths are configuration-driven.

site_name <- config$site_name
tables_path <- config$tables_path
file_type <- config$file_type
output_dir <- config$output_dir
diagnostic_dir <- config$cohort_diagnostic_dir

required_script01_directories <- c(output_dir, diagnostic_dir)

for (directory in required_script01_directories) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

missing_script01_directories <- required_script01_directories[
  !dir.exists(required_script01_directories)
]

if (length(missing_script01_directories) > 0L) {
  stop(
    "Script 01 could not create required output directory(s): ",
    paste(missing_script01_directories, collapse = "; ")
  )
}

message("Script 01 starting | site: ", site_name, " | input format: ", file_type)
message("Intermediate output directory: ", output_dir)
message("Cohort diagnostic directory: ", diagnostic_dir)

# =============================================================================
# Supplemental height source
# =============================================================================
# Supplemental height is activated only by config$supplemental_height.
# A null or omitted setting uses clif_vitals alone.
# Documentation fields in the JSON configuration never activate a source.
#
# Supported source formats:
#   michigan_crosswalk: Michigan encounter crosswalk stored as Parquet.
#   height_csv: Standard measurement-level CSV.
#
# Preserve optional-source behavior: an unavailable configured file does not
# prevent use of clif_vitals. The unavailable source is reported explicitly.
#
# Measurement dates are provenance fields, not replacements for hospitalization
# admission timestamps in the existing cross-encounter donor-window calculation.

supplemental_height_setting <- config[["supplemental_height"]]
supplemental_height_configured <- !is.null(supplemental_height_setting)

supplemental_height_format <- if (supplemental_height_configured) {
  supplemental_height_setting$format
} else {
  NA_character_
}

supplemental_height_path <- if (supplemental_height_configured) {
  supplemental_height_setting$path
} else {
  NA_character_
}

supplemental_height_available <- supplemental_height_configured &&
  file.exists(supplemental_height_path) &&
  !dir.exists(supplemental_height_path)

message("Supplemental height source configured: ", supplemental_height_configured)
message("Supplemental height source available: ", supplemental_height_available)

if (supplemental_height_configured) {
  message("Supplemental height format: ", supplemental_height_format)
  message("Supplemental height source: ", supplemental_height_path)

  if (!supplemental_height_available) {
    message("Configured supplemental height file unavailable; continuing with clif_vitals alone.")
  }
} else {
  message("Supplemental height disabled; using clif_vitals alone.")
}

# =============================================================================
# Load CLIF tables
# =============================================================================

table_names <- c("patient", "hospitalization", "adt", "resp_support", "vitals", "labs", "medication",
                 "patient_assessment", "patient_procedures")

supported_file_types <- c("parquet", "csv", "fst")

if (!file_type %in% supported_file_types) {
  stop("Unsupported file type: ", file_type, ". Supported types: ", paste(supported_file_types, collapse = ", "), ".")
}

table_paths <- set_names(file.path(tables_path, paste0("clif_", table_names, ".", file_type)), table_names)
missing_table_paths <- table_paths[!file.exists(table_paths)]

if (length(missing_table_paths) > 0L) {
  stop("Missing required CLIF table(s): ", paste(missing_table_paths, collapse = "; "))
}

read_clif_table <- function(path, type) {
  switch(
    type,
    parquet = arrow::read_parquet(path),
    csv = readr::read_csv(path, show_col_types = FALSE),
    fst = {
      if (!requireNamespace("fst", quietly = TRUE)) {
        stop("The fst package is required when file_type is 'fst'.")
      }

      fst::read_fst(path)
    }
  )
}

clif_tables <- map(table_paths, read_clif_table, type = file_type)
walk2(names(clif_tables), clif_tables, ~ message("Loaded ", .x, ": ", nrow(.y), " rows"))

clif_patient <- clif_tables$patient
clif_hospitalization <- clif_tables$hospitalization
clif_adt <- clif_tables$adt

clif_respiratory_support <- clif_tables$resp_support

clif_respiratory_support <- clif_respiratory_support %>%
  mutate(.resp_source_row_id = row_number())

clif_vitals <- clif_tables$vitals
clif_labs <- clif_tables$labs

clif_meds <- clif_tables$medication

clif_assessments <- clif_tables$patient_assessment

clif_patient_procedures <- clif_tables$patient_procedures

rm(clif_tables, table_paths, missing_table_paths, supported_file_types)

# =============================================================================
# Normalize and reconcile case-variant identifiers
# =============================================================================

normalize_clif_id <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  na_if(x, "")
}

# -----------------------------------------------------------------------------
# Review hospitalization IDs before normalization
# -----------------------------------------------------------------------------

hospitalization_id_case_qc <- clif_hospitalization %>%
  transmute(
    hospitalization_id_raw = str_trim(as.character(hospitalization_id)),
    hospitalization_id_normalized = normalize_clif_id(hospitalization_id),
    patient_id_normalized = normalize_clif_id(patient_id),
    admission_dttm, discharge_dttm, age_at_admission
  ) %>%
  group_by(hospitalization_id_normalized) %>%
  summarise(
    n_rows = n(),
    n_raw_spellings = n_distinct(hospitalization_id_raw, na.rm = TRUE),
    n_patients = n_distinct(patient_id_normalized, na.rm = FALSE),
    n_admission_datetimes = n_distinct(admission_dttm, na.rm = FALSE),
    n_discharge_datetimes = n_distinct(discharge_dttm, na.rm = FALSE),
    n_ages = n_distinct(age_at_admission, na.rm = FALSE),
    .groups = "drop"
  ) %>%
  summarise(
    n_source_rows = sum(n_rows),
    n_normalized_hospitalization_ids = n(),
    n_ids_with_case_variants = sum(n_raw_spellings > 1),
    n_ids_with_core_identity_conflict = sum(
      n_patients > 1 | n_admission_datetimes > 1 | n_discharge_datetimes > 1 | n_ages > 1
    )
  )

message("---- Hospitalization-ID capitalization QC ----")
print(hospitalization_id_case_qc, width = Inf)

if (hospitalization_id_case_qc$n_ids_with_core_identity_conflict > 0) {
  stop("At least one normalized hospitalization ID has conflicting patient, admission, discharge, or age information.")
}

# -----------------------------------------------------------------------------
# Normalize identifiers in every source table
# -----------------------------------------------------------------------------

clif_hospitalization <- clif_hospitalization %>%
  mutate(.source_row_order = row_number())

tables_to_normalize <- c("clif_patient", "clif_hospitalization", "clif_adt", "clif_respiratory_support",
                         "clif_vitals", "clif_labs", "clif_meds", "clif_assessments", "clif_patient_procedures")

for (object_name in tables_to_normalize) {
  current_table <- get(object_name, inherits = FALSE) %>%
    mutate(across(any_of(c("patient_id", "hospitalization_id")), normalize_clif_id))

  assign(object_name, current_table, envir = .GlobalEnv)
}

rm(current_table, object_name, tables_to_normalize)

# -----------------------------------------------------------------------------
# Verify patient-procedure schema and hospitalization linkage
# -----------------------------------------------------------------------------

required_patient_procedure_columns <- c("hospitalization_id")
missing_patient_procedure_columns <- setdiff(required_patient_procedure_columns, names(clif_patient_procedures))

if (length(missing_patient_procedure_columns) > 0L) {
  stop(
    "clif_patient_procedures is missing required column(s): ",
    paste(missing_patient_procedure_columns, collapse = ", ")
  )
}

patient_procedure_load_qc <- tibble(
  diagnostic = c(
    "Patient-procedure rows loaded", "Rows with missing hospitalization ID",
    "Rows linked to the hospitalization table", "Rows not linked to the hospitalization table"
  ),
  n = c(
    nrow(clif_patient_procedures), sum(is.na(clif_patient_procedures$hospitalization_id)),
    sum(clif_patient_procedures$hospitalization_id %in% clif_hospitalization$hospitalization_id, na.rm = TRUE),
    sum(
      !is.na(clif_patient_procedures$hospitalization_id) &
        !clif_patient_procedures$hospitalization_id %in% clif_hospitalization$hospitalization_id
    )
  )
)

message("---- Patient-procedure table load and linkage ----")
print(patient_procedure_load_qc, n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# Standardize equivalent discharge labels
# -----------------------------------------------------------------------------

clif_hospitalization <- clif_hospitalization %>%
  mutate(
    discharge_category = case_when(
      str_to_lower(str_trim(as.character(discharge_category))) %in% c("dead", "expired") ~ "Dead",
      TRUE ~ as.character(discharge_category)
    )
  )

# -----------------------------------------------------------------------------
# Verify analysis variables before selecting a canonical hospitalization row
# -----------------------------------------------------------------------------

hospitalization_analysis_conflicts <- clif_hospitalization %>%
  group_by(hospitalization_id) %>%
  summarise(
    n_rows = n(),
    n_patient_ids = n_distinct(patient_id, na.rm = FALSE),
    n_admission_datetimes = n_distinct(admission_dttm, na.rm = FALSE),
    n_discharge_datetimes = n_distinct(discharge_dttm, na.rm = FALSE),
    n_ages = n_distinct(age_at_admission, na.rm = FALSE),
    n_discharge_categories = n_distinct(discharge_category, na.rm = FALSE),
    .groups = "drop"
  ) %>%
  filter(
    n_rows > 1,
    n_patient_ids > 1 | n_admission_datetimes > 1 | n_discharge_datetimes > 1 |
      n_ages > 1 | n_discharge_categories > 1
  )

message(
  "Duplicate hospitalization IDs with conflicts in analysis variables: ",
  nrow(hospitalization_analysis_conflicts)
)

if (nrow(hospitalization_analysis_conflicts) > 0L) {
  print(hospitalization_analysis_conflicts %>% slice_head(n = 20), n = 20, width = Inf)
  stop("Duplicate hospitalization rows disagree on one or more variables used by this analysis.")
}

# -----------------------------------------------------------------------------
# Select one deterministic hospitalization row per normalized ID
# -----------------------------------------------------------------------------

n_hospitalization_rows_before_reconciliation <- nrow(clif_hospitalization)

clif_hospitalization <- clif_hospitalization %>%
  arrange(hospitalization_id, .source_row_order) %>%
  distinct(hospitalization_id, .keep_all = TRUE) %>%
  select(-.source_row_order)

n_hospitalization_rows_after_reconciliation <- nrow(clif_hospitalization)

hospitalization_reconciliation_summary <- tibble(
  diagnostic = c("Hospitalization rows before reconciliation", "Hospitalization rows after reconciliation",
                 "Rows removed during reconciliation", "IDs with case variants before reconciliation",
                 "IDs with conflicts in analysis variables", "Duplicate IDs remaining after reconciliation"),
  n = c(n_hospitalization_rows_before_reconciliation, n_hospitalization_rows_after_reconciliation,
        n_hospitalization_rows_before_reconciliation - n_hospitalization_rows_after_reconciliation,
        hospitalization_id_case_qc$n_ids_with_case_variants, nrow(hospitalization_analysis_conflicts),
        sum(duplicated(clif_hospitalization$hospitalization_id)))
)

message("---- Hospitalization-row reconciliation ----")
print(hospitalization_reconciliation_summary, n = Inf)

if (anyDuplicated(clif_hospitalization$hospitalization_id)) {
  stop("Hospitalization IDs remain duplicated after reconciliation.")
}

if (anyNA(clif_hospitalization$hospitalization_id)) {
  stop("The hospitalization table contains missing hospitalization IDs.")
}

# =============================================================================
# Reconcile patient-table duplicates
# =============================================================================

if (anyNA(clif_patient$patient_id)) {
  stop("The patient table contains missing patient IDs.")
}

n_patient_rows_before_exact_deduplication <- nrow(clif_patient)
clif_patient <- distinct(clif_patient)
n_patient_rows_after_exact_deduplication <- nrow(clif_patient)

duplicate_patient_ids_before_reconciliation <- clif_patient %>%
  count(patient_id, name = "n_rows") %>%
  filter(n_rows > 1)

duplicate_patient_id_values <- duplicate_patient_ids_before_reconciliation$patient_id

message("Patient IDs requiring reconciliation: ", length(duplicate_patient_id_values))

resolve_patient_character <- function(x, conflict_value = NA_character_) {
  values <- as.character(x)
  values <- unique(values[!is.na(values) & str_trim(values) != ""])

  if (length(values) == 0L) {
    return(NA_character_)
  }

  if (length(values) == 1L) {
    return(values[[1]])
  }

  conflict_value
}

patient_reconciliation_values <- clif_patient %>%
  filter(patient_id %in% duplicate_patient_id_values) %>%
  group_by(patient_id) %>%
  summarise(
    death_dttm_conflict = n_distinct(death_dttm, na.rm = TRUE) > 1,
    death_dttm_missingness_difference = any(is.na(death_dttm)) & any(!is.na(death_dttm)),
    language_conflict = n_distinct(language_category, na.rm = TRUE) > 1 |
      n_distinct(language_name, na.rm = TRUE) > 1,
    race_conflict = n_distinct(race_category, na.rm = TRUE) > 1 | n_distinct(race_name, na.rm = TRUE) > 1,
    birth_date_conflict = n_distinct(birth_date, na.rm = TRUE) > 1,
    death_dttm = resolve_patient_character(death_dttm),
    language_category = resolve_patient_character(language_category, "Conflicting source values"),
    language_name = resolve_patient_character(language_name, "Conflicting source values"),
    race_category = resolve_patient_character(race_category, "Multiple/Conflicting"),
    race_name = resolve_patient_character(race_name, "Multiple/Conflicting"),
    birth_date = resolve_patient_character(birth_date),
    .groups = "drop"
  )

patient_duplicate_rows_reconciled <- clif_patient %>%
  filter(patient_id %in% duplicate_patient_id_values) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  select(-any_of(c("death_dttm", "language_category", "language_name", "race_category", "race_name", "birth_date"))) %>%
  left_join(
    patient_reconciliation_values %>%
      select(patient_id, death_dttm, language_category, language_name, race_category, race_name, birth_date),
    by = "patient_id",
    relationship = "one-to-one"
  )

n_patient_rows_before_reconciliation <- nrow(clif_patient)

clif_patient <- bind_rows(
  clif_patient %>% filter(!patient_id %in% duplicate_patient_id_values),
  patient_duplicate_rows_reconciled
)

n_patient_rows_after_reconciliation <- nrow(clif_patient)

duplicate_patient_ids_after_reconciliation <- clif_patient %>%
  count(patient_id, name = "n_rows") %>%
  filter(n_rows > 1)

patient_id_qc <- tibble(
  diagnostic = c("Patient rows before exact deduplication", "Patient rows after exact deduplication",
                 "Exact patient rows removed", "Patient IDs requiring reconciliation",
                 "Patient rows before reconciliation", "Patient rows after reconciliation",
                 "Rows removed during reconciliation", "Patients with conflicting death datetimes",
                 "Patients with missing-versus-populated death datetime", "Patients with conflicting language",
                 "Patients with conflicting race", "Patients with conflicting birth dates",
                 "Duplicate patient IDs remaining"),
  n = c(n_patient_rows_before_exact_deduplication, n_patient_rows_after_exact_deduplication,
        n_patient_rows_before_exact_deduplication - n_patient_rows_after_exact_deduplication,
        length(duplicate_patient_id_values), n_patient_rows_before_reconciliation,
        n_patient_rows_after_reconciliation,
        n_patient_rows_before_reconciliation - n_patient_rows_after_reconciliation,
        sum(patient_reconciliation_values$death_dttm_conflict),
        sum(patient_reconciliation_values$death_dttm_missingness_difference),
        sum(patient_reconciliation_values$language_conflict),
        sum(patient_reconciliation_values$race_conflict),
        sum(patient_reconciliation_values$birth_date_conflict),
        nrow(duplicate_patient_ids_after_reconciliation))
)

message("---- Patient-ID QC and reconciliation ----")
print(patient_id_qc, n = Inf)

if (nrow(duplicate_patient_ids_after_reconciliation) > 0L) {
  stop("Patient IDs remain duplicated after reconciliation.")
}

# =============================================================================
# Identify respiratory intubation anchors independently of procedure evidence
# =============================================================================

required_procedure_columns <- c("hospitalization_id", "procedure_code", "procedure_billed_dttm")
required_resp_columns <- c("hospitalization_id", "recorded_dttm", "support", ".resp_source_row_id")
required_adt_columns <- c("hospitalization_id", "location_category", "in_dttm", "out_dttm")

prefix_missing_columns <- function(required_columns, available_columns, table_name) {
  missing_columns <- setdiff(required_columns, available_columns)

  if (length(missing_columns) == 0L) {
    return(character())
  }

  paste0(table_name, "$", missing_columns)
}

missing_intubation_columns <- c(
  prefix_missing_columns(required_procedure_columns, names(clif_patient_procedures), "patient_procedures"),
  prefix_missing_columns(required_resp_columns, names(clif_respiratory_support), "respiratory_support"),
  prefix_missing_columns(required_adt_columns, names(clif_adt), "adt")
)

if (length(missing_intubation_columns) > 0L) {
  stop("Missing required intubation column(s): ", paste(missing_intubation_columns, collapse = ", "))
}

coerce_to_utc <- function(x) {
  if (inherits(x, "POSIXct")) {
    return(lubridate::with_tz(x, "UTC"))
  }

  suppressWarnings(as.POSIXct(as.character(x), tz = "UTC"))
}

# -----------------------------------------------------------------------------
# Retain valid 31500 procedures as evidence, not as an anchor prerequisite
# -----------------------------------------------------------------------------

intubation_procedures_raw <- clif_patient_procedures %>%
  transmute(
    hospitalization_id,
    procedure_code_clean = str_trim(as.character(procedure_code)),
    procedure_billed_dttm = coerce_to_utc(procedure_billed_dttm)
  ) %>%
  filter(procedure_code_clean == "31500")

intubation_procedures <- intubation_procedures_raw %>%
  filter(!is.na(hospitalization_id), !is.na(procedure_billed_dttm)) %>%
  mutate(procedure_billed_date = as.Date(procedure_billed_dttm, tz = "UTC")) %>%
  distinct(hospitalization_id, procedure_billed_date, .keep_all = TRUE) %>%
  transmute(
    hospitalization_id, procedure_billed_date,
    procedure_billed_dttm = as.POSIXct(procedure_billed_date, tz = "UTC")
  )

intubation_hospitalization_ids <- unique(intubation_procedures$hospitalization_id)

# -----------------------------------------------------------------------------
# Reconcile respiratory states across all hospitalizations
# -----------------------------------------------------------------------------

resp_support_timestamp_rows <- clif_respiratory_support %>%
  transmute(
    hospitalization_id,
    recorded_dttm = coerce_to_utc(recorded_dttm),
    support_recorded = na_if(str_to_lower(str_trim(as.character(support))), "")
  ) %>%
  filter(!is.na(hospitalization_id), !is.na(recorded_dttm))

intubation_processing_scope <- tibble(
  diagnostic = c(
    "Hospitalizations with a valid 31500 procedure", "Respiratory rows in the full table",
    "Respiratory rows retained for procedure-independent processing",
    "Respiratory rows excluded for missing hospitalization ID or recorded datetime"
  ),
  n = c(
    length(intubation_hospitalization_ids), nrow(clif_respiratory_support), nrow(resp_support_timestamp_rows),
    nrow(clif_respiratory_support) - nrow(resp_support_timestamp_rows)
  )
)

resp_support_dt <- as.data.table(resp_support_timestamp_rows)

resp_support_dt[, support_rank := fcase(
  support_recorded == "invasive", 3L,
  support_recorded == "noninvasive", 2L,
  !is.na(support_recorded), 1L,
  default = 0L
)]

setorder(resp_support_dt, hospitalization_id, recorded_dttm, -support_rank, support_recorded)

resp_support_timestamp_dt <- resp_support_dt[, {
  nonmissing_support <- support_recorded[!is.na(support_recorded)]
  distinct_support <- sort(unique(nonmissing_support))

  list(
    n_source_rows = .N,
    n_nonmissing_support_rows = length(nonmissing_support),
    n_distinct_nonmissing_support = length(distinct_support),
    support_values = paste(distinct_support, collapse = "; "),
    selected_support = if (max(support_rank) > 0L) support_recorded[[1L]] else NA_character_
  )
}, by = .(hospitalization_id, recorded_dttm)]

resp_support_timestamp_dt[, `:=`(
  duplicate_timestamp = n_source_rows > 1L,
  conflicting_support = n_distinct_nonmissing_support > 1L,
  support_resolution = fcase(
    n_distinct_nonmissing_support <= 1L, "No nonmissing conflict",
    selected_support == "invasive", "Conflict resolved to invasive",
    selected_support == "noninvasive", "Conflict resolved to noninvasive",
    default = "Conflict resolved to other known support"
  )
)]

resp_support_timestamp_summary <- as_tibble(resp_support_timestamp_dt) %>%
  summarise(
    n_timestamp_keys = n(),
    n_duplicate_timestamp_keys = sum(duplicate_timestamp),
    n_source_rows_at_duplicate_timestamps = sum(n_source_rows[duplicate_timestamp]),
    n_timestamp_keys_with_no_support = sum(n_nonmissing_support_rows == 0L),
    n_timestamp_keys_with_conflicting_support = sum(conflicting_support)
  )

resp_support_conflicting_timestamps <- as_tibble(resp_support_timestamp_dt) %>%
  filter(conflicting_support)

resp_support_conflict_pattern_summary <- resp_support_conflicting_timestamps %>%
  count(support_values, name = "n_timestamp_keys") %>%
  arrange(desc(n_timestamp_keys), support_values)

setorder(resp_support_timestamp_dt, hospitalization_id, recorded_dttm)

resp_support_timestamp_dt[, support_was_missing_at_timestamp := is.na(selected_support)]

resp_support_timestamp_dt[, support_recorded := {
  support_position <- fifelse(!is.na(selected_support), seq_len(.N), 0L)
  last_known_position <- cummax(support_position)
  result <- rep(NA_character_, .N)
  has_known_support <- last_known_position > 0L
  result[has_known_support] <- selected_support[last_known_position[has_known_support]]
  result
}, by = hospitalization_id]

resp_support_timestamp_dt[, support_carried_forward :=
                            support_was_missing_at_timestamp & !is.na(support_recorded)]

resp_support_timestamp_dt[, support_previous := shift(support_recorded), by = hospitalization_id]

resp_support_timestamp_dt[, invasive_transition :=
                            !is.na(support_recorded) & support_recorded == "invasive" &
                            !is.na(support_previous) & support_previous != "invasive"]

resp_support_state_timeline <- resp_support_timestamp_dt[, .(
  hospitalization_id, recorded_dttm, support_at_timestamp = selected_support, support_resolution,
  support_was_missing_at_timestamp, support_recorded, support_carried_forward, support_previous,
  invasive_transition
)] %>%
  as_tibble()

resp_support_resolution_summary <- resp_support_state_timeline %>%
  count(support_resolution, support_carried_forward, name = "n_timestamp_keys") %>%
  mutate(percent = round(100 * n_timestamp_keys / sum(n_timestamp_keys), 3))

intubation_transition_source_qc <- resp_support_state_timeline %>%
  filter(invasive_transition) %>%
  summarise(
    n_invasive_transitions = n(),
    n_directly_documented = sum(!support_carried_forward),
    n_created_from_carried_support = sum(support_carried_forward)
  )

rm(resp_support_dt, resp_support_timestamp_dt, resp_support_timestamp_rows)
invisible(gc())

# -----------------------------------------------------------------------------
# Construct ICU intervals for initially-invasive anchors
# -----------------------------------------------------------------------------
# Entry is inclusive; exit is exclusive.
# A missing exit is filled from hospitalization discharge when available.

icu_intervals_raw <- clif_adt %>%
  filter(str_to_lower(str_trim(as.character(location_category))) == "icu") %>%
  transmute(
    hospitalization_id,
    icu_in_dttm = coerce_to_utc(in_dttm),
    icu_out_dttm_recorded = coerce_to_utc(out_dttm)
  ) %>%
  left_join(
    clif_hospitalization %>%
      transmute(hospitalization_id, hospitalization_discharge_dttm = coerce_to_utc(discharge_dttm)),
    by = "hospitalization_id", relationship = "many-to-one"
  ) %>%
  mutate(
    icu_out_dttm = coalesce(icu_out_dttm_recorded, hospitalization_discharge_dttm),
    exit_filled_from_discharge = is.na(icu_out_dttm_recorded) & !is.na(hospitalization_discharge_dttm),
    valid_interval = !is.na(hospitalization_id) & !is.na(icu_in_dttm) &
      !is.na(icu_out_dttm) & icu_out_dttm > icu_in_dttm
  )

icu_interval_qc <- icu_intervals_raw %>%
  summarise(
    n_icu_rows = n(),
    n_missing_entry = sum(is.na(icu_in_dttm)),
    n_missing_recorded_exit = sum(is.na(icu_out_dttm_recorded)),
    n_exit_filled_from_discharge = sum(exit_filled_from_discharge),
    n_invalid_intervals = sum(!valid_interval),
    n_valid_intervals = sum(valid_interval)
  )

icu_intervals <- icu_intervals_raw %>%
  filter(valid_interval) %>%
  distinct(hospitalization_id, icu_in_dttm, icu_out_dttm)

initial_respiratory_state <- resp_support_state_timeline %>%
  filter(!is.na(support_at_timestamp)) %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  distinct(hospitalization_id, .keep_all = TRUE) %>%
  transmute(hospitalization_id, initial_support_dttm = recorded_dttm, initial_support = support_at_timestamp)

initially_invasive_ids <- initial_respiratory_state %>%
  filter(initial_support == "invasive") %>%
  pull(hospitalization_id)

initially_invasive_icu_candidates <- resp_support_state_timeline %>%
  filter(hospitalization_id %in% initially_invasive_ids, support_at_timestamp == "invasive") %>%
  inner_join(
    icu_intervals,
    by = join_by(hospitalization_id, recorded_dttm >= icu_in_dttm, recorded_dttm < icu_out_dttm),
    relationship = "many-to-many"
  ) %>%
  arrange(hospitalization_id, recorded_dttm, icu_in_dttm, icu_out_dttm) %>%
  distinct(hospitalization_id, recorded_dttm, .keep_all = TRUE)

initially_invasive_icu_anchors <- initially_invasive_icu_candidates %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  distinct(hospitalization_id, .keep_all = TRUE) %>%
  transmute(
    hospitalization_id, intubation_dttm = recorded_dttm,
    support_before_intubation = NA_character_, support_at_intubation = support_at_timestamp,
    intubation_identification_method = "Initially invasive; first documented invasive support in ICU"
  )

transition_anchors <- resp_support_state_timeline %>%
  filter(invasive_transition) %>%
  transmute(
    hospitalization_id, intubation_dttm = recorded_dttm,
    support_before_intubation = support_previous, support_at_intubation = support_at_timestamp,
    intubation_identification_method = "Respiratory transition"
  )

# -----------------------------------------------------------------------------
# Create one event per unique respiratory anchor
# -----------------------------------------------------------------------------
# Billed-date fields are retained as typed missing values for schema compatibility.
# Procedure evidence is hospitalization-level; no billing time anchors an event.

intubation_events <- bind_rows(
  transition_anchors %>% mutate(anchor_priority = 1L),
  initially_invasive_icu_anchors %>% mutate(anchor_priority = 2L)
) %>%
  arrange(hospitalization_id, intubation_dttm, anchor_priority) %>%
  distinct(hospitalization_id, intubation_dttm, .keep_all = TRUE) %>%
  mutate(
    intubation_event_id = row_number(),
    has_31500_in_hospitalization = hospitalization_id %in% intubation_hospitalization_ids,
    procedure_billed_date = as.Date(NA_character_),
    procedure_billed_dttm = as.POSIXct(NA_character_, tz = "UTC")
  ) %>%
  select(-anchor_priority)

initially_invasive_anchor_validation <- initially_invasive_icu_anchors %>%
  left_join(
    initially_invasive_icu_candidates %>%
      select(hospitalization_id, intubation_dttm = recorded_dttm, icu_in_dttm, icu_out_dttm),
    by = c("hospitalization_id", "intubation_dttm"), relationship = "one-to-one"
  ) %>%
  mutate(
    initial_state_invasive = hospitalization_id %in% initially_invasive_ids,
    anchor_within_icu = !is.na(icu_in_dttm) & !is.na(icu_out_dttm) &
      intubation_dttm >= icu_in_dttm & intubation_dttm < icu_out_dttm
  )

procedure_without_anchor <- intubation_procedures %>%
  anti_join(intubation_events %>% distinct(hospitalization_id), by = "hospitalization_id")

intubation_identification_summary <- tibble(
  diagnostic = c(
    "Valid hospitalization-date 31500 procedure events", "Respiratory anchor events",
    "Hospitalizations with a respiratory anchor", "Procedure events in hospitalizations without an anchor"
  ),
  n = c(
    nrow(intubation_procedures), nrow(intubation_events),
    n_distinct(intubation_events$hospitalization_id), nrow(procedure_without_anchor)
  )
)

intubation_identification_method_summary <- intubation_events %>%
  group_by(intubation_identification_method, has_31500_in_hospitalization) %>%
  summarise(n_events = n(), n_hospitalizations = n_distinct(hospitalization_id), .groups = "drop")

intubation_identification_checks <- list(
  respiratory_timestamp_keys_unique =
    anyDuplicated(resp_support_state_timeline[c("hospitalization_id", "recorded_dttm")]) == 0L,
  timestamp_conflicts_quantified =
    nrow(resp_support_conflicting_timestamps) ==
    resp_support_timestamp_summary$n_timestamp_keys_with_conflicting_support,
  timestamp_conflicts_resolved =
    sum(resp_support_state_timeline$support_resolution != "No nonmissing conflict") ==
    nrow(resp_support_conflicting_timestamps),
  transitions_directly_documented = intubation_transition_source_qc$n_created_from_carried_support == 0L,
  event_ids_unique = anyDuplicated(intubation_events$intubation_event_id) == 0L,
  anchor_keys_unique = anyDuplicated(intubation_events[c("hospitalization_id", "intubation_dttm")]) == 0L,
  anchors_have_timestamps = !anyNA(intubation_events$intubation_dttm),
  anchors_documented_invasive = all(intubation_events$support_at_intubation == "invasive"),
  transitions_have_known_prior_noninvasive_state =
    all(!is.na(transition_anchors$support_before_intubation) & transition_anchors$support_before_intubation != "invasive"),
  initially_invasive_state_valid = all(initially_invasive_anchor_validation$initial_state_invasive),
  initially_invasive_anchors_within_icu = all(initially_invasive_anchor_validation$anchor_within_icu),
  carried_support_has_known_value =
    all(!resp_support_state_timeline$support_carried_forward | !is.na(resp_support_state_timeline$support_recorded))
)

if (any(lengths(intubation_identification_checks) != 1L)) {
  stop("At least one intubation integrity check did not return a scalar.")
}

intubation_identification_checks_table <- tibble(
  check = names(intubation_identification_checks),
  passed = as.logical(unlist(intubation_identification_checks, use.names = FALSE))
)

message("---- Respiratory-anchor identification ----")
print(intubation_identification_summary, n = Inf, width = Inf)
print(intubation_identification_method_summary, n = Inf, width = Inf)
print(intubation_identification_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(intubation_identification_checks_table$passed)))

# -----------------------------------------------------------------------------
# Join the expanded resolved support states back to the original respiratory rows
# -----------------------------------------------------------------------------

n_resp_rows_before_state_join <- nrow(clif_respiratory_support)

clif_respiratory_support <- clif_respiratory_support %>%
  select(-any_of(c("support_resolved", "support_at_timestamp", "support_carried_forward"))) %>%
  mutate(recorded_dttm = coerce_to_utc(recorded_dttm)) %>%
  left_join(
    resp_support_state_timeline %>%
      select(
        hospitalization_id, recorded_dttm, support_resolved = support_recorded,
        support_at_timestamp, support_carried_forward
      ),
    by = c("hospitalization_id", "recorded_dttm"), relationship = "many-to-one"
  )

resp_support_state_join_qc <- clif_respiratory_support %>%
  summarise(
    n_source_rows = n(),
    n_rows_missing_recorded_dttm = sum(is.na(recorded_dttm)),
    n_rows_with_resolved_support = sum(!is.na(support_resolved)),
    n_rows_without_resolved_support = sum(is.na(support_resolved)),
    n_rows_with_carried_support = sum(support_carried_forward, na.rm = TRUE)
  )

if (nrow(clif_respiratory_support) != n_resp_rows_before_state_join ||
    anyDuplicated(clif_respiratory_support$.resp_source_row_id)) {
  stop("The support-state join changed the respiratory row count or duplicated source-row identifiers.")
}

# =============================================================================
# Identify tidal volume during the first hour after each intubation event
# =============================================================================

INTUBATION_TIDAL_VOLUME_WINDOW_HOURS <- 1

matched_intubation_events <- intubation_events %>%
  filter(!is.na(intubation_dttm)) %>%
  mutate(intubation_window_end_dttm = intubation_dttm + hours(INTUBATION_TIDAL_VOLUME_WINDOW_HOURS))

respiratory_rows_for_window_matching <- clif_respiratory_support %>%
  filter(hospitalization_id %in% matched_intubation_events$hospitalization_id, !is.na(recorded_dttm))

intubation_window_matches <- matched_intubation_events %>%
  select(
    intubation_event_id, hospitalization_id, procedure_billed_date, procedure_billed_dttm,
    intubation_dttm, intubation_window_end_dttm, intubation_identification_method,
    has_31500_in_hospitalization
  ) %>%
  inner_join(
    respiratory_rows_for_window_matching,
    by = join_by(
      hospitalization_id, intubation_dttm <= recorded_dttm, intubation_window_end_dttm >= recorded_dttm
    ),
    relationship = "many-to-many"
  ) %>%
  mutate(
    hours_after_intubation = as.numeric(difftime(recorded_dttm, intubation_dttm, units = "hours")),
    invasive_at_timepoint = support_resolved == "invasive",
    has_usable_tidal_volume_within_one_hour = !is.na(tidal_volume_set) & tidal_volume_set > 0,
    has_positive_observed_tidal_volume = !is.na(tidal_volume_obs) & tidal_volume_obs > 0,
    tidal_volume_value = case_when(
      invasive_at_timepoint & has_usable_tidal_volume_within_one_hour ~ as.numeric(tidal_volume_set),
      invasive_at_timepoint & !has_usable_tidal_volume_within_one_hour & has_positive_observed_tidal_volume ~
        as.numeric(tidal_volume_obs),
      TRUE ~ NA_real_
    ),
    tidal_volume_source = case_when(
      invasive_at_timepoint & has_usable_tidal_volume_within_one_hour ~ "set",
      invasive_at_timepoint & !has_usable_tidal_volume_within_one_hour & has_positive_observed_tidal_volume ~ "observed",
      TRUE ~ NA_character_
    )
  )

if (any(
  !is.na(intubation_window_matches$tidal_volume_value) &
  intubation_window_matches$invasive_at_timepoint != TRUE,
  na.rm = TRUE
)) {
  stop("A tidal volume was retained without confirmed invasive support.")
}

if (any(
  intubation_window_matches$hours_after_intubation < 0 |
  intubation_window_matches$hours_after_intubation > INTUBATION_TIDAL_VOLUME_WINDOW_HOURS,
  na.rm = TRUE
)) {
  stop("At least one respiratory row fell outside its assigned one-hour intubation window.")
}

intubation_window_event_diagnostic <- matched_intubation_events %>%
  select(
    intubation_event_id, hospitalization_id, procedure_billed_date, procedure_billed_dttm,
    intubation_dttm, intubation_window_end_dttm, intubation_identification_method
  ) %>%
  left_join(
    intubation_window_matches %>%
      group_by(intubation_event_id) %>%
      summarise(
        n_window_rows = n(),
        n_invasive_rows = sum(invasive_at_timepoint, na.rm = TRUE),
        n_noninvasive_rows = sum(!is.na(support_resolved) & support_resolved != "invasive"),
        n_unknown_support_rows = sum(is.na(support_resolved)),
        n_usable_tidal_volume_rows = sum(!is.na(tidal_volume_value)),
        n_set_tidal_volume_rows = sum(tidal_volume_source == "set", na.rm = TRUE),
        n_observed_tidal_volume_rows = sum(tidal_volume_source == "observed", na.rm = TRUE),
        .groups = "drop"
      ),
    by = "intubation_event_id",
    relationship = "one-to-one"
  ) %>%
  mutate(
    across(
      c(
        n_window_rows, n_invasive_rows, n_noninvasive_rows, n_unknown_support_rows,
        n_usable_tidal_volume_rows, n_set_tidal_volume_rows, n_observed_tidal_volume_rows
      ),
      ~ replace_na(.x, 0L)
    ),
    event_tidal_volume_status = case_when(
      n_window_rows == 0L ~ "No respiratory rows in one-hour window",
      n_invasive_rows == 0L ~ "No confirmed invasive respiratory rows in one-hour window",
      n_usable_tidal_volume_rows == 0L ~ "No positive set or observed tidal volume in one-hour window",
      n_set_tidal_volume_rows > 0L ~ "Eligible; positive set tidal volume available",
      n_observed_tidal_volume_rows > 0L ~ "Eligible; observed tidal volume fallback available",
      TRUE ~ "Unexpected"
    )
  )

intubation_method_tidal_volume_summary <- intubation_window_event_diagnostic %>%
  group_by(intubation_identification_method) %>%
  summarise(
    n_identified_events = n(),
    n_with_window_rows = sum(n_window_rows > 0L),
    n_with_invasive_rows = sum(n_invasive_rows > 0L),
    n_with_eligible_tidal_volume = sum(n_usable_tidal_volume_rows > 0L),
    n_without_eligible_tidal_volume = sum(n_usable_tidal_volume_rows == 0L),
    n_with_set_tidal_volume = sum(n_set_tidal_volume_rows > 0L),
    n_with_observed_tidal_volume = sum(n_observed_tidal_volume_rows > 0L),
    .groups = "drop"
  )

intubation_window_overlap <- intubation_window_matches %>%
  count(.resp_source_row_id, name = "n_event_windows") %>%
  filter(n_event_windows > 1L)

intubation_window_overlap_summary <- tibble(
  diagnostic = c(
    "Matched intubation events", "Events with any respiratory row in the one-hour window",
    "Events with a usable tidal volume", "Respiratory source rows matched to multiple event windows",
    "Extra event-window matches attributable to overlap"
  ),
  n = c(
    nrow(matched_intubation_events),
    sum(intubation_window_event_diagnostic$n_window_rows > 0L),
    sum(intubation_window_event_diagnostic$n_usable_tidal_volume_rows > 0L),
    nrow(intubation_window_overlap),
    sum(intubation_window_overlap$n_event_windows - 1L)
  )
)

usable_intubation_window_rows <- intubation_window_matches %>%
  filter(!is.na(tidal_volume_value)) %>%
  arrange(
    .resp_source_row_id, abs(hours_after_intubation), intubation_dttm,
    intubation_event_id
  ) %>%
  distinct(.resp_source_row_id, .keep_all = TRUE)

if (anyDuplicated(usable_intubation_window_rows$.resp_source_row_id)) {
  stop("Respiratory source rows remain duplicated after overlapping intubation windows were reconciled.")
}

tidal_volume_source_summary <- usable_intubation_window_rows %>%
  count(tidal_volume_source, name = "n_respiratory_rows") %>%
  mutate(percent = round(100 * n_respiratory_rows / sum(n_respiratory_rows), 1))

one_hour_tidal_volume_checks <- list(
  matched_event_ids_unique = anyDuplicated(matched_intubation_events$intubation_event_id) == 0L,
  one_event_diagnostic_row_per_matched_event =
    nrow(intubation_window_event_diagnostic) == nrow(matched_intubation_events),
  usable_rows_confirmed_invasive =
    all(usable_intubation_window_rows$support_resolved == "invasive"),
  usable_values_positive = all(usable_intubation_window_rows$tidal_volume_value > 0),
  usable_source_valid = all(usable_intubation_window_rows$tidal_volume_source %in% c("set", "observed")),
  usable_source_rows_unique = anyDuplicated(usable_intubation_window_rows$.resp_source_row_id) == 0L,
  usable_rows_within_window = all(
    usable_intubation_window_rows$hours_after_intubation >= 0 &
      usable_intubation_window_rows$hours_after_intubation <= INTUBATION_TIDAL_VOLUME_WINDOW_HOURS
  )
)

one_hour_tidal_volume_checks_table <- tibble(
  check = names(one_hour_tidal_volume_checks),
  passed = as.logical(unlist(one_hour_tidal_volume_checks, use.names = FALSE))
)

message("---- One-hour tidal-volume integrity checks ----")
print(one_hour_tidal_volume_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(one_hour_tidal_volume_checks_table$passed)))

resp_diagnostic <- intubation_window_event_diagnostic %>%
  group_by(hospitalization_id) %>%
  summarise(
    n_window_intubation_events = n(),
    n_events_with_window_rows = sum(n_window_rows > 0L),
    n_events_with_usable_tidal_volume = sum(n_usable_tidal_volume_rows > 0L),
    has_usable_tidal_volume_within_one_hour = any(n_usable_tidal_volume_rows > 0L),
    has_set_tidal_volume_within_one_hour = any(n_set_tidal_volume_rows > 0L),
    has_observed_fallback_within_one_hour = any(n_observed_tidal_volume_rows > 0L),
    .groups = "drop"
  )

imv_ids <- resp_diagnostic %>%
  filter(has_usable_tidal_volume_within_one_hour) %>%
  pull(hospitalization_id)

message("---- One-hour post-intubation tidal-volume eligibility ----")
print(
  intubation_window_event_diagnostic %>%
    count(event_tidal_volume_status, name = "n_intubation_events"),
  n = Inf, width = Inf
)

message("---- One-hour tidal-volume eligibility by identification method ----")
print(intubation_method_tidal_volume_summary, n = Inf, width = Inf)

message("---- Overlapping intubation-window reconciliation ----")
print(intubation_window_overlap_summary, n = Inf, width = Inf)

message("---- Tidal-volume source in eligible one-hour rows ----")
print(tidal_volume_source_summary, n = Inf, width = Inf)

# =============================================================================
# Summarize ADT eligibility
# =============================================================================

adt_diagnostic <- clif_adt %>%
  group_by(hospitalization_id) %>%
  summarise(
    n_adt_records = n(),
    n_icu_records = sum(str_to_lower(as.character(location_category)) == "icu", na.rm = TRUE),
    has_icu_record = any(str_to_lower(as.character(location_category)) == "icu", na.rm = TRUE),
    all_location_missing = all(is.na(location_category)),
    .groups = "drop"
  )

icu_ids <- adt_diagnostic %>%
  filter(has_icu_record) %>%
  pull(hospitalization_id)

# =============================================================================
# Select the latest eligible hospitalization per patient
# =============================================================================

cohort_patient_and_hospitalization_ids <- clif_hospitalization %>%
  filter(
    age_at_admission >= 3, age_at_admission < 18, hospitalization_id %in% imv_ids,
    hospitalization_id %in% icu_ids
  ) %>%
  arrange(desc(admission_dttm)) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  select(patient_id, hospitalization_id)

eligible_patients <- cohort_patient_and_hospitalization_ids$patient_id
eligible_hospitalizations <- cohort_patient_and_hospitalization_ids$hospitalization_id

if (anyDuplicated(eligible_patients) || anyDuplicated(eligible_hospitalizations)) {
  stop("The selected cohort does not contain exactly one unique hospitalization per patient.")
}

message("Cohort size: ", length(eligible_hospitalizations), " patients/index hospitalizations")

# -----------------------------------------------------------------------------
# Restrict post-intubation tidal-volume rows to the selected cohort
# -----------------------------------------------------------------------------

post_intubation_tidal_volume_output <- usable_intubation_window_rows %>%
  filter(hospitalization_id %in% eligible_hospitalizations)

post_intubation_output_qc <- tibble(
  diagnostic = c(
    "Usable tidal-volume rows before final cohort restriction",
    "Usable tidal-volume rows after final cohort restriction",
    "Final cohort hospitalizations with at least one usable tidal-volume row",
    "Final cohort hospitalizations without a usable tidal-volume row"
  ),
  n = c(
    nrow(usable_intubation_window_rows),
    nrow(post_intubation_tidal_volume_output),
    n_distinct(post_intubation_tidal_volume_output$hospitalization_id),
    sum(!eligible_hospitalizations %in% post_intubation_tidal_volume_output$hospitalization_id)
  )
)

message("---- Post-intubation tidal-volume output scope ----")
print(post_intubation_output_qc, n = Inf, width = Inf)

if (any(!post_intubation_tidal_volume_output$hospitalization_id %in% eligible_hospitalizations)) {
  stop("The post-intubation tidal-volume output contains hospitalization IDs outside the selected cohort.")
}

if (any(!eligible_hospitalizations %in% post_intubation_tidal_volume_output$hospitalization_id)) {
  stop("At least one selected cohort hospitalization lacks a usable post-intubation tidal-volume row.")
}

# Confirm that each hospitalization maps to no more than one patient.
hospitalization_patient_lookup <- clif_hospitalization %>%
  distinct(patient_id, hospitalization_id)

duplicate_hospitalization_mappings <- hospitalization_patient_lookup %>%
  count(hospitalization_id, name = "n_patients") %>%
  filter(n_patients > 1)

if (nrow(duplicate_hospitalization_mappings) > 0L) {
  stop("At least one hospitalization ID maps to multiple patient IDs.")
}

# =============================================================================
# Create combined hospitalization-level cohort diagnostic
# =============================================================================

all_hospitalization_ids <- bind_rows(
  clif_hospitalization %>% distinct(hospitalization_id),
  clif_respiratory_support %>% distinct(hospitalization_id),
  clif_adt %>% distinct(hospitalization_id),
  intubation_procedures %>% distinct(hospitalization_id)
) %>%
  distinct(hospitalization_id)

procedure_diagnostic <- intubation_procedures %>%
  group_by(hospitalization_id) %>%
  summarise(n_31500_events = n(), .groups = "drop")

identified_intubation_diagnostic <- intubation_events %>%
  group_by(hospitalization_id) %>%
  summarise(
    n_respiratory_anchor_events = n(),
    n_identified_intubation_events = sum(!is.na(intubation_dttm)),
    has_identified_intubation = any(!is.na(intubation_dttm)),
    .groups = "drop"
  )

cohort_diagnostic <- all_hospitalization_ids %>%
  left_join(
    clif_hospitalization %>%
      select(patient_id, hospitalization_id, admission_dttm, age_at_admission) %>%
      mutate(in_hospitalization_table = TRUE),
    by = "hospitalization_id",
    relationship = "one-to-one"
  ) %>%
  left_join(procedure_diagnostic, by = "hospitalization_id", relationship = "one-to-one") %>%
  left_join(identified_intubation_diagnostic, by = "hospitalization_id", relationship = "one-to-one") %>%
  left_join(resp_diagnostic, by = "hospitalization_id", relationship = "one-to-one") %>%
  left_join(adt_diagnostic, by = "hospitalization_id", relationship = "one-to-one") %>%
  mutate(
    in_hospitalization_table = replace_na(in_hospitalization_table, FALSE),
    n_31500_events = replace_na(n_31500_events, 0L),
    n_respiratory_anchor_events = replace_na(n_respiratory_anchor_events, 0L),
    n_identified_intubation_events = replace_na(n_identified_intubation_events, 0L),
    n_window_intubation_events = replace_na(n_window_intubation_events, 0L),
    has_identified_intubation = replace_na(has_identified_intubation, FALSE),
    n_events_with_window_rows = replace_na(n_events_with_window_rows, 0L),
    n_events_with_usable_tidal_volume = replace_na(n_events_with_usable_tidal_volume, 0L),
    has_usable_tidal_volume_within_one_hour = replace_na(has_usable_tidal_volume_within_one_hour, FALSE),
    has_set_tidal_volume_within_one_hour = replace_na(has_set_tidal_volume_within_one_hour, FALSE),
    has_observed_fallback_within_one_hour = replace_na(has_observed_fallback_within_one_hour, FALSE),
    n_adt_records = replace_na(n_adt_records, 0L),
    n_icu_records = replace_na(n_icu_records, 0L),
    has_icu_record = replace_na(has_icu_record, FALSE),
    all_location_missing = replace_na(all_location_missing, FALSE),
    age_eligible = !is.na(age_at_admission) & age_at_admission >= 3 & age_at_admission < 18,
    age_status = case_when(
      !in_hospitalization_table ~ "Not present in hospitalization table",
      is.na(age_at_admission) ~ "Excluded; age missing",
      age_at_admission < 3 ~ "Excluded; age under 3",
      age_at_admission >= 18 ~ "Excluded; age 18 or older",
      TRUE ~ "Eligible; age 3 to under 18"
    ),
    intubation_status = case_when(
      !has_identified_intubation ~ "Excluded; no qualifying respiratory intubation anchor",
      n_events_with_window_rows == 0L ~ "Excluded; no respiratory rows within one hour after intubation",
      !has_usable_tidal_volume_within_one_hour ~
        "Excluded; no positive set or observed tidal volume within one hour after intubation",
      has_set_tidal_volume_within_one_hour ~ "Eligible; positive set tidal volume within one hour",
      has_observed_fallback_within_one_hour ~ "Eligible; observed tidal volume fallback within one hour",
      TRUE ~ "Unexpected"
    ),
    icu_status = case_when(
      has_icu_record ~ "Eligible; ICU location recorded",
      n_adt_records == 0L ~ "Excluded; no ADT records",
      all_location_missing ~ "Excluded; all location categories missing",
      TRUE ~ "Excluded; ADT present but no ICU location"
    ),
    fully_eligible = age_eligible & has_usable_tidal_volume_within_one_hour & has_icu_record,
    selected_index_hospitalization = hospitalization_id %in% eligible_hospitalizations,
    cohort_status = case_when(
      !in_hospitalization_table ~ "Not present in hospitalization table",
      !age_eligible ~ age_status,
      !has_usable_tidal_volume_within_one_hour ~ intubation_status,
      !has_icu_record ~ icu_status,
      selected_index_hospitalization ~ "Included; selected index hospitalization",
      TRUE ~ "Eligible; later hospitalization selected"
    )
  )

if (anyDuplicated(cohort_diagnostic$hospitalization_id)) {
  stop("The combined cohort diagnostic contains duplicate hospitalization IDs.")
}

# =============================================================================
# Generate cohort summary tables
# =============================================================================

tidal_volume_status_summary <- cohort_diagnostic %>%
  filter(in_hospitalization_table) %>%
  count(intubation_status, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

icu_status_summary <- cohort_diagnostic %>%
  filter(in_hospitalization_table) %>%
  count(icu_status, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 3))

age_status_summary <- cohort_diagnostic %>%
  filter(in_hospitalization_table) %>%
  group_by(age_status) %>%
  summarise(n_hospitalizations = n(), n_patients = n_distinct(patient_id), .groups = "drop") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

age_eligible_tidal_volume_summary <- cohort_diagnostic %>%
  filter(in_hospitalization_table, age_eligible) %>%
  group_by(intubation_status) %>%
  summarise(n_hospitalizations = n(), n_patients = n_distinct(patient_id), .groups = "drop") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

summarize_funnel_step <- function(data, step) {
  data %>%
    summarise(filter_step = step, n_hospitalizations = n_distinct(hospitalization_id),
              n_patients = n_distinct(patient_id))
}

cohort_filter_diagnostic <- bind_rows(
  cohort_diagnostic %>%
    filter(in_hospitalization_table, age_eligible) %>%
    summarize_funnel_step("1. Age 3 to under 18"),
  cohort_diagnostic %>%
    filter(in_hospitalization_table, age_eligible, has_usable_tidal_volume_within_one_hour) %>%
    summarize_funnel_step("2. Invasive tidal volume within one hour after intubation"),
  cohort_diagnostic %>%
    filter(fully_eligible) %>%
    summarize_funnel_step("3. ICU record"),
  cohort_diagnostic %>%
    filter(selected_index_hospitalization) %>%
    summarize_funnel_step("4. Latest eligible hospitalization per patient")
) %>%
  mutate(
    hospitalizations_removed = lag(n_hospitalizations) - n_hospitalizations,
    patients_removed = lag(n_patients) - n_patients
  )

linkage_summary <- tibble(
  diagnostic = c("ICU IDs in ADT", "ICU IDs also in hospitalization table",
                 "ICU IDs absent from hospitalization table"),
  n_hospitalizations = c(length(icu_ids), sum(icu_ids %in% clif_hospitalization$hospitalization_id),
                         sum(!icu_ids %in% clif_hospitalization$hospitalization_id))
)

message("---- One-hour post-intubation tidal-volume event status ----")
print(
  intubation_window_event_diagnostic %>%
    count(event_tidal_volume_status, name = "n_intubation_events"),
  n = Inf, width = Inf
)

message("---- Tidal volume among age-eligible hospitalizations ----")
print(age_eligible_tidal_volume_summary, n = Inf)

message("---- Cumulative cohort filter ----")
print(cohort_filter_diagnostic, n = Inf)

# =============================================================================
# Initialize diagnostic RDS object
# =============================================================================

cohort_diagnostics <- list(
  hospitalization_level = cohort_diagnostic,
  patient_procedure_load_qc = patient_procedure_load_qc,
  intubation_identification = list(
    configuration = list(
      procedure_required = FALSE,
      anchor_methods = c("Respiratory transition", "Initially invasive; first documented invasive support in ICU"),
      initial_state_definition = "First nonmissing reconciled support state in the hospitalization",
      icu_interval_boundary = "Entry inclusive; exit exclusive",
      missing_icu_exit_rule = "Use hospitalization discharge when available; otherwise exclude interval",
      procedure_evidence = "Valid 31500 presence within the hospitalization; not an event timestamp"
    ),
    processing_scope = intubation_processing_scope,
    event_level = intubation_events,
    identification_summary = intubation_identification_summary,
    method_summary = intubation_identification_method_summary,
    procedure_without_anchor = procedure_without_anchor,
    initial_respiratory_state = initial_respiratory_state,
    initially_invasive_anchor_validation = initially_invasive_anchor_validation,
    icu_interval_summary = icu_interval_qc,
    transition_source_qc = intubation_transition_source_qc,
    respiratory_timestamp_summary = resp_support_timestamp_summary,
    respiratory_timestamp_conflicts = resp_support_conflicting_timestamps,
    respiratory_conflict_patterns = resp_support_conflict_pattern_summary,
    respiratory_resolution_summary = resp_support_resolution_summary,
    integrity_checks = intubation_identification_checks_table,
    one_hour_tidal_volume = list(
      configuration = list(
        window_hours = INTUBATION_TIDAL_VOLUME_WINDOW_HOURS,
        window_inclusive = TRUE,
        support_requirement = "Resolved support equals invasive",
        tidal_volume_preference = "Positive set, otherwise positive observed"
      ),
      event_level = intubation_window_event_diagnostic,
      method_summary = intubation_method_tidal_volume_summary,
      overlap_summary = intubation_window_overlap_summary,
      overlap_detail = intubation_window_overlap,
      source_summary = tidal_volume_source_summary,
      support_state_join = resp_support_state_join_qc,
      output_scope = post_intubation_output_qc,
      checks = one_hour_tidal_volume_checks_table
    )
  ),
  tidal_volume_status_summary = tidal_volume_status_summary,
  age_status_summary = age_status_summary,
  icu_status_summary = icu_status_summary,
  age_eligible_tidal_volume_summary = age_eligible_tidal_volume_summary,
  cumulative_filter_summary = cohort_filter_diagnostic,
  linkage_summary = linkage_summary,
  identifier_qc = list(
    hospitalization_case_summary = hospitalization_id_case_qc,
    hospitalization_reconciliation = hospitalization_reconciliation_summary,
    hospitalization_analysis_conflicts = hospitalization_analysis_conflicts,
    patient_summary = patient_id_qc,
    patient_reconciliation = patient_reconciliation_values
  )
)

# =============================================================================
# Attrition log: steps 1-3
# =============================================================================
# All attrition steps count distinct hospitalizations so the counting unit remains consistent with the final cohort.

funnel <- clif_hospitalization %>%
  filter(age_at_admission >= 3, age_at_admission < 18)

n_age <- n_distinct(funnel$hospitalization_id)

funnel <- funnel %>%
  filter(hospitalization_id %in% icu_ids)

n_icu <- n_distinct(funnel$hospitalization_id)

funnel <- funnel %>%
  filter(hospitalization_id %in% imv_ids)

n_imv <- n_distinct(funnel$hospitalization_id)

attrition <- attrition_init() %>%
  attrition_add(ATTRITION_STEPS[1], n_age) %>%
  attrition_add(ATTRITION_STEPS[2], n_icu, exclusion_reason = "No ICU admission") %>%
  attrition_add(
    ATTRITION_STEPS[3], n_imv,
    exclusion_reason = "No positive set or observed tidal volume during invasive support within one hour after a qualifying respiratory anchor"
  )

attrition_step_1_to_3_checks <- list(
  age_eligible_hospitalizations_nonempty = n_age > 0L,
  icu_count_does_not_exceed_age_count = n_icu <= n_age,
  imv_count_does_not_exceed_icu_count = n_imv <= n_icu,
  partial_attrition_has_three_steps = nrow(attrition) == 3L,
  partial_attrition_labels_match_protocol = identical(attrition$step_label, ATTRITION_STEPS[1:3]),
  partial_attrition_counts_match = identical(as.integer(attrition$n_remaining), c(n_age, n_icu, n_imv)),
  partial_attrition_exclusions_nonnegative = all(attrition$n_excluded[-1] >= 0L)
)

print(attrition, n = Inf, width = Inf)

stopifnot(
  isTRUE(all(unlist(attrition_step_1_to_3_checks)))
)

message(
  "Attrition steps 1-3, hospitalization level: age=", n_age,
  ", ICU=", n_icu, ", invasive ventilation=", n_imv
)

# =============================================================================
# Run respiratory-support waterfall
# =============================================================================

waterfall_utility_path <- file.path(config$utils_dir, "process_resp_support_waterfall.R")

if (!file.exists(waterfall_utility_path)) {
  stop("Respiratory waterfall utility not found: ", waterfall_utility_path)
}

source(waterfall_utility_path)

resp_support_cohort <- clif_respiratory_support %>%
  filter(hospitalization_id %in% eligible_hospitalizations) %>%
  mutate(recorded_dttm = ymd_hms(recorded_dttm, quiet = TRUE))

waterfall_result <- process_resp_support_waterfall(resp_support_cohort, return_diagnostics = TRUE)
resp_waterfall <- waterfall_result$data
waterfall_process_diagnostics <- waterfall_result$diagnostics

rm(waterfall_result)

# =============================================================================
# QC: respiratory-waterfall retention
# =============================================================================

waterfall_diagnostic <- cohort_patient_and_hospitalization_ids %>%
  left_join(
    resp_support_cohort %>% count(hospitalization_id, name = "n_input_rows"),
    by = "hospitalization_id"
  ) %>%
  left_join(
    resp_waterfall %>%
      group_by(hospitalization_id) %>%
      summarise(
        n_output_rows = n(),
        n_real_rows = sum(!is_scaffold, na.rm = TRUE),
        n_invasive_rows = sum(str_to_lower(str_trim(as.character(support))) == "invasive", na.rm = TRUE),
        .groups = "drop"
      ),
    by = "hospitalization_id"
  ) %>%
  mutate(
    across(c(n_input_rows, n_output_rows, n_real_rows, n_invasive_rows), ~ replace_na(.x, 0L)),
    waterfall_status = case_when(
      n_input_rows == 0 ~ "No waterfall input",
      n_output_rows == 0 ~ "No waterfall output",
      n_real_rows == 0 ~ "Only scaffold rows",
      n_invasive_rows == 0 ~ "Output present, but no invasive-support rows",
      TRUE ~ "Invasive-support output present"
    )
  )

if (anyDuplicated(waterfall_diagnostic$hospitalization_id)) {
  stop("The waterfall-retention diagnostic contains duplicate hospitalization IDs.")
}

waterfall_retention_summary <- waterfall_diagnostic %>%
  count(waterfall_status, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

message("---- Respiratory-waterfall retention ----")
print(waterfall_retention_summary, n = Inf)

cohort_diagnostics$waterfall_hospitalization_level <- waterfall_diagnostic
cohort_diagnostics$waterfall_retention_summary <- waterfall_retention_summary
cohort_diagnostics$waterfall_process <- waterfall_process_diagnostics

rm(waterfall_process_diagnostics)

# =============================================================================
# Extract and impute height
# =============================================================================
# Height determines both PBW and PFVC, so provenance is retained for downstream sensitivity analysis.
#
# Source preference:
#   1. clif_vitals at the index hospitalization
#   2. supplemental encounter-level source at the index hospitalization
#   3. nearest in-window height from another hospitalization, when enabled
#
# Values are screened for plausibility before within-encounter aggregation and summarized using the median.
# Cross-encounter imputation is restricted to a growth-aware time window.

HEIGHT_AGGREGATION_MIN_CM <- 30
HEIGHT_AGGREGATION_MAX_CM <- 250
HEIGHT_ALLOW_DONOR_IMPUTATION <- TRUE
HEIGHT_WINDOW_STANDARD_DAYS <- 365
HEIGHT_WINDOW_OLDER_FUTURE_DAYS <- 730
HEIGHT_WINDOW_OLDER_AGE_CUTOFF <- 16

plausible_for_aggregation <- function(x) {
  !is.na(x) & x >= HEIGHT_AGGREGATION_MIN_CM & x <= HEIGHT_AGGREGATION_MAX_CM
}

height_within_window <- function(age_years, donor_gap_days) {
  tolerance <- if_else(
    !is.na(age_years) & age_years >= HEIGHT_WINDOW_OLDER_AGE_CUTOFF & !is.na(donor_gap_days) & donor_gap_days < 0,
    HEIGHT_WINDOW_OLDER_FUTURE_DAYS, HEIGHT_WINDOW_STANDARD_DAYS
  )

  !is.na(donor_gap_days) & abs(donor_gap_days) <= tolerance
}

coerce_to_utc <- function(x) {
  if (inherits(x, "POSIXct")) {
    return(lubridate::with_tz(x, "UTC"))
  }

  suppressWarnings(as.POSIXct(as.character(x), tz = "UTC"))
}

# This object is also used by the weight section below. The added admission timestamp is ignored there.
all_hospitalizations_for_imputation <- clif_hospitalization %>%
  filter(patient_id %in% eligible_patients) %>%
  distinct(patient_id, hospitalization_id, admission_dttm) %>%
  mutate(encounter_admit_dttm = coerce_to_utc(admission_dttm)) %>%
  select(-admission_dttm)

# -----------------------------------------------------------------------------
# Primary source: clif_vitals
# -----------------------------------------------------------------------------

vitals_height_rows <- clif_vitals %>%
  filter(vital_category == "height_cm", !is.na(vital_value)) %>%
  transmute(hospitalization_id, height_cm = suppressWarnings(as.numeric(vital_value))) %>%
  filter(!is.na(height_cm))

vitals_height_by_hospitalization <- vitals_height_rows %>%
  filter(plausible_for_aggregation(height_cm)) %>%
  summarise(vitals_height_cm = median(height_cm), n_vitals_height_values = n(), .by = hospitalization_id)

# -----------------------------------------------------------------------------
# Supplemental source: validated format-specific reader and encounter linkage
# -----------------------------------------------------------------------------
# Reader eligibility requires encounter linkage and no supplied patient-ID conflict.
# Apply the shared plausibility bounds here, before the existing median aggregation.
# Measurement-date text is retained for local provenance only; donor timing continues
# to use hospitalization admission timestamps.

supplemental_height_reader_path <- file.path(config$utils_dir, "read_supplemental_height.R")

if (!file.exists(supplemental_height_reader_path)) {
  stop("Supplemental height reader not found: ", supplemental_height_reader_path)
}

source(supplemental_height_reader_path)

supplemental_height_result <- read_supplemental_height(
  setting = config[["supplemental_height"]],
  hospitalization_lookup = clif_hospitalization %>%
    select(hospitalization_id, patient_id)
)

supplemental_height_rows <- supplemental_height_result$rows

supplemental_height_by_hospitalization <- supplemental_height_rows %>%
  filter(plausible_for_aggregation(height_cm)) %>%
  summarise(
    supplemental_height_cm = median(height_cm),
    n_supplemental_height_values = n(),
    supplemental_height_spread_cm = max(height_cm) - min(height_cm),
    .by = hospitalization_id
  )

supplemental_height_row_diagnostics <- supplemental_height_result$row_diagnostics %>%
  mutate(
    height_within_aggregation_bounds = plausible_for_aggregation(height_cm),
    used_in_encounter_aggregation = eligible_for_aggregation & height_within_aggregation_bounds,
    aggregation_status = case_when(
      !eligible_for_aggregation ~ row_status,
      !height_within_aggregation_bounds ~ "Excluded; height outside shared aggregation bounds",
      TRUE ~ "Included in supplemental encounter median"
    )
  )

supplemental_height_linkage <- bind_rows(
  supplemental_height_result$linkage_summary,
  tibble(
    diagnostic = c(
      "Linked numeric rows excluded by shared height bounds",
      "Supplemental rows included in encounter aggregation",
      "Linked encounters with a plausible supplemental height"
    ),
    n = c(
      sum(
        supplemental_height_row_diagnostics$eligible_for_aggregation &
          !supplemental_height_row_diagnostics$height_within_aggregation_bounds
      ),
      sum(supplemental_height_row_diagnostics$used_in_encounter_aggregation),
      nrow(supplemental_height_by_hospitalization)
    )
  )
)

supplemental_height_checks <- list(
  one_diagnostic_row_per_source_row =
    anyDuplicated(supplemental_height_row_diagnostics$source_row_id) == 0L,
  reader_rows_match_eligible_diagnostics =
    nrow(supplemental_height_rows) ==
    sum(supplemental_height_row_diagnostics$eligible_for_aggregation),
  aggregation_counts_match_used_rows =
    sum(supplemental_height_by_hospitalization$n_supplemental_height_values) ==
    sum(supplemental_height_row_diagnostics$used_in_encounter_aggregation),
  aggregated_encounter_ids_unique =
    anyDuplicated(supplemental_height_by_hospitalization$hospitalization_id) == 0L,
  aggregated_encounters_linked =
    all(
      supplemental_height_by_hospitalization$hospitalization_id %in%
        clif_hospitalization$hospitalization_id
    ),
  aggregated_heights_within_shared_bounds =
    all(plausible_for_aggregation(supplemental_height_by_hospitalization$supplemental_height_cm)),
  used_rows_have_no_patient_conflict =
    all(
      !supplemental_height_row_diagnostics$used_in_encounter_aggregation |
        (
          !supplemental_height_row_diagnostics$patient_identifier_conflict &
            !supplemental_height_row_diagnostics$patient_lookup_conflict
        )
    )
)

supplemental_height_checks_table <- tibble(
  check = names(supplemental_height_checks),
  passed = as.logical(unlist(supplemental_height_checks, use.names = FALSE))
)

message("---- Supplemental height source linkage and aggregation ----")
message("Supplemental height reader status: ", supplemental_height_result$status)
print(supplemental_height_linkage, n = Inf, width = Inf)
print(supplemental_height_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(supplemental_height_checks_table$passed)))

# Local diagnostic detail only; no new export or identifier-copying action.
# This entry is preserved by the existing complete diagnostic RDS save/read-back.
cohort_diagnostics$supplemental_height <- list(
  configuration = list(
    configured = supplemental_height_configured,
    available = supplemental_height_available,
    format = supplemental_height_format,
    path = supplemental_height_path,
    reader_status = supplemental_height_result$status,
    measurement_date_policy = "Retain source text without parsing; precision policy pending",
    linkage_policy = "Encounter linkage required; supplied patient-ID conflicts excluded",
    aggregation_bounds_cm = c(
      minimum = HEIGHT_AGGREGATION_MIN_CM,
      maximum = HEIGHT_AGGREGATION_MAX_CM
    ),
    within_encounter_reducer = "median",
    donor_timing = "Hospitalization admission timestamps; measurement dates do not determine donor windows"
  ),
  row_level = supplemental_height_row_diagnostics,
  linkage_summary = supplemental_height_linkage,
  checks = supplemental_height_checks_table
)

rm(
  supplemental_height_result,
  supplemental_height_rows,
  supplemental_height_row_diagnostics,
  supplemental_height_checks,
  supplemental_height_checks_table
)

# -----------------------------------------------------------------------------
# Encounter-level height from both sources
# -----------------------------------------------------------------------------

encounter_height <- all_hospitalizations_for_imputation %>%
  left_join(vitals_height_by_hospitalization, by = "hospitalization_id", relationship = "one-to-one") %>%
  left_join(supplemental_height_by_hospitalization, by = "hospitalization_id", relationship = "one-to-one") %>%
  mutate(
    across(c(n_vitals_height_values, n_supplemental_height_values), ~ replace_na(.x, 0L)),
    encounter_height_cm = coalesce(vitals_height_cm, supplemental_height_cm)
  )

# -----------------------------------------------------------------------------
# Index-hospitalization height, then optional windowed imputation
# -----------------------------------------------------------------------------

index_encounter_height <- cohort_patient_and_hospitalization_ids %>%
  left_join(
    clif_hospitalization %>%
      transmute(hospitalization_id, index_admit_dttm = coerce_to_utc(admission_dttm), age_at_admission),
    by = "hospitalization_id", relationship = "one-to-one"
  ) %>%
  left_join(
    encounter_height %>%
      select(
        hospitalization_id, index_vitals_height_cm = vitals_height_cm,
        index_supplemental_height_cm = supplemental_height_cm
      ),
    by = "hospitalization_id", relationship = "one-to-one"
  )

donor_height_candidates <- index_encounter_height %>%
  select(patient_id, hospitalization_id, index_admit_dttm, age_at_admission) %>%
  inner_join(
    encounter_height %>%
      select(
        patient_id, donor_hospitalization_id = hospitalization_id,
        donor_admit_dttm = encounter_admit_dttm, donor_height_cm = encounter_height_cm
      ),
    by = "patient_id", relationship = "many-to-many"
  ) %>%
  filter(hospitalization_id != donor_hospitalization_id, plausible_for_aggregation(donor_height_cm)) %>%
  mutate(
    donor_gap_days = as.numeric(difftime(index_admit_dttm, donor_admit_dttm, units = "days")),
    donor_in_window = height_within_window(age_at_admission, donor_gap_days)
  )

donor_datetime_qc <- donor_height_candidates %>%
  summarise(
    n_candidates = n(), n_missing_index_admit_dttm = sum(is.na(index_admit_dttm)),
    n_missing_donor_admit_dttm = sum(is.na(donor_admit_dttm)), n_missing_donor_gap_days = sum(is.na(donor_gap_days))
  )

selected_donor_height <- if (HEIGHT_ALLOW_DONOR_IMPUTATION) {
  donor_height_candidates %>%
    filter(donor_in_window) %>%
    arrange(hospitalization_id, abs(donor_gap_days), donor_hospitalization_id) %>%
    distinct(hospitalization_id, .keep_all = TRUE) %>%
    select(hospitalization_id, donor_height_cm, donor_gap_days, donor_hospitalization_id)
} else {
  tibble(
    hospitalization_id = character(), donor_height_cm = numeric(),
    donor_gap_days = numeric(), donor_hospitalization_id = character()
  )
}

cohort_height_assembled <- index_encounter_height %>%
  left_join(selected_donor_height, by = "hospitalization_id", relationship = "one-to-one") %>%
  mutate(
    height_cm = coalesce(index_vitals_height_cm, index_supplemental_height_cm, donor_height_cm),
    height_source = case_when(
      plausible_for_aggregation(index_vitals_height_cm) ~ "Index encounter, vitals",
      plausible_for_aggregation(index_supplemental_height_cm) ~ "Index encounter, supplemental source",
      plausible_for_aggregation(donor_height_cm) ~ "Another encounter within window",
      TRUE ~ "Unavailable"
    ),
    height_gap_days = if_else(height_source == "Another encounter within window", donor_gap_days, 0)
  )

cohort_heights <- cohort_height_assembled %>%
  filter(!is.na(height_cm)) %>%
  select(hospitalization_id, height_cm, height_source, height_gap_days)

current_height_usable <- cohort_height_assembled %>%
  transmute(
    hospitalization_id, height_cm,
    current_usable = !is.na(height_cm) & height_cm >= 80 & height_cm <= 200
  )

if (anyDuplicated(cohort_heights$hospitalization_id)) {
  stop("The final height table contains duplicate hospitalization IDs.")
}

if (any(cohort_heights$height_source == "Unavailable")) {
  stop("A retained height row is labelled as unavailable.")
}

message(
  "Heights extracted: ", nrow(cohort_heights), " of ",
  nrow(cohort_patient_and_hospitalization_ids), " hospitalizations"
)

# =============================================================================
# QC: height availability, provenance, and source agreement
# =============================================================================

height_source_summary <- cohort_height_assembled %>%
  count(height_source, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

# Counterfactual reproducing the prior implementation: an unfiltered mean within each hospitalization, followed by an
# unwindowed median across all hospitalizations for that patient, using clif_vitals only.
legacy_height_by_hospitalization <- vitals_height_rows %>%
  summarise(legacy_height_cm = mean(height_cm), .by = hospitalization_id)

legacy_patient_median_height <- all_hospitalizations_for_imputation %>%
  left_join(legacy_height_by_hospitalization, by = "hospitalization_id", relationship = "one-to-one") %>%
  summarise(
    legacy_patient_height_cm = if (any(!is.na(legacy_height_cm))) median(legacy_height_cm, na.rm = TRUE) else NA_real_,
    .by = patient_id
  )

legacy_cohort_height <- cohort_patient_and_hospitalization_ids %>%
  left_join(legacy_height_by_hospitalization, by = "hospitalization_id", relationship = "one-to-one") %>%
  left_join(legacy_patient_median_height, by = "patient_id", relationship = "many-to-one") %>%
  mutate(
    legacy_final_height_cm = coalesce(legacy_height_cm, legacy_patient_height_cm),
    legacy_usable = !is.na(legacy_final_height_cm) & legacy_final_height_cm >= 80 & legacy_final_height_cm <= 200
  )

height_recovery_summary <- tibble(
  diagnostic = c(
    "Index hospitalizations", "Usable height, previous implementation", "Usable height, current implementation",
    "Net hospitalizations recovered", "From index vitals", "From supplemental source at index",
    "From another encounter within window", "Still unavailable"
  ),
  n = c(
    nrow(cohort_height_assembled), sum(legacy_cohort_height$legacy_usable), sum(current_height_usable$current_usable),
    sum(current_height_usable$current_usable) - sum(legacy_cohort_height$legacy_usable),
    sum(cohort_height_assembled$height_source == "Index encounter, vitals"),
    sum(cohort_height_assembled$height_source == "Index encounter, supplemental source"),
    sum(cohort_height_assembled$height_source == "Another encounter within window"),
    sum(cohort_height_assembled$height_source == "Unavailable")
  )
)

height_source_agreement <- cohort_height_assembled %>%
  filter(plausible_for_aggregation(index_vitals_height_cm), plausible_for_aggregation(index_supplemental_height_cm)) %>%
  mutate(difference_cm = index_supplemental_height_cm - index_vitals_height_cm)

height_source_agreement_summary <- if (nrow(height_source_agreement) > 0L) {
  height_source_agreement %>%
    summarise(
      n_both_sources = n(), median_difference_cm = median(difference_cm),
      n_within_1cm = sum(abs(difference_cm) <= 1), percent_within_1cm = round(100 * mean(abs(difference_cm) <= 1), 1),
      n_beyond_5cm = sum(abs(difference_cm) > 5), n_beyond_15cm = sum(abs(difference_cm) > 15)
    )
} else {
  tibble(
    n_both_sources = 0L, median_difference_cm = NA_real_, n_within_1cm = 0L,
    percent_within_1cm = NA_real_, n_beyond_5cm = 0L, n_beyond_15cm = 0L
  )
}

height_source_disagreements <- height_source_agreement %>%
  filter(abs(difference_cm) > 5) %>%
  transmute(
    hospitalization_id, age_at_admission, index_vitals_height_cm,
    index_supplemental_height_cm, difference_cm
  ) %>%
  arrange(desc(abs(difference_cm)))

supplemental_height_dispersion_rows <- supplemental_height_by_hospitalization %>%
  semi_join(cohort_patient_and_hospitalization_ids, by = "hospitalization_id") %>%
  filter(n_supplemental_height_values > 1)

supplemental_height_dispersion <- if (nrow(supplemental_height_dispersion_rows) > 0L) {
  supplemental_height_dispersion_rows %>%
    summarise(
      n_encounters_with_multiple_values = n(), median_spread_cm = median(supplemental_height_spread_cm),
      n_spread_within_2cm = sum(supplemental_height_spread_cm <= 2),
      n_spread_beyond_10cm = sum(supplemental_height_spread_cm > 10),
      maximum_spread_cm = max(supplemental_height_spread_cm)
    )
} else {
  tibble(
    n_encounters_with_multiple_values = 0L, median_spread_cm = NA_real_,
    n_spread_within_2cm = 0L, n_spread_beyond_10cm = 0L, maximum_spread_cm = NA_real_
  )
}

donor_growth_plausibility_rows <- donor_height_candidates %>%
  filter(donor_in_window, abs(donor_gap_days) > 30) %>%
  left_join(
    cohort_height_assembled %>%
      transmute(
        hospitalization_id,
        reference_height_cm = coalesce(index_vitals_height_cm, index_supplemental_height_cm)
      ),
    by = "hospitalization_id", relationship = "many-to-one"
  ) %>%
  filter(plausible_for_aggregation(reference_height_cm)) %>%
  mutate(implied_cm_per_year = (reference_height_cm - donor_height_cm) / (donor_gap_days / 365.25))

donor_growth_plausibility <- if (nrow(donor_growth_plausibility_rows) > 0L) {
  donor_growth_plausibility_rows %>%
    summarise(
      n_donor_links_evaluated = n(), median_implied_cm_per_year = median(implied_cm_per_year),
      q25_implied_cm_per_year = quantile(implied_cm_per_year, 0.25, names = FALSE),
      q75_implied_cm_per_year = quantile(implied_cm_per_year, 0.75, names = FALSE),
      n_implying_shrinkage = sum(implied_cm_per_year < -1),
      n_implying_above_12_cm_per_year = sum(implied_cm_per_year > 12)
    )
} else {
  tibble(
    n_donor_links_evaluated = 0L, median_implied_cm_per_year = NA_real_,
    q25_implied_cm_per_year = NA_real_, q75_implied_cm_per_year = NA_real_,
    n_implying_shrinkage = 0L, n_implying_above_12_cm_per_year = 0L
  )
}

height_diagnostic <- cohort_height_assembled %>%
  select(
    patient_id, hospitalization_id, age_at_admission, height_cm, height_source, height_gap_days,
    index_vitals_height_cm, index_supplemental_height_cm, donor_height_cm, donor_hospitalization_id
  )

height_checks <- list(
  one_row_per_index_hospitalization = nrow(cohort_height_assembled) == nrow(cohort_patient_and_hospitalization_ids),
  no_duplicate_height_rows = anyDuplicated(cohort_heights$hospitalization_id) == 0L,
  retained_heights_within_aggregation_bounds =
    nrow(cohort_heights) == 0L || all(plausible_for_aggregation(cohort_heights$height_cm)),
  provenance_complete = !any(is.na(cohort_height_assembled$height_source)),
  recovery_is_nonnegative =
    sum(current_height_usable$current_usable) >= sum(legacy_cohort_height$legacy_usable),
  donor_gaps_respect_window = nrow(selected_donor_height) == 0L ||
    all(
      height_within_window(
        cohort_height_assembled$age_at_admission[
          cohort_height_assembled$height_source == "Another encounter within window"
        ],
        cohort_height_assembled$height_gap_days[
          cohort_height_assembled$height_source == "Another encounter within window"
        ]
      )
    ),
  donor_rows_have_provenance = nrow(selected_donor_height) == 0L ||
    all(!is.na(selected_donor_height$donor_hospitalization_id) & !is.na(selected_donor_height$donor_gap_days)),
  no_index_encounter_used_as_donor = nrow(selected_donor_height) == 0L ||
    all(selected_donor_height$hospitalization_id != selected_donor_height$donor_hospitalization_id)
)

height_checks_table <- tibble(
  check = names(height_checks),
  passed = as.logical(unlist(height_checks, use.names = FALSE))
)

message("---- Height source provenance ----")
print(height_source_summary, n = Inf, width = Inf)

message("---- Height recovery relative to previous implementation ----")
print(height_recovery_summary, n = Inf)

message("---- Agreement between height sources at the index hospitalization ----")
print(height_source_agreement_summary, width = Inf)

if (nrow(height_source_disagreements) > 0L) {
  message("---- Index-hospitalization source disagreements beyond 5 cm ----")
  print(height_source_disagreements, n = Inf, width = Inf)
}

if (supplemental_height_dispersion$n_encounters_with_multiple_values > 0L) {
  message("---- Within-encounter dispersion in the supplemental source ----")
  print(supplemental_height_dispersion, width = Inf)
}

if (donor_growth_plausibility$n_donor_links_evaluated > 0L) {
  message("---- Implied growth velocity for admitted cross-encounter heights ----")
  print(donor_growth_plausibility, width = Inf)
}

message("---- Height integrity checks ----")
print(height_checks_table, n = Inf, width = Inf)

stopifnot(isTRUE(all(height_checks_table$passed)))

cohort_diagnostics$height <- list(
  configuration = list(
    supplemental_source_available = supplemental_height_available,
    supplemental_source_path = if (supplemental_height_available) supplemental_height_path else NA_character_,
    aggregation_bounds_cm = c(minimum = HEIGHT_AGGREGATION_MIN_CM, maximum = HEIGHT_AGGREGATION_MAX_CM),
    analytic_availability_bounds_cm = c(minimum = 80, maximum = 200),
    within_encounter_reducer = "median", donor_imputation_allowed = HEIGHT_ALLOW_DONOR_IMPUTATION,
    window_standard_days = HEIGHT_WINDOW_STANDARD_DAYS,
    window_older_future_days = HEIGHT_WINDOW_OLDER_FUTURE_DAYS,
    window_older_age_cutoff = HEIGHT_WINDOW_OLDER_AGE_CUTOFF,
    source_preference = "Index vitals, then index supplemental source, then nearest in-window encounter"
  ),
  hospitalization_level = height_diagnostic,
  current_analytic_availability = current_height_usable,
  source_summary = height_source_summary,
  recovery_summary = height_recovery_summary,
  source_agreement = height_source_agreement_summary,
  source_disagreements = height_source_disagreements,
  supplemental_linkage = supplemental_height_linkage,
  supplemental_dispersion = supplemental_height_dispersion,
  donor_datetime_qc = donor_datetime_qc,
  donor_growth_plausibility = donor_growth_plausibility,
  checks = height_checks_table
)

# =============================================================================
# Extract and impute weight
# =============================================================================

weight_by_hospitalization <- clif_vitals %>%
  filter(vital_category == "weight_kg", !is.na(vital_value)) %>%
  mutate(weight_kg = suppressWarnings(as.numeric(vital_value))) %>%
  filter(!is.na(weight_kg), weight_kg >= 5, weight_kg <= 250) %>%
  summarise(weight_kg = mean(weight_kg, na.rm = TRUE), .by = hospitalization_id)

patient_weight_availability <- all_hospitalizations_for_imputation %>%
  left_join(
    clif_vitals %>%
      filter(vital_category == "weight_kg") %>%
      transmute(hospitalization_id, weight_numeric = suppressWarnings(as.numeric(vital_value))),
    by = "hospitalization_id"
  ) %>%
  group_by(patient_id) %>%
  summarise(
    has_numeric_weight = any(!is.na(weight_numeric)),
    has_weight_within_5_250 = any(!is.na(weight_numeric) & weight_numeric >= 5 & weight_numeric <= 250),
    .groups = "drop"
  )

patient_median_weight <- all_hospitalizations_for_imputation %>%
  left_join(weight_by_hospitalization, by = "hospitalization_id", relationship = "one-to-one") %>%
  group_by(patient_id) %>%
  summarise(
    median_weight = if_else(any(!is.na(weight_kg)), median(weight_kg, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )

cohort_weights <- cohort_patient_and_hospitalization_ids %>%
  left_join(weight_by_hospitalization, by = "hospitalization_id", relationship = "one-to-one") %>%
  left_join(patient_median_weight, by = "patient_id", relationship = "many-to-one") %>%
  mutate(
    weight_source = case_when(
      !is.na(weight_kg) ~ "Index hospitalization",
      is.na(weight_kg) & !is.na(median_weight) ~ "Another hospitalization",
      TRUE ~ "Unavailable"
    ),
    weight_kg = coalesce(weight_kg, median_weight),
    weight_outside_expected_range = !is.na(weight_kg) & (weight_kg < 9 | weight_kg > 200)
  ) %>%
  select(hospitalization_id, weight_kg, weight_outside_expected_range, weight_source)

if (anyDuplicated(cohort_weights$hospitalization_id)) {
  stop("The final weight table contains duplicate hospitalization IDs.")
}

message(
  "Weights extracted: ", sum(!is.na(cohort_weights$weight_kg)), " of ", nrow(cohort_weights), " hospitalizations"
)

# =============================================================================
# QC: weight availability, imputation, and range
# =============================================================================

weight_diagnostic <- cohort_patient_and_hospitalization_ids %>%
  left_join(cohort_weights, by = "hospitalization_id", relationship = "one-to-one") %>%
  left_join(patient_weight_availability, by = "patient_id", relationship = "many-to-one") %>%
  mutate(
    across(c(has_numeric_weight, has_weight_within_5_250, weight_outside_expected_range), ~ replace_na(.x, FALSE)),
    weight_status = case_when(
      weight_source == "Index hospitalization" ~ "Available in index hospitalization",
      weight_source == "Another hospitalization" ~ "Imputed from another hospitalization",
      has_weight_within_5_250 & is.na(weight_kg) ~ "Unexpected processing loss",
      has_numeric_weight & !has_weight_within_5_250 ~ "Numeric weight outside 5-250 kg bounds",
      TRUE ~ "Unavailable across all hospitalizations"
    )
  )

weight_summary <- weight_diagnostic %>%
  count(weight_status, name = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / sum(n_hospitalizations), 1))

weight_range_summary <- weight_diagnostic %>%
  summarise(
    n_index_hospitalizations = n(),
    n_with_final_weight = sum(!is.na(weight_kg)),
    n_outside_expected_9_200 = sum(weight_outside_expected_range),
    n_without_final_weight = sum(is.na(weight_kg)),
    n_unexpected_processing_loss = sum(weight_status == "Unexpected processing loss")
  )

message("---- Weight availability ----")
print(weight_summary, n = Inf)

message("---- Weight-range QC ----")
print(weight_range_summary, width = Inf)

cohort_diagnostics$weight <- list(
  hospitalization_level = weight_diagnostic,
  availability_summary = weight_summary,
  range_summary = weight_range_summary,
  accepted_bounds_kg = c(minimum = 5, maximum = 250),
  expected_range_kg = c(minimum = 9, maximum = 200)
)

rm(
  vitals_height_rows, vitals_height_by_hospitalization, supplemental_height_by_hospitalization,
  encounter_height, index_encounter_height, donor_height_candidates, selected_donor_height,
  legacy_height_by_hospitalization, legacy_patient_median_height, legacy_cohort_height,
  current_height_usable, supplemental_height_dispersion_rows, donor_growth_plausibility_rows,
  donor_datetime_qc, weight_by_hospitalization, patient_weight_availability, patient_median_weight
)

# =============================================================================
# Extract clinical data
# =============================================================================

cohort_spo2 <- clif_vitals %>%
  filter(hospitalization_id %in% eligible_hospitalizations, vital_category == "spo2") %>%
  mutate(vital_value = suppressWarnings(as.numeric(vital_value))) %>%
  transmute(hospitalization_id, recorded_dttm, spo2_value = vital_value) %>%
  filter(!is.na(spo2_value), spo2_value <= 97)

cohort_map <- clif_vitals %>%
  filter(hospitalization_id %in% eligible_hospitalizations, vital_category == "map") %>%
  mutate(vital_value = suppressWarnings(as.numeric(vital_value))) %>%
  transmute(hospitalization_id, recorded_dttm, map_value = vital_value) %>%
  filter(!is.na(map_value), map_value > 0)

cohort_vitals <- bind_rows(
  cohort_spo2 %>%
    transmute(hospitalization_id, recorded_dttm, vital_category = "spo2", vital_value = spo2_value),
  cohort_map %>%
    transmute(hospitalization_id, recorded_dttm, vital_category = "map", vital_value = map_value)
)

lab_categories_needed <- c("po2_arterial", "pco2_arterial", "creatinine", "bilirubin_total", "platelet_count")

cohort_labs <- clif_labs %>%
  filter(hospitalization_id %in% eligible_hospitalizations, lab_category %in% lab_categories_needed) %>%
  mutate(lab_value_numeric = suppressWarnings(as.numeric(lab_value)))

message("Labs extracted: ", nrow(cohort_labs), " rows")

cohort_meds <- clif_meds %>%
  filter(
    hospitalization_id %in% eligible_hospitalizations,
    med_category %in% c("norepinephrine", "epinephrine", "vasopressin", "dopamine", "phenylephrine", "dobutamine")
  )

cohort_assessments <- clif_assessments %>%
  filter(hospitalization_id %in% eligible_hospitalizations, assessment_category == "gcs_total")

# =============================================================================
# Merge demographics
# =============================================================================

cohort_demographics <- clif_hospitalization %>%
  filter(hospitalization_id %in% eligible_hospitalizations) %>%
  left_join(clif_patient, by = "patient_id", relationship = "many-to-one") %>%
  mutate(
    race_category = case_when(
      str_to_lower(str_trim(race_category)) == "white" ~ "WHITE",
      str_to_lower(str_trim(race_category)) %in% c("black or african american", "black or african-american") ~ "BLACK",
      TRUE ~ "OTHER"
    ),
    deceased = case_when(
      str_to_lower(str_trim(discharge_category)) %in% c("dead", "expired") ~ 1L,
      !is.na(discharge_category) ~ 0L,
      TRUE ~ NA_integer_
    )
  )

if (
  nrow(cohort_demographics) != length(eligible_hospitalizations) ||
  n_distinct(cohort_demographics$hospitalization_id) != length(eligible_hospitalizations)
) {
  stop("The demographics join did not produce exactly one row per eligible hospitalization.")
}

message("Demographics: ", nrow(cohort_demographics), " rows")
message(
  "Mortality rate: ", round(mean(cohort_demographics$deceased, na.rm = TRUE) * 100, 1), "% (",
  sum(is.na(cohort_demographics$deceased)), " missing outcome)"
)

# =============================================================================
# QC: clinical-data availability
# =============================================================================

vitals_qc <- clif_vitals %>%
  filter(hospitalization_id %in% eligible_hospitalizations, vital_category %in% c("spo2", "map")) %>%
  mutate(value_numeric = suppressWarnings(as.numeric(vital_value))) %>%
  group_by(hospitalization_id) %>%
  summarise(
    has_spo2 = any(vital_category == "spo2" & !is.na(value_numeric)),
    has_spo2_at_or_below_97 = any(
      vital_category == "spo2" & !is.na(value_numeric) & value_numeric >= 0 & value_numeric <= 97
    ),
    has_spo2_outside_0_100 = any(
      vital_category == "spo2" & !is.na(value_numeric) & (value_numeric < 0 | value_numeric > 100)
    ),
    has_map = any(vital_category == "map" & !is.na(value_numeric)),
    has_valid_map = any(vital_category == "map" & !is.na(value_numeric) & value_numeric > 0),
    has_nonpositive_map = any(vital_category == "map" & !is.na(value_numeric) & value_numeric <= 0),
    .groups = "drop"
  )

labs_qc <- clif_labs %>%
  filter(hospitalization_id %in% eligible_hospitalizations, lab_category %in% lab_categories_needed) %>%
  mutate(value_numeric = suppressWarnings(as.numeric(lab_value))) %>%
  group_by(hospitalization_id) %>%
  summarise(
    has_pao2 = any(lab_category == "po2_arterial" & !is.na(value_numeric)),
    has_paco2 = any(lab_category == "pco2_arterial" & !is.na(value_numeric)),
    has_creatinine = any(lab_category == "creatinine" & !is.na(value_numeric)),
    has_bilirubin = any(lab_category == "bilirubin_total" & !is.na(value_numeric)),
    has_platelets = any(lab_category == "platelet_count" & !is.na(value_numeric)),
    .groups = "drop"
  )

clinical_data_diagnostic <- cohort_patient_and_hospitalization_ids %>%
  left_join(vitals_qc, by = "hospitalization_id", relationship = "one-to-one") %>%
  left_join(labs_qc, by = "hospitalization_id", relationship = "one-to-one") %>%
  mutate(
    has_vasopressor_record = hospitalization_id %in% cohort_meds$hospitalization_id,
    has_gcs_record = hospitalization_id %in% cohort_assessments$hospitalization_id,
    has_demographics = hospitalization_id %in% cohort_demographics$hospitalization_id
  ) %>%
  left_join(
    cohort_demographics %>% select(hospitalization_id, discharge_category, deceased),
    by = "hospitalization_id",
    relationship = "one-to-one"
  ) %>%
  mutate(across(starts_with("has_"), ~ replace_na(.x, FALSE)))

if (
  nrow(clinical_data_diagnostic) != length(eligible_hospitalizations) ||
  anyDuplicated(clinical_data_diagnostic$hospitalization_id)
) {
  stop("Clinical-data QC did not produce exactly one row per eligible hospitalization.")
}

availability_variables <- c("has_spo2", "has_spo2_at_or_below_97", "has_map", "has_valid_map", "has_pao2",
                            "has_paco2", "has_creatinine", "has_bilirubin", "has_platelets",
                            "has_vasopressor_record", "has_gcs_record", "has_demographics")

clinical_data_availability_summary <- clinical_data_diagnostic %>%
  summarise(across(all_of(availability_variables), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "data_element", values_to = "n_hospitalizations") %>%
  mutate(percent = round(100 * n_hospitalizations / length(eligible_hospitalizations), 1))

clinical_data_quality_summary <- clinical_data_diagnostic %>%
  summarise(
    n_index_hospitalizations = n(),
    n_without_numeric_spo2 = sum(!has_spo2),
    n_with_spo2_but_none_at_or_below_97 = sum(has_spo2 & !has_spo2_at_or_below_97),
    n_with_spo2_outside_0_100 = sum(has_spo2_outside_0_100),
    n_without_valid_map = sum(!has_valid_map),
    n_with_nonpositive_map = sum(has_nonpositive_map),
    n_without_gcs = sum(!has_gcs_record),
    n_without_demographics = sum(!has_demographics),
    n_missing_discharge_category = sum(is.na(discharge_category))
  )

message("---- Clinical-data availability ----")
print(clinical_data_availability_summary, n = Inf)

message("---- Clinical-data quality ----")
print(clinical_data_quality_summary, width = Inf)

cohort_diagnostics$clinical_data <- list(
  hospitalization_level = clinical_data_diagnostic,
  availability_summary = clinical_data_availability_summary,
  quality_summary = clinical_data_quality_summary
)

# =============================================================================
# Save diagnostics
# =============================================================================

normalization_analysis_diagnostics <- cohort_diagnostics

diagnostic_rds_path <- file.path(diagnostic_dir, "cohort_diagnostics.rds")

if (file.exists(diagnostic_rds_path)) {
  prior_diagnostics <- readRDS(diagnostic_rds_path)

  if (!is.list(prior_diagnostics)) {
    stop("The existing diagnostic RDS does not contain a list.")
  }

  if (!is.null(prior_diagnostics$intubation_inclusion_sensitivity)) {
    normalization_analysis_diagnostics$intubation_inclusion_sensitivity <-
      prior_diagnostics$intubation_inclusion_sensitivity
  }

  rm(prior_diagnostics)
}

saveRDS(normalization_analysis_diagnostics, diagnostic_rds_path)

if (!file.exists(diagnostic_rds_path)) {
  stop("The cohort diagnostic RDS was not created: ", diagnostic_rds_path)
}

diagnostic_readback <- readRDS(diagnostic_rds_path)

if (!identical(diagnostic_readback, normalization_analysis_diagnostics)) {
  stop("The Script 01 diagnostic RDS failed its complete read-back check.")
}

if (!identical(diagnostic_readback$height$checks, normalization_analysis_diagnostics$height$checks)) {
  stop("The saved Script 01 height diagnostics failed their read-back check.")
}

message("Cohort diagnostics saved to: ", diagnostic_rds_path)

rm(diagnostic_readback, vitals_qc, labs_qc, availability_variables)

required_post_intubation_columns <- c(
  "intubation_event_id", "hospitalization_id", "intubation_dttm", "recorded_dttm",
  "intubation_identification_method", "tidal_volume_source", "tidal_volume_value"
)

missing_post_intubation_columns <- setdiff(
  required_post_intubation_columns,
  names(post_intubation_tidal_volume_output)
)

if (length(missing_post_intubation_columns) > 0L) {
  stop(
    "The post-intubation tidal-volume output is missing required column(s): ",
    paste(missing_post_intubation_columns, collapse = ", ")
  )
}

allowed_intubation_identification_methods <- c(
  "Respiratory transition",
  "Initially invasive; first documented invasive support in ICU"
)

if (any(is.na(post_intubation_tidal_volume_output$intubation_identification_method))) {
  stop("The post-intubation tidal-volume output contains a missing intubation-identification method.")
}

if (any(
  !post_intubation_tidal_volume_output$intubation_identification_method %in%
  allowed_intubation_identification_methods
)) {
  stop("The post-intubation tidal-volume output contains an unexpected intubation-identification method.")
}

if (any(!post_intubation_tidal_volume_output$tidal_volume_source %in% c("set", "observed"))) {
  stop("The post-intubation tidal-volume output contains an unexpected tidal-volume source.")
}

if (any(
  is.na(post_intubation_tidal_volume_output$tidal_volume_value) |
  post_intubation_tidal_volume_output$tidal_volume_value <= 0
)) {
  stop("The post-intubation tidal-volume output contains a missing or nonpositive selected value.")
}

if (any(
  post_intubation_tidal_volume_output$recorded_dttm < post_intubation_tidal_volume_output$intubation_dttm |
  post_intubation_tidal_volume_output$recorded_dttm >
  post_intubation_tidal_volume_output$intubation_dttm + hours(INTUBATION_TIDAL_VOLUME_WINDOW_HOURS)
)) {
  stop("The post-intubation tidal-volume output contains a row outside its eligible intubation window.")
}

if (anyDuplicated(post_intubation_tidal_volume_output$.resp_source_row_id)) {
  stop("The post-intubation tidal-volume output contains duplicate respiratory source rows.")
}

# =============================================================================
# Save intermediate data
# =============================================================================

saveRDS(eligible_hospitalizations, file.path(output_dir, "cohort_hospitalization_ids.rds"))
write_parquet(resp_waterfall, file.path(output_dir, "resp_support_waterfall.parquet"))
write_parquet(cohort_demographics, file.path(output_dir, "cohort_demographics.parquet"))
write_parquet(cohort_vitals, file.path(output_dir, "cohort_vitals.parquet"))
write_parquet(cohort_labs, file.path(output_dir, "cohort_labs.parquet"))
write_parquet(cohort_meds, file.path(output_dir, "cohort_meds.parquet"))
write_parquet(cohort_assessments, file.path(output_dir, "cohort_assessments.parquet"))
write_parquet(cohort_heights, file.path(output_dir, "cohort_heights.parquet"))
write_parquet(cohort_weights, file.path(output_dir, "cohort_weights.parquet"))
write_parquet(intubation_events, file.path(output_dir, "intubation_events.parquet"))
write_parquet(post_intubation_tidal_volume_output, file.path(output_dir, "post_intubation_tidal_volume_rows.parquet"))

attrition_partial_path <- file.path(output_dir, "attrition_log_partial.csv")

write_csv(attrition, attrition_partial_path)

attrition_partial_readback <- read_csv(attrition_partial_path, show_col_types = FALSE)

attrition_partial_readback_matches <- identical(names(attrition_partial_readback), names(attrition)) &&
  nrow(attrition_partial_readback) == nrow(attrition) &&
  identical(as.integer(attrition_partial_readback$step_order), as.integer(attrition$step_order)) &&
  identical(attrition_partial_readback$step_label, attrition$step_label) &&
  identical(as.integer(attrition_partial_readback$n_remaining), as.integer(attrition$n_remaining)) &&
  identical(as.integer(attrition_partial_readback$n_excluded), as.integer(attrition$n_excluded)) &&
  identical(attrition_partial_readback$exclusion_reason, attrition$exclusion_reason)

attrition_partial_output_checks <- list(
  attrition_partial_file_exists = file.exists(attrition_partial_path),
  attrition_partial_readback_matches = attrition_partial_readback_matches,
  attrition_partial_has_three_steps = nrow(attrition_partial_readback) == 3L,
  attrition_partial_labels_match_protocol = identical(
    attrition_partial_readback$step_label,
    ATTRITION_STEPS[1:3]
  ),
  attrition_partial_counts_nonincreasing = all(diff(attrition_partial_readback$n_remaining) <= 0L),
  attrition_partial_exclusions_nonnegative = all(
    attrition_partial_readback$n_excluded[-1] >= 0,
    na.rm = TRUE
  )
)

attrition_partial_output_qc <- tibble(
  check = names(attrition_partial_output_checks),
  passed = as.logical(unlist(attrition_partial_output_checks, use.names = FALSE))
)

print(attrition_partial_readback, n = Inf, width = Inf)
print(attrition_partial_output_qc, n = Inf, width = Inf)

stopifnot(
  isTRUE(all(attrition_partial_output_qc$passed))
)

message("Partial hospitalization-level attrition log saved to: ", attrition_partial_path)

required_output_files <- c("cohort_hospitalization_ids.rds", "resp_support_waterfall.parquet",
                           "cohort_demographics.parquet", "cohort_vitals.parquet", "cohort_labs.parquet",
                           "cohort_meds.parquet", "cohort_assessments.parquet", "cohort_heights.parquet",
                           "cohort_weights.parquet", "attrition_log_partial.csv", "intubation_events.parquet",
                           "post_intubation_tidal_volume_rows.parquet")

missing_output_files <- required_output_files[
  !file.exists(file.path(output_dir, required_output_files))
]

if (length(missing_output_files) > 0L) {
  stop("Script 01 failed to create required output(s): ", paste(missing_output_files, collapse = ", "))
}

message("All intermediates saved to: ", output_dir)
message("Script 01 complete.")
