# Read supplemental height without aggregation, imputation, or output-file creation.
read_supplemental_height <- function(setting, hospitalization_lookup) {
  required_packages <- c("dplyr", "tibble", "readr")
  
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Required supplemental-height dependency is unavailable: ", package)
    }
  }
  
  normalize_id <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x[!is.na(x) & x == ""] <- NA_character_
    x
  }
  
  clean_text <- function(x) {
    x <- trimws(as.character(x))
    missing_value <- !is.na(x) & x %in% c("", "NULL", "null", "NA", "N/A", ".")
    x[missing_value] <- NA_character_
    x
  }
  
  empty_rows <- tibble::tibble(
    hospitalization_id = character(),
    height_cm = numeric()
  )
  
  empty_detail <- tibble::tibble(
    source_row_id = integer(),
    source_format = character(),
    hospitalization_id = character(),
    supplied_patient_id = character(),
    linked_patient_id = character(),
    height_text = character(),
    height_cm = numeric(),
    height_dttm_raw = character(),
    encounter_identifier_conflict = logical(),
    patient_identifier_conflict = logical(),
    encounter_linked = logical(),
    patient_lookup_conflict = logical(),
    eligible_for_aggregation = logical(),
    row_status = character()
  )
  
  unavailable_result <- function(status) {
    list(
      rows = empty_rows,
      row_diagnostics = empty_detail,
      linkage_summary = tibble::tibble(diagnostic = status, n = NA_integer_),
      status = status
    )
  }
  
  if (is.null(setting)) {
    return(unavailable_result("Supplemental height source not configured"))
  }
  
  if (
    !is.list(setting) ||
    is.null(names(setting)) ||
    anyDuplicated(names(setting)) ||
    !setequal(names(setting), c("format", "path"))
  ) {
    stop("Supplemental height setting must contain exactly format and path.")
  }
  
  for (key in c("format", "path")) {
    value <- setting[[key]]
    
    if (
      !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !nzchar(trimws(value))
    ) {
      stop("Supplemental height ", key, " must be one nonempty character value.")
    }
  }
  
  if (!setting$format %in% c("michigan_crosswalk", "height_csv")) {
    stop("Unsupported supplemental height format: ", setting$format)
  }
  
  if (dir.exists(setting$path)) {
    stop("Supplemental height path identifies a directory: ", setting$path)
  }
  
  if (!file.exists(setting$path)) {
    return(unavailable_result("Configured supplemental height file unavailable"))
  }
  
  required_lookup_columns <- c("hospitalization_id", "patient_id")
  missing_lookup_columns <- setdiff(required_lookup_columns, names(hospitalization_lookup))
  
  if (length(missing_lookup_columns) > 0L) {
    stop("Missing hospitalization lookup columns: ", paste(missing_lookup_columns, collapse = ", "))
  }
  
  lookup <- tibble::tibble(
    hospitalization_id = normalize_id(hospitalization_lookup$hospitalization_id),
    linked_patient_id = normalize_id(hospitalization_lookup$patient_id)
  )
  lookup <- dplyr::distinct(lookup)
  
  if (anyNA(lookup$hospitalization_id) || anyDuplicated(lookup$hospitalization_id)) {
    stop("Supplemental height requires a unique, nonmissing hospitalization-ID lookup.")
  }
  
  lookup$encounter_linked <- TRUE
  
  if (setting$format == "michigan_crosswalk") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("arrow is required for the Michigan supplemental height crosswalk.")
    }
    
    raw <- arrow::read_parquet(setting$path)
    required_columns <- "height_cm"
    
    if (!any(c("EncounterID", "encounter_id") %in% names(raw))) {
      stop("Michigan supplemental height requires EncounterID or encounter_id.")
    }
  } else {
    raw <- readr::read_csv(
      setting$path,
      col_types = readr::cols(.default = readr::col_character()),
      na = character(),
      trim_ws = FALSE,
      name_repair = "minimal",
      show_col_types = FALSE,
      progress = FALSE
    )
    
    if (nrow(readr::problems(raw)) > 0L) {
      stop("Supplemental height CSV has parsing problems; correct the source file before use.")
    }
    
    required_columns <- c("encounter_id", "patient_id", "height_cm", "height_dttm")
  }
  
  if (anyDuplicated(names(raw))) {
    stop("Supplemental height source contains duplicate column names.")
  }
  
  missing_columns <- setdiff(required_columns, names(raw))
  
  if (length(missing_columns) > 0L) {
    stop("Missing supplemental height columns: ", paste(missing_columns, collapse = ", "))
  }
  
  id_column <- function(name) {
    if (name %in% names(raw)) {
      normalize_id(raw[[name]])
    } else {
      rep(NA_character_, nrow(raw))
    }
  }
  
  encounter_lower <- id_column("encounter_id")
  encounter_upper <- if (setting$format == "michigan_crosswalk") {
    id_column("EncounterID")
  } else {
    rep(NA_character_, nrow(raw))
  }
  
  encounter_conflict <- !is.na(encounter_upper) & !is.na(encounter_lower) &
    encounter_upper != encounter_lower
  
  if (any(encounter_conflict)) {
    stop(
      "The supplemental height source contains ", sum(encounter_conflict),
      " row(s) where EncounterID and encounter_id disagree."
    )
  }
  
  encounter_id <- dplyr::coalesce(encounter_upper, encounter_lower)
  
  if (setting$format == "michigan_crosswalk" && all(is.na(encounter_id))) {
    stop("The supplemental height source contains no usable encounter identifier.")
  }
  
  patient_lower <- id_column("patient_id")
  patient_upper <- if (setting$format == "michigan_crosswalk") {
    id_column("PatientID")
  } else {
    rep(NA_character_, nrow(raw))
  }
  
  patient_conflict <- !is.na(patient_upper) & !is.na(patient_lower) &
    patient_upper != patient_lower
  
  height_text <- clean_text(raw[["height_cm"]])
  
  # Preserve source representation only; do not infer measurement precision or timezone.
  height_dttm_raw <- if ("height_dttm" %in% names(raw)) {
    as.character(raw[["height_dttm"]])
  } else {
    rep(NA_character_, nrow(raw))
  }
  
  detail <- tibble::tibble(
    source_row_id = seq_len(nrow(raw)),
    source_format = setting$format,
    hospitalization_id = encounter_id,
    supplied_patient_id = dplyr::coalesce(patient_lower, patient_upper),
    height_text = height_text,
    height_cm = suppressWarnings(as.numeric(height_text)),
    height_dttm_raw = height_dttm_raw,
    encounter_identifier_conflict = encounter_conflict,
    patient_identifier_conflict = patient_conflict
  )
  
  detail <- dplyr::left_join(
    detail,
    lookup,
    by = "hospitalization_id",
    relationship = "many-to-one",
    na_matches = "never"
  )
  
  detail <- dplyr::mutate(
    detail,
    encounter_linked = dplyr::coalesce(encounter_linked, FALSE),
    patient_lookup_conflict = encounter_linked &
      !is.na(supplied_patient_id) & !is.na(linked_patient_id) &
      supplied_patient_id != linked_patient_id,
    eligible_for_aggregation = encounter_linked &
      !patient_identifier_conflict & !patient_lookup_conflict & is.finite(height_cm),
    row_status = dplyr::case_when(
      is.na(hospitalization_id) ~ "Excluded; missing encounter identifier",
      !encounter_linked ~ "Excluded; encounter not linked",
      patient_identifier_conflict ~ "Excluded; conflicting supplied patient identifiers",
      patient_lookup_conflict ~ "Excluded; patient disagrees with hospitalization lookup",
      !is.finite(height_cm) ~ "Excluded; missing or nonfinite numeric height",
      TRUE ~ "Linked numeric height; pending existing plausibility screen"
    )
  )
  
  rows <- dplyr::filter(detail, eligible_for_aggregation)
  rows <- dplyr::select(rows, hospitalization_id, height_cm)
  
  linkage_summary <- tibble::tibble(
    diagnostic = c(
      "Supplemental rows read",
      "Rows missing encounter identifier",
      "Rows with an encounter identifier not linked to hospitalization",
      "Rows with conflicting supplied patient identifiers",
      "Rows with patient disagreement against hospitalization lookup",
      "Rows with missing or nonfinite numeric height",
      "Linked numeric rows eligible for existing plausibility screen",
      "Rows with a supplied measurement-date value"
    ),
    n = c(
      nrow(detail),
      sum(is.na(detail$hospitalization_id)),
      sum(!is.na(detail$hospitalization_id) & !detail$encounter_linked),
      sum(detail$patient_identifier_conflict),
      sum(detail$patient_lookup_conflict),
      sum(!is.finite(detail$height_cm)),
      nrow(rows),
      sum(!is.na(clean_text(detail$height_dttm_raw)))
    )
  )
  
  list(
    rows = rows,
    row_diagnostics = detail,
    linkage_summary = linkage_summary,
    status = "Supplemental height source read"
  )
}