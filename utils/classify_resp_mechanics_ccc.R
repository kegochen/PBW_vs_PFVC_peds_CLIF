# Classify custom respiratory-mechanics groups using index-hospitalization POA diagnoses.
# Defines a function only; sourcing does not run classification or write files.
# Source diagnosis codes remain unchanged. Matching uses a normalized copy.
# Code usability is structural, not validation against an annual ICD-10-CM code dictionary.

classify_resp_mechanics_ccc <- function(hospitalization_ids, diagnosis_rows, mapping) {
  required_columns <- c(
    "hospitalization_id", "diagnosis_code", "diagnosis_code_format", "poa_present"
  )
  
  missing_columns <- setdiff(required_columns, names(diagnosis_rows))
  
  if (length(missing_columns) > 0L) {
    stop("Diagnosis input is missing: ", paste(missing_columns, collapse = ", "))
  }
  
  hospitalization_ids <- as.character(hospitalization_ids)
  
  if (anyNA(hospitalization_ids) || anyDuplicated(hospitalization_ids)) {
    stop("Selected hospitalization IDs must be nonmissing and unique.")
  }
  
  prefixes <- mapping$prefixes
  
  if (
    !is.list(prefixes) || length(prefixes) != 6L ||
    is.null(names(prefixes)) || anyDuplicated(names(prefixes))
  ) {
    stop("Mapping must contain six uniquely named groups.")
  }
  
  normalize_code <- function(x) {
    gsub(".", "", toupper(trimws(as.character(x))), fixed = TRUE)
  }
  
  prefixes <- lapply(prefixes, normalize_code)
  
  if (any(vapply(prefixes, function(x) anyNA(x) || any(!nzchar(x)), logical(1)))) {
    stop("Mapping contains missing or empty code prefixes.")
  }
  
  rows <- diagnosis_rows[
    diagnosis_rows$hospitalization_id %in% hospitalization_ids &
      !is.na(diagnosis_rows$poa_present) & diagnosis_rows$poa_present == 1,
    required_columns,
    drop = FALSE
  ]
  
  code <- normalize_code(rows$diagnosis_code)
  format_supported <- !is.na(rows$diagnosis_code_format) &
    rows$diagnosis_code_format == "ICD10CM"
  
  code_structurally_usable <- !is.na(code) & grepl("^[A-Z][0-9][A-Z0-9]{1,5}$", code)
  usable <- format_supported & code_structurally_usable
  
  result <- data.frame(
    hospitalization_id = hospitalization_ids,
    stringsAsFactors = FALSE
  )
  
  row_index <- match(rows$hospitalization_id, hospitalization_ids)
  
  n_poa_rows <- tabulate(row_index, nbins = length(hospitalization_ids))
  n_usable_rows <- tabulate(row_index[usable], nbins = length(hospitalization_ids))
  n_gap_rows <- tabulate(row_index[!usable], nbins = length(hospitalization_ids))
  
  absence_evaluable <- n_usable_rows > 0L & n_gap_rows == 0L
  
  for (group in names(prefixes)) {
    positive_row <- rep(FALSE, nrow(rows))
    
    for (prefix in prefixes[[group]]) {
      positive_row[usable] <- positive_row[usable] | startsWith(code[usable], prefix)
    }
    
    positive_hospitalization <- tabulate(
      row_index[positive_row],
      nbins = length(hospitalization_ids)
    ) > 0L
    
    indicator <- rep(NA_integer_, length(hospitalization_ids))
    indicator[absence_evaluable] <- 0L
    indicator[positive_hospitalization] <- 1L
    
    result[[group]] <- indicator
  }
  
  availability <- data.frame(
    hospitalization_id = hospitalization_ids,
    n_poa_diagnosis_rows = n_poa_rows,
    n_usable_poa_icd10cm_rows = n_usable_rows,
    n_poa_coding_gap_rows = n_gap_rows,
    absence_evaluable = absence_evaluable,
    stringsAsFactors = FALSE
  )
  
  distribution <- do.call(
    rbind,
    lapply(names(prefixes), function(group) {
      indicator <- result[[group]]
      
      data.frame(
        group = group,
        n_present = sum(indicator == 1L, na.rm = TRUE),
        n_absent = sum(indicator == 0L, na.rm = TRUE),
        n_missing = sum(is.na(indicator)),
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(distribution) <- NULL
  
  source_summary <- data.frame(
    diagnostic = c(
      "Selected hospitalizations",
      "Selected hospitalizations with POA diagnosis rows",
      "Selected hospitalizations with usable POA ICD10CM rows",
      "Selected hospitalizations with POA coding gaps",
      "Selected hospitalizations with evaluable absence",
      "POA diagnosis rows",
      "Usable POA ICD10CM diagnosis rows",
      "POA rows with missing or unsupported code format",
      "POA rows with structurally unusable diagnosis codes"
    ),
    n = c(
      length(hospitalization_ids),
      sum(n_poa_rows > 0L),
      sum(n_usable_rows > 0L),
      sum(n_gap_rows > 0L),
      sum(absence_evaluable),
      nrow(rows),
      sum(usable),
      sum(!format_supported),
      sum(!code_structurally_usable)
    ),
    stringsAsFactors = FALSE
  )
  
  list(
    data = result,
    availability = availability,
    diagnostics = list(
      source_summary = source_summary,
      distribution = distribution
    ),
    provenance = list(
      mapping = mapping,
      analysis_unit = "Selected index hospitalization",
      diagnosis_scope = "poa_present == 1; all diagnosis positions",
      supported_source_format = "ICD10CM",
      procedure_codes_used = FALSE,
      code_normalization = "Trim, uppercase, and remove decimal points in a matching copy",
      code_usability = "Structural pattern only; no annual code-dictionary validation",
      positive_rule = "Observed qualifying evidence remains 1 despite other coding gaps",
      negative_rule = "0 only with usable POA diagnoses and no POA coding gaps",
      missing_rule = "NA when absence cannot be established",
      multiple_group_membership = TRUE
    )
  )
}