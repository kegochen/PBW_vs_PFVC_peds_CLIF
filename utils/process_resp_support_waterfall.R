#' Clean + waterfall-fill the CLIF resp_support table
#'
#' R port of the Python reference pipeline. Builds an hourly scaffold,
#' applies device/mode heuristics, creates hierarchical episode IDs,
#' and performs numeric waterfall fill inside each mode_name_id block.
#'
#' @param resp_support A data.frame / tibble of the raw CLIF respiratory-support
#'   table (timestamps assumed UTC).
#' @param id_col Character. Encounter-level identifier column.
#'   Default \code{"hospitalization_id"}.
#' @param bfill Logical. If TRUE numeric setters are back-filled after
#'   forward-fill; if FALSE (default) only forward-fill is used.
#' @param verbose Logical. Print progress banners when TRUE.
#' @param return_diagnostics Logical. If TRUE, return a list containing the
#'   processed waterfall data and concise diagnostics. If FALSE, return only
#'   the processed waterfall data.
#'   
#' @return If `return_diagnostics = FALSE`, a processed tibble. If TRUE, a list
#'   with elements `data` and `diagnostics`.

library(tidyverse)
library(data.table)
library(lubridate)

process_resp_support_waterfall <- function(resp_support,
                                               id_col = "hospitalization_id",
                                               bfill = FALSE,
                                               verbose = TRUE,
                                               return_diagnostics = FALSE
) {

  p <- if (verbose) message else function(...) invisible(NULL)

  # Helper: forward-fill (+ optional back-fill) a vector
  fb <- function(x) {
    out <- vctrs::vec_fill_missing(x, direction = "down")
    if (bfill) out <- vctrs::vec_fill_missing(out, direction = "up")
    out
  }

  # Helper: cumulative transition IDs calculated separately within each encounter.
  change_id <- function(col, id) {
    dt <- data.table(
      .id = id,
      .val = replace_na(col, "missing")
    )
    
    dt[, .result := cumsum(c(TRUE, .val[-1] != .val[-.N])), by = .id]
    
    as.integer(dt$.result)
  }

  # ------------------------------------------------------------------ #
  # Phase 0 - set-up & hourly scaffold                                 #
  # ------------------------------------------------------------------ #
  p("Phase 0: initialise & create hourly scaffold")
  rs <- as_tibble(resp_support)

  # Lower-case categorical strings
  cat_cols <- intersect(
    c("device_category", "device_name", "mode_category", "mode_name", "support"),
    names(rs)
  )
  rs <- rs |> mutate(across(all_of(cat_cols), tolower))

  # Numeric coercion
  num_cols <- intersect(
    c("tracheostomy", "fio2_set", "lpm_set", "peep_set",
      "tidal_volume_set", "resp_rate_set", "resp_rate_obs",
      "pressure_support_set", "peak_inspiratory_pressure_set", "pressure_control_set"),
    names(rs)
  )
  rs <- rs |> mutate(across(all_of(num_cols), \(x) suppressWarnings(as.numeric(x))))

  # FiO2 scaling: if documented as 40 -> 0.40

  if ("fio2_set" %in% names(rs)) {
    fio2_mean <- mean(rs$fio2_set, na.rm = TRUE)
    if (!is.na(fio2_mean) && fio2_mean > 1.0) {
      rs <- rs |> mutate(fio2_set = if_else(fio2_set > 1, fio2_set / 100, fio2_set))
      p("  Scaled FiO2 values > 1 down by /100")
    }
  }

  # Build hourly scaffold
  p("  Building hourly scaffold")
  bounds <- rs |>
    filter(!is.na(recorded_dttm)) |>
    group_by(.data[[id_col]]) |>
    summarise(
      tmin_h = floor_date(min(recorded_dttm), "hour"),
      tmax_h = floor_date(max(recorded_dttm), "hour"),
      .groups = "drop"
    )

  scaffold <- bounds |>
    rowwise() |>
    reframe(
      !!id_col := .data[[id_col]],
      recorded_dttm = seq(tmin_h, tmax_h, by = "1 hour") + seconds(59 * 60 + 59)
    ) |>
    mutate(
      recorded_date = as.Date(recorded_dttm),
      recorded_hour = hour(recorded_dttm),
      is_scaffold = TRUE
    )

  if (verbose) p(sprintf("  Scaffold rows created: %s", format(nrow(scaffold), big.mark = ",")))

  rs <- rs |>
    mutate(
      recorded_date = as.Date(recorded_dttm),
      recorded_hour = hour(recorded_dttm)
    )

  # ------------------------------------------------------------------ #
  # Phase 1 - heuristic device / mode inference                        #
  # ------------------------------------------------------------------ #
  p("Phase 1: heuristic inference of device & mode")

  # Most-frequent fall-back labels
  device_counts <- rs |>
    filter(!is.na(device_name), !is.na(device_category)) |>
    count(device_name, device_category, sort = TRUE)

  most_common_imv_name <- device_counts |>
    filter(device_category == "imv") |>
    slice(1) |>
    pull(device_name)
  most_common_imv_name <- if (length(most_common_imv_name) == 0) "ventilator" else most_common_imv_name

  most_common_nippv_name <- device_counts |>
    filter(device_category == "nippv") |>
    slice(1) |>
    pull(device_name)
  most_common_nippv_name <- if (length(most_common_nippv_name) == 0) "bipap" else most_common_nippv_name

  mode_counts <- rs |>
    filter(!is.na(mode_name), !is.na(mode_category)) |>
    count(mode_name, mode_category, sort = TRUE)

  most_common_cmv_name <- mode_counts |>
    filter(mode_category == "assist control-volume control") |>
    slice(1) |>
    pull(mode_name)
  most_common_cmv_name <- if (length(most_common_cmv_name) == 0) "AC/VC" else most_common_cmv_name

  # 1-a: IMV from mode_category
  imv_mode_pattern <- "assist control-volume control|simv|pressure control"
  mask_1a <- is.na(rs$device_category) &
    is.na(rs$device_name) &
    str_detect(replace_na(rs$mode_category, ""), imv_mode_pattern)
  rs$device_category[mask_1a] <- "imv"
  rs$device_name[mask_1a] <- most_common_imv_name

  # 1-b: IMV look-behind/ahead
  rs <- rs |> arrange(.data[[id_col]], recorded_dttm)
  rs <- rs |>
    group_by(.data[[id_col]]) |>
    mutate(
      prev_cat = lag(device_category),
      next_cat = lead(device_category)
    ) |>
    ungroup()

  imv_like <- is.na(rs$device_category) &
    (replace_na(rs$prev_cat, "") == "imv" | replace_na(rs$next_cat, "") == "imv") &
    replace_na(rs$peep_set > 1, FALSE) &
    replace_na(rs$resp_rate_set > 1, FALSE) &
    replace_na(rs$tidal_volume_set > 1, FALSE)
  rs$device_category[imv_like] <- "imv"
  rs$device_name[imv_like] <- most_common_imv_name

  # 1-c: NIPPV heuristics (recompute prev/next after 1-b changes)
  rs <- rs |>
    group_by(.data[[id_col]]) |>
    mutate(
      prev_cat = lag(device_category),
      next_cat = lead(device_category)
    ) |>
    ungroup()

  nippv_like <- is.na(rs$device_category) &
    (replace_na(rs$prev_cat, "") == "nippv" | replace_na(rs$next_cat, "") == "nippv") &
    replace_na(rs$peak_inspiratory_pressure_set > 1, FALSE) &
    replace_na(rs$pressure_support_set > 1, FALSE)
  rs$device_category[nippv_like] <- "nippv"
  rs$device_name[nippv_like & is.na(rs$device_name)] <- most_common_nippv_name

  rs <- rs |> select(-prev_cat, -next_cat)

  # 1-d: Clean duplicates & empty rows
  rs <- rs |> arrange(.data[[id_col]], recorded_dttm)
  rs <- rs |>
    group_by(.data[[id_col]], recorded_dttm) |>
    mutate(dup_count = n()) |>
    ungroup() |>
    filter(!(dup_count > 1 & device_category == "nippv"))

  rs <- rs |>
    group_by(.data[[id_col]], recorded_dttm) |>
    mutate(dup_count = n()) |>
    ungroup() |>
    filter(!(dup_count > 1 & is.na(device_category))) |>
    select(-dup_count)

  # 1-e: Nasal-cannula rows must never carry PEEP
  if ("peep_set" %in% names(rs)) {
    mask_bad_nc <- replace_na(rs$device_category == "nasal cannula" & rs$peep_set > 0, FALSE)
    if (any(mask_bad_nc)) {
      rs$device_category[mask_bad_nc] <- NA_character_
      p(sprintf("  %s rows had PEEP>0 on nasal cannula — device_category reset",
                format(sum(mask_bad_nc), big.mark = ",")))
    }
  }

  # 1-f: Infer mode_category from ventilator settings
  #
  # ORDERING MATTERS: each rule only fires on rows where mode_category is still
  # NA, and writes before the next rule is evaluated. Rules are ordered from
  # most-specific to least-specific so that stronger evidence wins:
  #   1. PRVC   (BOTH volume target AND pressure target present)
  #   2. PCV    (pressure_control_set present, no volume target)
  #   2b PCV    (PIP-based; OPTIONAL — see clinical caveat)
  #   3. VCV    (tidal_volume_set present, no pressure target)
  #   3b VCV    (relaxed: TV present, rate missing / PIP tolerated)
  #   4. PS/CPAP (pressure support without a mandatory rate)
  #   4b CPAP   (PEEP only, no PS/rate/volume/pressure target)
  #   5. NIPPV  (settings-based, then broad default)
  #
  # Shared study inference rules; these are not site-adjustable clinical settings.
  #
  # peak_inspiratory_pressure_set must represent a prescribed inspiratory-pressure
  # target. Measured peak pressure must not be mapped into this field.
  # Michigan's prescribed-target interpretation has been confirmed; other sites
  # must verify equivalent source mapping before execution.
  #
  # When mode remains missing, the existing noninvasive-positive-pressure fallback
  # assigns pressure support/cpap. This is inferred mode, not documented mode.
  # The existing predicate includes NIPPV device or noninvasive support labels.
  USE_PIP_FOR_PCV <- TRUE
  USE_NIPPV_DEFAULT <- TRUE
  
  n_missing_mode_before <- sum(is.na(rs$mode_category) &
                                 replace_na(rs$device_category == "imv", FALSE),
                               na.rm = TRUE)
  
  # convenience predicates (all NA-safe) -------------------------------------
  is_imv <- replace_na(
    rs$device_category == "imv" |
      rs$support == "invasive",
    FALSE
  )
  
  is_nippv <- replace_na(
    rs$device_category == "nippv" |
      rs$support == "noninvasive",
    FALSE
  )
  
  has_pc  <- !is.na(rs$pressure_control_set)          & rs$pressure_control_set          > 0
  has_tv  <- !is.na(rs$tidal_volume_set)              & rs$tidal_volume_set              > 0
  has_rr  <- !is.na(rs$resp_rate_set)                 & rs$resp_rate_set                 > 0
  has_ps  <- !is.na(rs$pressure_support_set)          & rs$pressure_support_set          > 0
  has_pip <- !is.na(rs$peak_inspiratory_pressure_set) & rs$peak_inspiratory_pressure_set > 0
  has_peep<- !is.na(rs$peep_set)                      & rs$peep_set                      > 0
  
  no_pc <- is.na(rs$pressure_control_set) | rs$pressure_control_set == 0
  no_tv <- is.na(rs$tidal_volume_set)     | rs$tidal_volume_set     == 0
  no_rr <- is.na(rs$resp_rate_set)        | rs$resp_rate_set        == 0
  no_ps <- is.na(rs$pressure_support_set) | rs$pressure_support_set == 0
  
  # --- 1. Infer PRVC (Pressure-Regulated Volume Control) -------------------- #
  # BOTH tidal_volume_set AND pressure_control_set present: volume target
  # delivered via pressure. Must run FIRST or PCV/VCV would claim these rows.
  infer_prvc <- is.na(rs$mode_category) & is_imv & has_tv & has_pc & has_rr
  rs$mode_category[infer_prvc] <- "pressure-regulated volume control"
  rs$mode_name[infer_prvc & is.na(rs$mode_name)] <- "prvc"
  
  # --- 2. Infer Pressure Control (from pressure_control_set) ---------------- #
  # pressure_control_set is a SET value = definitive PCV indicator, no volume.
  infer_pcv <- is.na(rs$mode_category) & is_imv & has_pc & has_rr & no_tv & no_ps
  rs$mode_category[infer_pcv] <- "pressure control"
  rs$mode_name[infer_pcv & is.na(rs$mode_name)] <- "pressure control"
  
  # --- 2b. Infer Pressure Control from a prescribed PIP target ------------- #
  # Shared inference rule. Source mapping must preserve prescribed-target
  # semantics for peak_inspiratory_pressure_set; measured PIP is not eligible.
  infer_pcv_pip <- rep(FALSE, nrow(rs))
  if (USE_PIP_FOR_PCV) {
    infer_pcv_pip <- is.na(rs$mode_category) & is_imv &
      has_pip & has_rr & no_tv & no_pc
    rs$mode_category[infer_pcv_pip] <- "pressure control"
    rs$mode_name[infer_pcv_pip & is.na(rs$mode_name)] <- "pressure control"
  }
  
  # --- 3. Infer Volume Control (strict) ------------------------------------- #
  # tidal_volume_set present, no pressure target, rate set = classic AC/VC.
  infer_vcv <- is.na(rs$mode_category) & is_imv & has_tv & has_rr & no_pc & no_ps
  rs$mode_category[infer_vcv] <- "assist control-volume control"
  rs$mode_name[infer_vcv & is.na(rs$mode_name)] <- most_common_cmv_name
  
  # --- 3b. Infer Volume Control (relaxed: rate missing) --------------------- #
  # A set tidal volume is strong VCV evidence even if resp_rate_set is absent.
  infer_vcv_norate <- is.na(rs$mode_category) & is_imv & has_tv & no_pc & no_ps
  rs$mode_category[infer_vcv_norate] <- "assist control-volume control"
  rs$mode_name[infer_vcv_norate & is.na(rs$mode_name)] <- most_common_cmv_name
  
  # --- 3c. Infer Volume Control (TV + RR present, PIP tolerated) ------------ #
  # Set TV dominates a measured PIP; still volume targeting (no PC setter).
  infer_vcv_withpip <- is.na(rs$mode_category) & is_imv & has_tv & has_rr & no_pc
  rs$mode_category[infer_vcv_withpip] <- "assist control-volume control"
  rs$mode_name[infer_vcv_withpip & is.na(rs$mode_name)] <- most_common_cmv_name
  
  # --- 4. Infer Pressure Support / CPAP (PS present, no mandatory rate) ----- #
  # Patient triggering all breaths.
  infer_ps <- is.na(rs$mode_category) & is_imv & has_ps & no_rr
  rs$mode_category[infer_ps] <- "pressure support/cpap"
  rs$mode_name[infer_ps & is.na(rs$mode_name)] <- "pressure support"
  
  # --- 4b. Infer Pressure Support (relaxed: PS present, no volume/PC target) #
  # Catches PS rows that also carry a rate but no volume/pressure-control set.
  infer_ps2 <- is.na(rs$mode_category) & is_imv & has_ps & no_tv & no_pc
  rs$mode_category[infer_ps2] <- "pressure support/cpap"
  rs$mode_name[infer_ps2 & is.na(rs$mode_name)] <- "pressure support"
  
  # --- 4c. Infer CPAP (PEEP only, no PS/rate/volume/pressure-control) ------- #
  infer_cpap <- is.na(rs$mode_category) & is_imv &
    has_peep & no_ps & no_rr & no_tv & no_pc
  rs$mode_category[infer_cpap] <- "pressure support/cpap"
  rs$mode_name[infer_cpap & is.na(rs$mode_name)] <- "cpap"
  
  # --- 5. Infer NIPPV mode (settings-based) --------------------------------- #
  # NIPPV device + pressure support OR peak inspiratory pressure documented.
  infer_nippv <- is.na(rs$mode_category) & is_nippv &
    (!is.na(rs$pressure_support_set) | !is.na(rs$peak_inspiratory_pressure_set))
  rs$mode_category[infer_nippv] <- "pressure support/cpap"
  rs$mode_name[infer_nippv & is.na(rs$mode_name)] <- "pressure support"
  
  # --- 5b. Shared broad noninvasive missing-mode fallback ------------------ #
  # Assign an inferred pressure support/cpap category to remaining rows
  # meeting the existing NIPPV-device or noninvasive-support predicate,
  # including rows without documented PS/PIP targets.
  infer_nippv_default <- rep(FALSE, nrow(rs))
  if (USE_NIPPV_DEFAULT) {
    infer_nippv_default <- is.na(rs$mode_category) & is_nippv
    rs$mode_category[infer_nippv_default] <- "pressure support/cpap"
    rs$mode_name[infer_nippv_default & is.na(rs$mode_name)] <- "pressure support"
  }
  
  # --- Report --------------------------------------------------------------- #
  n_missing_mode_after <- sum(is.na(rs$mode_category) &
                                replace_na(rs$device_category == "imv", FALSE),
                              na.rm = TRUE)

  # Drop rows with nothing useful
  all_na_cols <- intersect(
    c("device_category", "device_name", "mode_category", "mode_name",
      "tracheostomy", "fio2_set", "lpm_set", "peep_set", "tidal_volume_set",
      "resp_rate_set", "resp_rate_obs", "pressure_support_set",
      "peak_inspiratory_pressure_set"),
    names(rs)
  )
  rs <- rs |> filter(!if_all(all_of(all_na_cols), is.na))

  # Unique per timestamp
  rs <- rs |> distinct(.data[[id_col]], recorded_dttm, .keep_all = TRUE)

  # Merge scaffold
  rs <- rs |> mutate(is_scaffold = FALSE)
  rs <- bind_rows(rs, scaffold) |>
    arrange(.data[[id_col]], recorded_dttm, recorded_date, recorded_hour)

  # ------------------------------------------------------------------ #
  # Phase 2 - hierarchical IDs                                         #
  # ------------------------------------------------------------------ #
  p("Phase 2: build hierarchical IDs")

  # Forward-fill device_category per encounter
  rs <- rs |>
    arrange(.data[[id_col]], recorded_dttm) |>
    group_by(.data[[id_col]]) |>
    mutate(device_category = vctrs::vec_fill_missing(device_category, direction = "down")) |>
    ungroup()

  rs$device_cat_id <- change_id(rs$device_category, rs[[id_col]])

  # Forward(/back)-fill device_name within (id, device_cat_id)
  rs <- rs |>
    arrange(recorded_dttm) |>
    group_by(.data[[id_col]], device_cat_id) |>
    mutate(device_name = fb(device_name)) |>
    ungroup()

  rs$device_id <- change_id(rs$device_name, rs[[id_col]])

  # Forward(/back)-fill mode_category within (id, device_id)
  rs <- rs |>
    arrange(.data[[id_col]], recorded_dttm) |>
    group_by(.data[[id_col]], device_id) |>
    mutate(mode_category = fb(mode_category)) |>
    ungroup()

  rs$mode_cat_id <- change_id(replace_na(rs$mode_category, "missing"), rs[[id_col]])

  # Forward(/back)-fill mode_name within (id, mode_cat_id)
  rs <- rs |>
    group_by(.data[[id_col]], mode_cat_id) |>
    mutate(mode_name = fb(mode_name)) |>
    ungroup()

  rs$mode_name_id <- change_id(replace_na(rs$mode_name, "missing"), rs[[id_col]])

  # ------------------------------------------------------------------ #
  # Phase 3 - numeric waterfall                                        #
  # ------------------------------------------------------------------ #
  fill_type <- if (bfill) "bi-directional" else "forward-only"
  p(sprintf("Phase 3: %s numeric fill inside mode_name_id blocks", fill_type))

  # FiO2 default for room-air
  if ("fio2_set" %in% names(rs)) {
    rs <- rs |>
      mutate(fio2_set = if_else(
        replace_na(device_category == "room air", FALSE) & is.na(fio2_set),
        0.21, fio2_set
      ))
  }

  # Tidal-volume clean-up
  if ("tidal_volume_set" %in% names(rs)) {
    bad_tv <- (
      (replace_na(rs$mode_category == "pressure support/cpap", FALSE) &
         !is.na(rs$pressure_support_set)) |
      (is.na(rs$mode_category) &
         replace_na(str_detect(rs$device_name, "trach"), FALSE)) |
      (replace_na(rs$mode_category == "pressure support/cpap", FALSE) &
         replace_na(str_detect(rs$device_name, "trach"), FALSE))
    )
    rs$tidal_volume_set[bad_tv] <- NA_real_
  }

  num_cols_fill <- intersect(
    c("fio2_set", "lpm_set", "peep_set", "tidal_volume_set",
      "pressure_support_set", "resp_rate_set", "resp_rate_obs",
      "peak_inspiratory_pressure_set", "pressure_control_set"),
    names(rs)
  )

  n_enc <- n_distinct(rs[[id_col]])
  p(sprintf("  applying waterfall fill to %s encounters", format(n_enc, big.mark = ",")))

  # Waterfall fill within mode_name_id blocks
  # Use data.table for performance
  dt <- as.data.table(rs)

  for (col_name in num_cols_fill) {
    # Within each (id, mode_name_id) block, handle trach collar sub-breaks
    dt[, (col_name) := {
      vals <- get(col_name)
      dc <- device_category
      if (any(dc == "trach collar", na.rm = TRUE)) {
        breaker <- cumsum(replace_na(dc == "trach collar", FALSE))
        unsplit(lapply(split(vals, breaker), fb), breaker)
      } else {
        fb(vals)
      }
    }, by = c(id_col, "mode_name_id")]
  }

  rs <- as_tibble(dt)

  # T-piece -> classify as blow by
  tpiece <- is.na(rs$mode_category) &
    replace_na(str_detect(rs$device_name, "t-piece"), FALSE)
  rs$mode_category[tpiece] <- "blow by"

  # Tracheostomy flag forward-fill per encounter
  if ("tracheostomy" %in% names(rs)) {
    rs <- rs |>
      group_by(.data[[id_col]]) |>
      mutate(tracheostomy = vctrs::vec_fill_missing(tracheostomy, direction = "down")) |>
      ungroup()
  }

  # ------------------------------------------------------------------ #
  # Phase 3b - post-fill mode reconciliation                           #
  # ------------------------------------------------------------------ #
  # After numeric fill, some rows labeled "pressure support/cpap" (often
  # inherited via block fill) carry a SET tidal volume, which PS/CPAP never
  # targets. Reconcile using the co-occurring settings:
  #   (a) TV + mandatory rate + PS  -> SIMV+PS hybrid (mixed flow; MP left NA
  #       downstream by design)
  #   (b) TV + mandatory rate, no PS/PC -> genuine volume control (AC/VC)
  #   (c) TV, no rate -> too ambiguous; left unchanged
  SIMV_MODE_CATEGORY <- "simv"
  SIMV_MODE_NAME     <- "pressure simv"
  
  is_imv_ps <- replace_na(rs$device_category == "imv", FALSE) &
    replace_na(rs$mode_category == "pressure support/cpap", FALSE) &
    !is.na(rs$tidal_volume_set) & rs$tidal_volume_set > 0
  
  reclass_simv <- is_imv_ps &
    !is.na(rs$resp_rate_set)        & rs$resp_rate_set        > 0 &
    !is.na(rs$pressure_support_set) & rs$pressure_support_set > 0
  rs$mode_category[reclass_simv] <- SIMV_MODE_CATEGORY
  rs$mode_name[reclass_simv & is.na(rs$mode_name)] <- SIMV_MODE_NAME
  
  reclass_vcv <- is_imv_ps &
    !is.na(rs$resp_rate_set) & rs$resp_rate_set > 0 &
    (is.na(rs$pressure_support_set) | rs$pressure_support_set == 0) &
    (is.na(rs$pressure_control_set) | rs$pressure_control_set == 0)
  rs$mode_category[reclass_vcv] <- "assist control-volume control"
  rs$mode_name[reclass_vcv & is.na(rs$mode_name)] <- most_common_cmv_name
  
  if (verbose) p(sprintf("  Phase 3b reconciliation: %s PS/CPAP+TV rows -> SIMV, %s -> VCV",
                         format(sum(reclass_simv, na.rm = TRUE), big.mark = ","),
                         format(sum(reclass_vcv,  na.rm = TRUE), big.mark = ",")))
  
  # ------------------------------------------------------------------ #
  # Phase 4 - final tidy-up                                            #
  # ------------------------------------------------------------------ #
  p("Phase 4: final dedup & ordering")
  rs <- rs |>
    distinct() |>
    arrange(.data[[id_col]], recorded_dttm) |>
    select(-any_of(c("recorded_date", "recorded_hour")))

  # ------------------------------------------------------------------ #
  # Phase 5 - concise processing diagnostics
  # ------------------------------------------------------------------ #
  
  mode_inference_summary <- tibble(
    rule = c(
      "PRVC",
      "PCV",
      "PCV from PIP",
      "VCV",
      "VCV without rate",
      "VCV with PIP",
      "Pressure support",
      "Pressure support relaxed",
      "CPAP",
      "NIPPV",
      "NIPPV default",
      "Post-fill reclassified as SIMV",
      "Post-fill reclassified as VCV"
    ),
    n_rows = c(
      sum(infer_prvc, na.rm = TRUE),
      sum(infer_pcv, na.rm = TRUE),
      sum(infer_pcv_pip, na.rm = TRUE),
      sum(infer_vcv, na.rm = TRUE),
      sum(infer_vcv_norate, na.rm = TRUE),
      sum(infer_vcv_withpip, na.rm = TRUE),
      sum(infer_ps, na.rm = TRUE),
      sum(infer_ps2, na.rm = TRUE),
      sum(infer_cpap, na.rm = TRUE),
      sum(infer_nippv, na.rm = TRUE),
      sum(infer_nippv_default, na.rm = TRUE),
      sum(reclass_simv, na.rm = TRUE),
      sum(reclass_vcv, na.rm = TRUE)
    )
  )
  
  mode_missing_summary <- rs %>%
    filter(device_category %in% c("imv", "nippv")) %>%
    mutate(
      row_type = if_else(is_scaffold, "scaffold", "real")
    ) %>%
    group_by(device_category, row_type) %>%
    summarise(
      n_rows = n(),
      n_mode_missing = sum(is.na(mode_category)),
      percent_mode_missing = round(
        100 * mean(is.na(mode_category)),
        1
      ),
      .groups = "drop"
    )
  
  mode_missing_encounter_summary <- rs %>%
    filter(device_category %in% c("imv", "nippv")) %>%
    group_by(.data[[id_col]]) %>%
    summarise(
      any_mode_missing = any(is.na(mode_category)),
      all_modes_missing = all(is.na(mode_category)),
      .groups = "drop"
    ) %>%
    summarise(
      n_encounters = n(),
      n_with_any_missing_mode = sum(any_mode_missing),
      n_with_all_modes_missing = sum(all_modes_missing),
      percent_with_any_missing_mode = round(
        100 * mean(any_mode_missing),
        1
      ),
      percent_with_all_modes_missing = round(
        100 * mean(all_modes_missing),
        1
      )
    )
  
  waterfall_process_summary <- tibble(
    measure = c(
      "Input rows",
      "Output rows",
      "Output encounters",
      "Scaffold rows",
      "Real rows",
      "IMV rows missing mode before Phase1f inference",
      "IMV rows missing mode after Phase 1f inference (before scaffold/fill)"
    ),
    value = c(
      nrow(resp_support),
      nrow(rs),
      n_distinct(rs[[id_col]]),
      sum(rs$is_scaffold, na.rm = TRUE),
      sum(!rs$is_scaffold, na.rm = TRUE),
      n_missing_mode_before,
      n_missing_mode_after
    )
  )
  
  diagnostics <- list(
    process_summary = waterfall_process_summary,
    mode_inference_summary = mode_inference_summary,
    mode_missing_summary = mode_missing_summary,
    mode_missing_encounter_summary =
      mode_missing_encounter_summary
  )
  
  if (verbose) {
    p(sprintf(
      "[OK] Waterfall complete: %s rows across %s encounters",
      format(nrow(rs), big.mark = ","),
      format(n_distinct(rs[[id_col]]), big.mark = ",")
    ))
  }
  
  if (return_diagnostics) {
    return(list(
      data = rs,
      diagnostics = diagnostics
    ))
  }
  
  rs
}
