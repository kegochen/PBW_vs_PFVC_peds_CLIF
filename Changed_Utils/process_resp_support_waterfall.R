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
#'
#' @return A tibble with hourly scaffold rows inserted, device/mode heuristics
#'   applied, hierarchical episode IDs, numeric waterfall fill, tracheostomy
#'   flag forward-filled, and one unique row per (id_col, recorded_dttm).

library(tidyverse)
library(data.table)
library(lubridate)

process_resp_support_waterfall <- function(resp_support,
                                            id_col = "hospitalization_id",
                                            bfill = FALSE,
                                            verbose = TRUE) {

  p <- if (verbose) message else function(...) invisible(NULL)

  # Helper: forward-fill (+ optional back-fill) a vector
  fb <- function(x) {
    out <- vctrs::vec_fill_missing(x, direction = "down")
    if (bfill) out <- vctrs::vec_fill_missing(out, direction = "up")
    out
  }

  # Helper: change-ID — cumulative sum of transitions within groups
  change_id <- function(col, id) {
    filled <- replace_na(col, "missing")
    tibble(.id = id, .val = filled) |>
      mutate(.changed = .val != lag(.val, default = ""),
             .changed = replace_na(.changed, TRUE),
             .block = cumsum(.changed),
             .data.table.aware = NULL,
             .id_change = .val != lag(.val, default = ""),
             .id_change = replace_na(.id_change, TRUE)) |>
      pull(.block) -> blocks
    # Redo properly: per-group cumsum
    dt <- data.table(.id = id, .val = replace_na(col, "missing"))
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
    c("device_category", "device_name", "mode_category", "mode_name"),
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

  # NEW CODE ====================================================================
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
  # TUNABLE FLAG: peak_inspiratory_pressure_set may be a SET pressure-control
  # target OR just a MEASURED peak pressure, depending on your site's data
  # dictionary. If it is only a measured peak, leave USE_PIP_FOR_PCV = FALSE,
  # otherwise VCV rows will be mislabeled as PCV.
  USE_PIP_FOR_PCV     <- TRUE    # set TRUE only if peak_inspiratory_pressure_set is a PCV setter at your site
  USE_NIPPV_DEFAULT   <- TRUE    # broad fallback: NIPPV -> pressure support/cpap
  
  n_missing_mode_before <- sum(is.na(rs$mode_category) &
                                 replace_na(rs$device_category == "imv", FALSE),
                               na.rm = TRUE)
  
  # convenience predicates (all NA-safe) -------------------------------------
  is_imv   <- replace_na(rs$device_category == "imv",   FALSE)
  is_nippv <- replace_na(rs$device_category == "nippv", FALSE)
  
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
  
  # --- 2b. Infer Pressure Control (from PIP) — OPTIONAL --------------------- #
  # Some EHRs record peak_inspiratory_pressure_set instead of
  # pressure_control_set. Only enable if PIP is a SET target at your site.
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
  
  # --- 5b. Infer NIPPV mode (broad default) — OPTIONAL ---------------------- #
  # NIPPV is almost always a pressure-support-family mode; assign a default
  # for remaining NIPPV rows even when PS/PIP were never charted.
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
  
  p(sprintf("  Inferred PRVC:          %s rows", format(sum(infer_prvc,        na.rm = TRUE), big.mark = ",")))
  p(sprintf("  Inferred PCV:           %s rows", format(sum(infer_pcv,         na.rm = TRUE), big.mark = ",")))
  p(sprintf("  Inferred PCV (PIP):     %s rows", format(sum(infer_pcv_pip,     na.rm = TRUE), big.mark = ",")))
  p(sprintf("  Inferred VCV:           %s rows", format(sum(infer_vcv,         na.rm = TRUE), big.mark = ",")))
  p(sprintf("  Inferred VCV (no rate): %s rows", format(sum(infer_vcv_norate,  na.rm = TRUE), big.mark = ",")))
  p(sprintf("  Inferred VCV (w/ PIP):  %s rows", format(sum(infer_vcv_withpip, na.rm = TRUE), big.mark = ",")))
  p(sprintf("  Inferred PS:            %s rows", format(sum(infer_ps,          na.rm = TRUE), big.mark = ",")))
  p(sprintf("  Inferred PS (relaxed):  %s rows", format(sum(infer_ps2,         na.rm = TRUE), big.mark = ",")))
  p(sprintf("  Inferred CPAP:          %s rows", format(sum(infer_cpap,        na.rm = TRUE), big.mark = ",")))
  p(sprintf("  Inferred NIPPV:         %s rows", format(sum(infer_nippv,       na.rm = TRUE), big.mark = ",")))
  p(sprintf("  Inferred NIPPV (dflt):  %s rows", format(sum(infer_nippv_default, na.rm = TRUE), big.mark = ",")))
  p(sprintf("  IMV rows missing mode:  %s -> %s",
            format(n_missing_mode_before, big.mark = ","),
            format(n_missing_mode_after,  big.mark = ",")))
  
  # ================================================================================
  # =========================== START NATIVE CLIF CODE =============================
  # ================================================================================
  
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

  # NEW CODE ====================================================================
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
  SIMV_MODE_CATEGORY <- "SIMV"
  SIMV_MODE_NAME     <- "Pressure SIMV"
  
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
  
  # ================================================================================
  # =========================== START NATIVE CLIF CODE =============================
  # ================================================================================
  
  # ------------------------------------------------------------------ #
  # Phase 4 - final tidy-up                                            #
  # ------------------------------------------------------------------ #
  p("Phase 4: final dedup & ordering")
  rs <- rs |>
    distinct() |>
    arrange(.data[[id_col]], recorded_dttm) |>
    select(-any_of(c("recorded_date", "recorded_hour")))

  # =============================================================================
  # DIAGNOSTIC: Why are mode_category / mode_name still NA?
  # =============================================================================
  p("DIAGNOSTIC: Generating mode NA diagnostics")
  
  mode_diag <- rs %>%
    mutate(
      # --- Why is mode_category NA? ---
      mode_cat_na_reason = case_when(
        !is.na(mode_category) ~ "Not missing",
        
        # Scaffold rows that never got filled
        is_scaffold == TRUE & is.na(device_category) ~
          "Scaffold row: no device category to anchor fill",
        is_scaffold == TRUE & device_category == "room air" ~
          "Scaffold row: room air (no mode expected)",
        is_scaffold == TRUE & device_category == "nasal cannula" ~
          "Scaffold row: nasal cannula (no mode expected)",
        is_scaffold == TRUE & device_category == "high flow nasal cannula" ~
          "Scaffold row: HFNC (no mode expected)",
        is_scaffold == TRUE & device_category == "trach collar" ~
          "Scaffold row: trach collar (no mode expected)",
        is_scaffold == TRUE & device_category %in% c("imv", "nippv") ~
          "Scaffold row: IMV/NIPPV but mode never documented in this block",
        
        # Real (non-scaffold) rows
        is_scaffold == FALSE & is.na(device_category) ~
          "Real row: no device category (unclassified row)",
        is_scaffold == FALSE & device_category == "room air" ~
          "Real row: room air (no mode expected)",
        is_scaffold == FALSE & device_category == "nasal cannula" ~
          "Real row: nasal cannula (no mode expected)",
        is_scaffold == FALSE & device_category == "high flow nasal cannula" ~
          "Real row: HFNC (no mode expected)",
        is_scaffold == FALSE & device_category == "trach collar" ~
          "Real row: trach collar (no mode expected)",
        is_scaffold == FALSE & device_category %in% c("imv", "nippv") ~
          "Real row: IMV/NIPPV but mode was never charted",
        is_scaffold == FALSE & device_category == "nippv" ~
          "Real row: NIPPV but mode not documented",
        
        TRUE ~ "Unknown reason"
      ),
      
      # --- Why is mode_name NA? ---
      mode_name_na_reason = case_when(
        !is.na(mode_name) ~ "Not missing",
        
        # If mode_category is also NA, inherit that reason
        is.na(mode_category) ~
          paste0("mode_category also NA: ", mode_cat_na_reason),
        
        # mode_category exists but mode_name doesn't
        !is.na(mode_category) & is_scaffold == TRUE ~
          "Scaffold row: mode_category filled but mode_name never documented in block",
        !is.na(mode_category) & is_scaffold == FALSE ~
          "Real row: mode_category exists but mode_name was not charted",
        
        TRUE ~ "Unknown reason"
      ),
      
      # Additional context columns for troubleshooting
      has_peep = !is.na(peep_set),
      has_tv = !is.na(tidal_volume_set),
      has_rr = !is.na(resp_rate_set),
      has_pip = !is.na(peak_inspiratory_pressure_set)
    )
  
  # --- Summarize mode_category NA reasons ---
  mode_cat_summary <- mode_diag %>%
    count(mode_cat_na_reason, sort = TRUE) %>%
    mutate(pct = round(n / sum(n) * 100, 1)) %>%
    rename(reason = mode_cat_na_reason)
  
  # --- Summarize mode_name NA reasons ---
  mode_name_summary <- mode_diag %>%
    count(mode_name_na_reason, sort = TRUE) %>%
    mutate(pct = round(n / sum(n) * 100, 1)) %>%
    rename(reason = mode_name_na_reason)
  
  # --- Breakdown by device_category (where are the NAs concentrated?) ---
  mode_cat_by_device <- mode_diag %>%
    filter(is.na(mode_category)) %>%
    count(device_category, is_scaffold, sort = TRUE) %>%
    mutate(pct = round(n / sum(n) * 100, 1))
  
  # --- For IMV/NIPPV rows specifically: what data IS present? ---
  imv_nippv_missing_mode <- mode_diag %>%
    filter(is.na(mode_category),
           device_category %in% c("imv", "nippv")) %>%
    summarise(
      total_rows = n(),
      n_scaffold = sum(is_scaffold),
      n_real = sum(!is_scaffold),
      pct_has_peep = round(mean(has_peep) * 100, 1),
      pct_has_tv = round(mean(has_tv) * 100, 1),
      pct_has_rr = round(mean(has_rr) * 100, 1),
      pct_has_pip = round(mean(has_pip) * 100, 1)
    )
  
  # --- Print to console ---
  p("\n========== MODE_CATEGORY NA SUMMARY ==========")
  print(mode_cat_summary, n = Inf)
  
  p("\n========== MODE_NAME NA SUMMARY ==========")
  print(mode_name_summary, n = Inf)
  
  p("\n========== MODE NA BY DEVICE CATEGORY ==========")
  print(mode_cat_by_device, n = Inf)
  
  p("\n========== IMV/NIPPV ROWS MISSING MODE — WHAT DATA IS PRESENT? ==========")
  print(imv_nippv_missing_mode)
  
  # --- Save to files ---
  diag_out <- list(
    mode_category_reasons = mode_cat_summary,
    mode_name_reasons = mode_name_summary,
    mode_na_by_device = mode_cat_by_device,
    imv_nippv_context = imv_nippv_missing_mode
  )
  
  # Save as single CSV (stacked)
  bind_rows(
    mode_cat_summary %>% mutate(table = "mode_category_reasons"),
    mode_name_summary %>% mutate(table = "mode_name_reasons"),
    mode_cat_by_device %>%
      mutate(reason = paste(device_category, "| scaffold:", is_scaffold),
             table = "mode_na_by_device") %>%
      select(table, reason, n, pct)
  ) %>%
    write_csv(file.path(output_dir, "mode_na_diagnostics.csv"))
  
  # Save detailed row-level file for IMV/NIPPV problem rows (first 1000)
  mode_diag %>%
    filter(is.na(mode_category),
           device_category %in% c("imv", "nippv")) %>%
    select(all_of(id_col), recorded_dttm, device_category, device_name,
           mode_category, mode_name, mode_cat_na_reason,
           is_scaffold, has_peep, has_tv, has_rr, has_pip) %>%
    slice_head(n = 1000) %>%
    write_csv(file.path(output_dir, "mode_na_imv_nippv_detail.csv"))
  
# ==================================================================================
# ============================= END DIAGONSTIC =====================================
# ==================================================================================
  
  p("[OK] Respiratory-support waterfall complete.")
  rs
}
