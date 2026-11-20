# Pediatric CLIF: PBW versus PFVC

Single-site R analysis and diagnostic reporting for a pediatric Common
Longitudinal ICU data Format (CLIF) study comparing predicted body weight
(PBW) and predicted forced vital capacity (PFVC) normalization of respiratory
physiology measures.

See the [Executive summary](EXECUTIVE_SUMMARY.md) for the study question,
cohort definitions, and analysis workflow.

Version 3a adds study-specific chronic-condition adjustment to the pediatric
PBW-versus-PFVC analysis. Version 3 is preserved as the baseline at commit
`26e023228ebe6a27c2cb8364c97921b051f81420`.

The updated analysis completed on Windows, followed by separate rendering and
verified publication of all three diagnostic PDFs. Reporting-status
synchronization and the normalization report-inclusive manifest passed the
renderer checks. All three reports were visually reviewed.

Local version 3a execution used the existing baseline project library through
an explicit development override. A separate version 3a library restoration
has not been tested. The prior version 3 short-path Windows restoration
validation remains baseline evidence, not a new version 3a setup test.

The chronic-condition groups are customized from pediatric Complex Chronic
Conditions version 3, with asthma added separately. Classification uses
present-on-admission diagnoses from the selected index hospitalization.
Groups may overlap, and unavailable diagnosis information is not recoded as
condition absence.

See [CCC adjustment policy](reference-data/ccc_adjustment_policy.md) for group
definitions, fixed model-family term sets, and missingness. Individual-modifier
interactions remain exploratory, with unadjusted P values across related
specifications. See [Deferred pooled analyses](reference-data/ccc_pooled_candidates.md)
for candidates not enabled for execution.

The runner uses separate R subprocesses and collects unused parent-process
memory after each successful analysis script. Reporting controls nested R
startup and library routing and calls `gc()` between validated report renders.

Version 3a is published on `main` at commit
`feaa203408e0e8688213c548590b3c3612fa7d34`. Version 3 is preserved on branch
`v3`. The repository contains shared code, configuration templates, references,
and the lockfile; site data and generated results remain local.

The [Neurologic interaction rationale](reference-data/ccc_interaction_rationale.md)
documents the clinical motivation and limits of that exploratory analysis.
Buddy-site testing remains deferred. Pooled exploratory analyses and the
separate growth-method comparison remain deferred until buddy-site testing
is complete.

## Repository structure

- `code/00_run_pipeline_v3a.R`: function-based local workflow runner.
- `code/01_cohort_identification_v3a.R`: cohort identification.
- `code/02_quality_checks_v3a.R`: input quality checks and cleaning.
- `code/03_variable_derivation_v3a.R`: clinical and analysis-variable derivation.
- `code/04_analysis_v3a.R`: primary analyses and model diagnostics.
- `code/05_normalization_analysis_v3a.R`: normalization discordance and exploratory
  prognostic analyses.
- `utils/`: configuration, clinical helpers, diagnostic QMDs, and the separate
  reporting utility.
- `config/config_template.json`: site configuration template.
- `outlier-thresholds/`: shared threshold files.
- `reference-data/growth/`: growth-reference files.
- `DESCRIPTION`: direct-dependency manifest; this is not an installable R package.
- `renv.lock`: recorded R package environment.

Preserve filename capitalization when moving to case-sensitive systems.

## Dependencies

Use R 4.6.1 and the recorded `renv.lock` environment. Opening the `.Rproj`
sources `renv/activate.R` through the root `.Rprofile`.

On Windows, use a short checkout path, for example
`C:/projects/PBW_vs_PFVC_peds_CLIF`. A long-path restoration test produced an
inaccessible 263-character package-header path; the short-path test restored
successfully.

`DESCRIPTION` declares 36 direct dependencies. `algorithmDiagnostics` is pinned
to commit `1cd593b613fd0aa68c937eedad97767ee1ea6f11`.

From the repository root, restore dependencies before running the pipeline:

```r
renv::restore(prompt = TRUE)
```

Restoration downloads packages and writes to the project library and package
cache. It is separate from analysis execution. Keep the lockfile unchanged;
do not substitute a package update or new snapshot for restoration.

This project supports Windows execution. The restoration test used an
existing R installation and restoration tools with a separate package
library and cache.

PDF reporting requires the Quarto command-line executable and a working TeX
installation in addition to the R package `quarto`. The tested environment
used Quarto 1.10.19 and TinyTeX. The runner does not install software.

## Site configuration

The development site's inputs are reported to comply with CLIF version 2.1,
with datetime strings expressed in UTC. Each participating site must confirm
its own schema, units, clinical mappings, and timestamp convention before
execution; filename compatibility alone does not establish those properties.

The patient-assessment input uses `numeric_value` for numeric assessments,
including GCS total. The pipeline does not assume a `numerical_value` column.

Copy `config/config_template.json` to `config/config.json`, then edit the copy:

```json
{
  "site_name": "YOUR_SITE",
  "tables_path": "/path/to/clif_tables",
  "file_type": "parquet",
  "output_root": "output",
  "supplemental_height": null
}
```

Replace the placeholder site label and input path before loading configuration.

- `site_name`: a non-placeholder label containing letters, digits, underscores,
  or hyphens.
- `tables_path`: an existing directory containing the site's CLIF input tables.
- `file_type`: exactly `parquet`, `csv`, or `fst`.
- `output_root`: the destination root for local artifacts. It must be separate
  from, and not nested within or around, `tables_path`.
- `supplemental_height`: `null`, or an object containing exactly `format` and
  `path`.

Relative paths resolve against the repository root. Loading configuration
validates settings but does not create output directories or read clinical
tables.

The template's `_comment` and `_supplemental_height_example` fields are
documentation only. They do not activate supplemental input.

An alternative configuration can be selected with `PBWPFVC_CONFIG_PATH` or the
runner's `config_path` argument. `PBWPFVC_REPO_ROOT` or the `repo_root` argument
can identify the repository explicitly.

The active `config/config.json` is excluded from version control.

### Accepted input filenames

Script 01 expects files named `clif_<table>.<file_type>` in `tables_path`.

The following table names have explicit filename alternatives:

- Respiratory support: `clif_resp_support` or `clif_respiratory_support`.
- Continuous medication administration: `clif_medication` or
  `clif_medication_admin_continuous`.
- Patient assessments: `clif_patient_assessment` or `clif_patient_assessments`.

Append the configured extension: `.parquet`, `.csv`, or `.fst`.

Exactly one alternative must exist for each of these tables. If both exist,
the loader stops rather than selecting one silently. It also stops if an
expected filename identifies a directory.

Other required table stems are `patient`, `hospitalization`, `adt`, `vitals`,
`labs`, `patient_procedures`, and `hospital_diagnosis`, each with the `clif_`
prefix and configured extension.

The diagnosis input requires:

- `hospitalization_id`
- `diagnosis_code`
- `diagnosis_code_format`
- `poa_present`

Only index-hospitalization rows with `poa_present == 1` and
`diagnosis_code_format == "ICD10CM"` contribute to classification. All diagnosis
positions are eligible. Missing code-system labels are not inferred from code
text. See the [CCC policy](reference-data/ccc_adjustment_policy.md) for the
indicator missingness and fixed adjustment rules.

The medication input must contain continuous administrations only.
Filename aliases do not rename columns, convert units, standardize clinical
categories, or prove schema compliance.

### Optional supplemental height

For a standard supplemental CSV, use:

```json
"supplemental_height": {
  "format": "height_csv",
  "path": "/path/to/height.csv"
}
```

Required CSV columns:

- `encounter_id`
- `patient_id`
- `height_cm`
- `height_dttm`

The other supported format is `michigan_crosswalk`, a Michigan encounter
crosswalk stored as Parquet. Its format is independent of the CLIF `file_type`.

A null setting disables supplemental input. An unavailable configured optional
file is reported, and cohort identification continues using CLIF vitals.

Encounter linkage and patient-conflict checks are retained in local diagnostics.
Patient-only unlinked records are not used as height donors by this reader.

Cross-encounter donor eligibility uses hospitalization admission timestamps
and the existing age-dependent windows. Supplemental measurement dates remain
unparsed source text for provenance and do not determine donor eligibility.
Do not invent measurement times.

### Supplemental source representation and precision

The development-site supplemental source is described as individual height
measurements in centimeters. `height_dttm` is the measurement timestamp,
formatted as `YYYY-MM-DD HH24:MI:SS`, in Eastern local time with daylight-saving
adjustment. No conversion, rounding, truncation, or reformatting before export
was reported. Each site should confirm its source contract.

Measurement dates remain unparsed provenance. Donor windows use hospitalization
admission timestamps. The supplemental Eastern-time convention is separate
from the CLIF UTC-export convention.

CSV tests passed for date-only, fractional-second, whitespace-containing,
invalid-date, and missing representations. Parquet tests passed for character
dates, typed timestamps, and an absent date column. Tests also covered row order,
numeric-height handling, linkage exclusions, and missing-date eligibility.

`height_text` is cleaned text. For Parquet, `height_dttm_raw` is a character
representation of the decoded column, not the original lexical field.
Future use of measurement timestamps for matching would require an explicit
timezone and daylight-saving transition policy.

## Execution

Open a fresh RStudio session through the repository `.Rproj`. Configure the site
and restore dependencies before execution.

### Full workflow

Run this block once from the repository root:

```r
source("code/00_run_pipeline_v3.R")
analysis_result <- NULL
reporting_result <- NULL

local({
  analysis_result <<- run_pbwpfvc_analysis()
  reporting_result <<- run_pbwpfvc_reports()
})
```

Sourcing the runner defines functions only. The block runs Scripts 01-05,
then diagnostic reporting. An analysis error prevents reporting from starting.

Scripts run sequentially in separate R processes using the project and base R
libraries, with system and user profiles disabled. The parent calls `gc()` after
each successful script return.

### Separate analysis or reporting

```r
source("code/00_run_pipeline_v3.R")
analysis_result <- run_pbwpfvc_analysis()
```

Once saved analysis prerequisites are complete, render reports without
rerunning analysis:

```r
reporting_result <- run_pbwpfvc_reports()
```

Reporting reads saved diagnostics and repository QMDs without refitting models.
It produces cohort, analysis, and normalization PDFs. The renderer passes the
selected libraries to nested R processes and disables their startup profiles
to avoid unintended project bootstrapping.

### Explicit script ranges

The default analysis call runs Scripts 01–05. An explicit range can reuse
completed prerequisites:

```r
analysis_result <- run_pbwpfvc_analysis(
  first_script = 4L,
  last_script = 5L
)
```

Use a partial range only when its saved prerequisite inputs are current.
The runner stops at the first failure and retains per-script warnings.

### Optional development library override

Normal execution uses this repository's project library plus R's base library.
For development, an existing compatible project library can be selected with
`PBWPFVC_PROJECT_LIBRARY`. This does not restore, install, or update packages.

```r
local({
  previous_library <- Sys.getenv("PBWPFVC_PROJECT_LIBRARY", unset = NA_character_)

  on.exit({
    if (is.na(previous_library)) {
      Sys.unsetenv("PBWPFVC_PROJECT_LIBRARY")
    } else {
      Sys.setenv(PBWPFVC_PROJECT_LIBRARY = previous_library)
    }
  }, add = TRUE)

  Sys.setenv(PBWPFVC_PROJECT_LIBRARY = "C:/projects/existing_project_library")

  analysis_result <<- run_pbwpfvc_analysis()
})
```

Replace the example with an existing R-version/platform-compatible project
library. The override is a development option, not proof of dependency-version
compatibility. The caller still needs `renv` and `callr` available.

### Warnings

Each analysis invocation retains an `analysis_warning_run_*` directory under
the model-diagnostic directory. Per-script CSVs contain the warning sequence,
UTC timestamp, script, condition class, message, and call. Scripts with no
captured warnings retain a header-only CSV.

Warnings remain visible. Synthetic tests passed for warnings before and after
global-workspace clearance, zero-warning scripts, and logging failures.
If warning retention fails, the child stops and later scripts do not run.

Warnings suppressed inside script code may not reach the handler. A hard crash
can interrupt a write. Previous child warnings cannot be recovered through
parent-session `warnings()`. This capture covers analysis, not report rendering.
Warning logs may contain paths or data values and are not part of the analysis
output manifest.

### Failure and recovery

The runner stops at the first failure and does not retry automatically.
Earlier and partial outputs remain; analysis scripts may overwrite existing
artifacts. Retain the error and last checkpoint before choosing recovery.
Do not rerun analysis just to inspect saved status.

Each reporting invocation retains a unique directory under `final/reports`,
including staged sources, render files, results, and backups. Existing final
PDFs are backed up before replacement. The utility performs no explicit deletion.

PDFs describe the pre-render diagnostic snapshot. The normalization PDF may
show reporting pending while the subsequently saved RDS records completion.

Publication, status updates, and model-QC synchronization are separate writes.
A failure can leave partial publication or disagreement between artifacts.
Inspect retained files and backups before retrying.

The runner does not install packages, pool sites, or use Git.

## Local output layout

Under `<output_root>/<site_name>_output/`:

```text
intermediate/
diagnostics/
  cohort/
  table1/
  models/
final/
  tables/
  figures/
  reports/
  optional_diagnostics/
```

Detailed diagnostic RDS files remain in the local diagnostic directories.
The workflow does not automatically copy them to an export bundle.

Script 05 validates an analysis-only inventory of 19 files:

- 14 CSV artifacts.
- 5 analysis PDFs; the two-page Part 1 PDF counts as one file.

The separate renderer preserves that inventory and validates a distinct
normalization report-inclusive inventory of 20 files: 14 CSVs and 6 PDFs.
This is not a manifest of every file produced by the entire pipeline or retained
reporting run directories.

Reporting completion is recorded after PDF and manifest validation, then the
normalization diagnostic section is synchronized into the latest model-QC RDS
while preserving unrelated sections.

## Analysis definitions

The cohort includes ages 3 to under 18 and the latest eligible hospitalization
per patient. Scripts and diagnostics record thresholds, derivations, model
specifications, and analysis-specific sample policies.

Key interpretation points:

- PBW/PFVC is a unit-bearing conversion factor, not measured lung volume.
- Terciles are cohort-relative categories, not clinical injury thresholds.
- PIP-derived fallback driving pressure differs from plateau-derived static
  driving pressure.
- AIC and apparent in-sample discrimination are exploratory comparisons.
- Odds ratios are not risk ratios.
- Fine-Gray subdistribution hazard ratios are not changes in ventilator-free days.
- Events-per-parameter criteria guide interpretation or curve display; they do
  not apply blanket cell suppression.

Sites must verify equivalent respiratory mappings, including the prescribed-target
interpretation of set PIP.

## Data and outputs

Intermediate datasets, diagnostic RDS files, and some report tables can contain
patient or hospitalization identifiers. Figures may include individual points
or small groups. Keep these artifacts local and review specific files before
sharing.

Do not commit site data, active configuration, generated outputs, or diagnostic
RDS files to GitHub. There is no automatic results-upload or cross-site pooling
workflow.

## Next steps

- Test the generic Windows workflow against another site's CLIF exports.
- Implement growth-method and interpolation comparisons as a separate script
  after sites supply their data. This script is planned, not part of the
  current pipeline.

The current growth-reference CSVs match the official CDC downloads.
Reference details are recorded in
[Growth references](reference-data/growth/README.md).

Study dates, participating sites, sample size, waivers, and access, retention,
and sharing policies are not specified by this repository.

