# Pediatric CLIF: PBW versus PFVC

Single-site R analysis and diagnostic reporting for a pediatric Common
Longitudinal ICU data Format (CLIF) study comparing predicted body weight
(PBW) and predicted forced vital capacity (PFVC) normalization of respiratory
physiology measures.

For a high-level description of the pediatric cohort, analysis workflow,
outputs, and interpretation limits, see the
[Executive summary](EXECUTIVE_SUMMARY.md). It describes the implementation,
locally verified execution status, and remaining validation limits.
It is not an approved study protocol.

## Current status

This repository supports local single-site execution with shared code
distributed through GitHub.

A complete single-session analysis-and-reporting invocation returned
successfully in the local Windows development environment. Scripts 01–05
completed sequentially, followed by all three diagnostic reports. The renderer
validated staged PDFs, published the final PDFs, validated the normalization
report-inclusive manifest, and read back and synchronized reporting status
into model QC. All three diagnostic-report layouts were previously visually
accepted.

The successful invocation used parent-session garbage collection after each
analysis subprocess returned. Earlier attempts experienced intermittent
subprocess failures. One successful invocation does not establish that
garbage collection caused the improvement or guarantee future reliability.

Dependency loading and lockfile coverage have been checked locally.
Cross-platform dependency restoration remains unvalidated.

Do not interpret code availability, process completion, or computational checks
as scientific validation, disclosure approval, or authorization to execute.

GitHub hosts shared code, configuration templates, reference files, and the
dependency lockfile—not site data or analysis results. Cross-site pooling is
not implemented.

## Repository structure

- `code/00_run_pipeline_v3.R`: function-based local workflow runner.
- `code/01_cohort_identification_v3.R`: cohort identification.
- `code/02_quality_checks_v3.R`: input quality checks and cleaning.
- `code/03_variable_derivation_v3.R`: clinical and analysis-variable derivation.
- `code/04_analysis_v3.R`: primary analyses and model diagnostics.
- `code/05_normalization_analysis_v3.R`: normalization discordance and exploratory
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

The recorded environment uses R 4.6.1 and `renv`. The root `.Rprofile` sources
`renv/activate.R`, so opening the project can activate its package environment.

`DESCRIPTION` declares 36 direct dependencies. `callr` provides isolated
subprocess execution. `algorithmDiagnostics` is pinned to GitHub commit:

`1cd593b613fd0aa68c937eedad97767ee1ea6f11`

A collaborator must restore the project environment before running the workflow.
Restoration can download packages and write to the project library and renv
cache. It is a separate setup action, not part of the runner.

After approving that setup action, run from the repository root:

```r
renv::restore(prompt = TRUE)
```

Do not substitute an unrestricted update or snapshot for restoration of the
recorded environment.

PDF reporting additionally requires the external Quarto command-line executable
and a working TeX installation. Installing the R package `quarto` alone does not
provide both external prerequisites. The reporting utility does not install
software.

The development machine has Quarto 1.10.19 and TinyTeX available.
All three diagnostic PDFs have rendered successfully in that environment.
This does not establish successful rendering on another system.

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
`labs`, and `patient_procedures`, each with the `clif_` prefix and configured
extension.

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

## Execution: analysis and reporting are separate

The examples below are for separately authorized local execution.
Loading definitions, running analysis, and rendering reports are distinct
actions. Do not rerun a completed workflow merely to inspect its status.

From the repository root, load the runner definitions:

```r
source("code/00_run_pipeline_v3.R")
```

Sourcing this file defines functions only. It does not start analysis or reporting.

### Standard end-to-end invocation

Use a fresh RStudio session opened through the repository `.Rproj`, with the
validated configuration and existing project dependencies. After local
analysis and reporting execution are authorized, run this block once:

```r
source("code/00_run_pipeline_v3.R")
analysis_result <- NULL
reporting_result <- NULL
local({
  analysis_result <<- run_pbwpfvc_analysis()
  reporting_result <<- run_pbwpfvc_reports()
})
```

An analysis error stops the block before reporting starts. Analysis and
reporting remain separate functions, but the standard block invokes them
in order within one parent session. A three-session split is not required
by this standard.

Existing analysis artifacts may be overwritten. Reporting retains its own
run directory and backups. Do not automatically repeat a failed invocation;
retain the error and inspect affected artifacts before choosing recovery.

### Analysis

After execution is authorized:

```r
run_pbwpfvc_analysis()
```

The runner executes Scripts 01–05 sequentially in separate R processes using
the project library and base R library. User and system profiles are disabled
in those child processes. After each successful subprocess return, the parent
R session calls `gc()` before continuing, including before reporting after
Script 05. This requests collection of unused parent-session R memory; it
does not remove live analysis objects, guarantee memory release to the
operating system, or establish the cause of earlier subprocess failures.

The runner stops at the first error. Earlier and partial outputs remain;
existing analysis artifacts may have been overwritten by the scripts. It does
not automatically retry, delete outputs, render diagnostic reports, install
packages, pool sites, or use Git.

Successful subprocess return means process completion, not independent
scientific validation.

### Diagnostic PDF reporting

After analysis prerequisites are complete and reporting is separately authorized:

```r
run_pbwpfvc_reports()
```

Reporting consumes saved diagnostic RDS files and repository QMDs without
refitting models. It produces cohort, analysis, and normalization diagnostic PDFs.

Each invocation retains a unique run directory under `final/reports`, including
render artifacts and backups. Existing final PDFs are backed up before
replacement. The utility performs no explicit file deletion.

The PDFs describe the diagnostic snapshot read before rendering. In particular,
a normalization PDF can show reporting pending while the subsequently saved
diagnostic RDS records successful rendering and validation.

Rendering, publication, diagnostic status updates, and model-QC synchronization
are multiple operations, not an atomic transaction. An error can leave partial
outputs or partially updated diagnostic artifacts. Inspect the reported error
and retained run files before considering a retry.

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

## Scientific interpretation

The cohort targets ages 3 to under 18 and retains the latest eligible
hospitalization per patient. Clinical thresholds, derivation rules, model
specifications, and sample policies are retained in the scripts and diagnostics.

Important interpretation limits include:

- PBW/PFVC is a unit-bearing conversion factor, not measured lung volume.
- Cohort-relative terciles are not validated clinical injury thresholds.
- PIP-derived driving-pressure fallbacks are not equivalent to plateau-derived
  static driving pressure.
- Apparent in-sample discrimination and AIC comparisons are exploratory.
- Odds ratios are not risk ratios.
- Fine–Gray subdistribution hazard ratios are not changes in ventilator-free days.
- Computational stability does not establish causal identification, clinical
  validity, or a preferred normalization policy.
- Events-per-parameter criteria are interpretation or visualization heuristics,
  not blanket disclosure rules or absolute validity boundaries.

Sites must verify that their respiratory mappings represent equivalent clinical
concepts, including the prescribed-target interpretation of set PIP.

## Disclosure and sharing

Outputs and detailed diagnostics require review before any sharing.
Some local diagnostics and report tables may contain patient or hospitalization
identifiers. No report or artifact bundle has been established as wholly
identifier-free or approved for disclosure.

There is no implemented cross-site pooling or automatic results-upload workflow.
Do not commit site data, active configuration, generated outputs, or diagnostic
RDS files to GitHub.

A prospective study or IRB summary describes intended behavior; it does not
establish protocol approval, completed implementation, or privacy or security
assurances.

## Outstanding validation

Successful local execution and reporting do not establish portability,
independent scientific validation, or disclosure clearance.
Remaining reviews include:

- Participating-site CLIF schema, units, clinical mappings, and timestamp compatibility.
- Supplemental measurement-date provenance and source precision.
- Full growth-reference comparison and provenance.
- Diagnostic preservation across future reruns.
- Reporting failure paths and recovery behavior.
- Cross-platform dependency restoration.
- Reliability across repeated end-to-end invocations and other environments.
- Review of analysis warning messages not captured in earlier runs.
- Independent scientific validation and artifact-specific disclosure review.

Study dates, participating sites, sample size, waivers, and access, retention,
and sharing policies are not specified by this repository.
