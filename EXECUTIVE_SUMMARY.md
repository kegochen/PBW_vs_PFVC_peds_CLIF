# Pediatric CLIF PBW versus PFVC — Executive Summary

## Study question

This retrospective pediatric study examines how predicted body weight (PBW)
and predicted forced vital capacity (PFVC) normalization differ across
demographic groups and how the resulting ventilation and respiratory-mechanic
metrics associate with clinical outcomes.

PBW and PFVC are predicted quantities. Their relationship does not establish
measured lung volume, prediction accuracy, causal injury, or the clinical
superiority of either normalization method.

## Implementation status

The repository has been refactored for local, single-site execution with shared
code suitable for distribution through GitHub.

A complete single-session invocation of Scripts 01–05 followed by diagnostic
reporting returned successfully in the local Windows development environment.
The renderer validated and published all three PDFs, validated the
normalization report-inclusive manifest, and read back and synchronized
reporting status into model QC. All three diagnostic-report layouts were
previously visually accepted.

The successful invocation used parent-session garbage collection between
analysis subprocesses. Earlier attempts experienced intermittent subprocess
failures; their cause remains unproved. This success does not establish a
causal garbage-collection fix or reliability across future invocations.

Dependency loading has been checked locally; cross-platform restoration
remains unvalidated. Process completion, file integrity, and visual acceptance
do not establish independent scientific validity or disclosure approval.

The workflow below describes the implementation and its remaining limits.
GitHub hosts code and reference material, not site results.

## Cohort and scientific definitions

The cohort includes children aged 3 to under 18 and selects the latest eligible
hospitalization per patient.

Respiratory intubation anchors are identified independently of procedure-code
31500 evidence. Eligibility requires a positive set tidal volume, or positive
observed fallback, during confirmed invasive support within one hour of a
qualifying respiratory anchor, together with ICU eligibility.

Subsequent derivation requires an eligible complete invasive timepoint meeting
at least one hypoxemia criterion: OI >= 4, OSI >= 5, PF ratio <= 300, or
SF ratio <= 250. VT/PBW is not restricted to a 6–8 mL/kg eligibility band.

The index-selection rule prefers the first qualifying timepoint with valid
driving pressure within the existing six-hour window; otherwise it selects
the first qualifying timepoint.

## What the pipeline does

### Script 01 — Cohort identification

Loads the required CLIF tables, reconciles identifiers, identifies respiratory
anchors and one-hour tidal-volume eligibility, selects index hospitalizations,
processes respiratory support, extracts clinical inputs, and starts the
hospitalization-level attrition log.

Height uses index-hospitalization vitals first, optional linked supplemental
height second, then an eligible cross-encounter donor. Donor eligibility uses
hospitalization admission timestamps and the existing age-dependent windows.
Supplemental measurement dates remain unparsed provenance text.

### Script 02 — Quality checks

Applies retained laboratory, vital-sign, height, and respiratory-support
thresholds. Preserves detailed availability, deduplication, provenance,
and cleaning diagnostics.

### Script 03 — Variable derivation

Calculates McLaren PBW, NCHS predicted weight, GLI-based PFVC, actual-weight
comparisons, normalized tidal volumes, respiratory mechanics, oxygenation,
pediatric SOFA, survival variables, and ventilator-free-day outcomes.

It produces the cross-sectional and retained-timepoint datasets, completes
the six-step hospitalization-level attrition log, and writes distribution
summaries using fixed shared bins.

Selected tidal volumes are cleaned without substituting another source after
an invalid selected value becomes missing. PIP-derived driving-pressure
fallbacks remain distinguished from plateau-derived static driving pressure.

### Script 04 — Analysis

Produces cohort characteristics, mortality and respiratory-mechanic models,
Fine–Gray liberation analyses with death as a competing event, survival
analyses, demographic associations, PFVC-versus-weight comparisons, model
diagnostics, figures, and unified effect-estimate tables.

Models retain their existing adjustment specifications and analysis-specific
complete-case policies. A subdistribution hazard ratio is not a change in
ventilator-free days, and an odds ratio is not a risk ratio.

### Script 05 — Normalization analysis

Examines normalization discordance, demographic and age patterns, and movement
between independently constructed, tie-preserving cohort-relative terciles.

Exploratory prognostic analyses compare mechanic-only, PBW-locked, PFVC-locked,
and separate mechanic-plus-size models across seven alternative adjustment
terms. Additional analyses evaluate model-form comparisons, demographic
modification, age/PFVC interaction ladders, balanced modifier models, and
encompassing size models.

Common samples are retained within the specified mechanic families. Apparent
in-sample discrimination and AIC comparisons do not establish external
validation or a preferred clinical normalization policy.

## Execution and dependencies

`code/00_run_pipeline_v3.R` defines an analysis runner and a separate reporting
runner. Loading that file does not start either workflow.

The standard workflow uses a fresh repository `.Rproj` session and invokes
analysis followed by reporting in one parent session. See the README for the
complete execution block; an analysis error prevents reporting from starting.

Analysis scripts execute sequentially in isolated R subprocesses. The parent
session calls `gc()` after each successful return. The runner stops at the
first failure and retains existing or partial artifacts; it does not retry
automatically. Garbage collection does not remove live objects or establish
the cause of earlier failures.

Dependency restoration is a separate, explicitly approved setup action using
`renv`; the runner does not install packages.

Diagnostic PDFs are rendered separately from saved RDS diagnostics without
refitting models. Reports describe the pre-render diagnostic snapshot.
Reporting completion is recorded afterward following file and manifest
validation.

Reporting operations are not atomic. A failure can leave partial publication
or disagreement between diagnostic artifacts; retained backups and the
original error must be reviewed before retrying.

## Outputs and disclosure

Local outputs are organized under:

- `intermediate/`: row-level working and analysis datasets.
- `diagnostics/cohort`, `diagnostics/table1`, and `diagnostics/models`:
  detailed saved diagnostics.
- `final/tables`, `final/figures`, and `final/reports`: analysis artifacts
  and human-readable outputs.
- `final/optional_diagnostics`: reserved optional diagnostic destination;
  no automatic RDS export is implied.

Script 05 maintains an analysis-only manifest of 19 files: 14 CSVs and
5 analysis PDFs. Separate reporting validates a normalization report-inclusive
manifest of 20 files: 14 CSVs and 6 PDFs. These are not whole-pipeline inventories.

Minimum-group filters in selected distribution outputs are not blanket
cell-suppression or disclosure clearance. Some diagnostics and report tables
may contain identifiers. Neither the final directory nor the reports have
been established as wholly identifier-free or approved for sharing.

Cross-site pooling, control cohorts, and the adult study's longitudinal
injury-stage workflows are not implemented in this pediatric repository.

## Outstanding decisions and validation

Remaining work includes source precision and provenance review, full growth-file
comparison, participating-site input column/unit and timestamp compatibility,
cross-platform dependency restoration, reliability across repeated end-to-end
invocations, reporting failure-path review, review of uncaptured analysis
warnings, independent scientific validation, and artifact-specific disclosure
review.

Local analysis execution and the complete reporting workflow have succeeded;
these remaining reviews must not be represented as already completed.

Study dates, participating sites, sample size, waivers, and access, retention,
and sharing policies remain unspecified. This document does not establish
IRB approval, execution authorization, or privacy or security assurances.
