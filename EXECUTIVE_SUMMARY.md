# Pediatric CLIF PBW versus PFVC — Executive Summary

## Study question

This retrospective pediatric study compares predicted body weight (PBW) and
predicted forced vital capacity (PFVC) normalization across demographic groups
and examines associations of normalized ventilation and respiratory-mechanic
metrics with clinical outcomes. PBW and PFVC are formula-derived predictions,
not measured lung volumes.

## Implementation status

The full single-site analysis and three-report workflow completed in one Windows
RStudio session. Final PDF publication, the normalization report-inclusive
manifest, and reporting-status synchronization passed readback checks. Report
layouts were visually reviewed.

A short-path Windows restoration installed 186 packages; eight lockfile packages
were supplied by R's system library. All 194 versions and the pinned GitHub
commit matched, and all 36 declared dependencies loaded. The test used an
existing R installation and restoration tools. This project supports Windows
execution.

The runner uses separate child processes and calls `gc()` between analysis
scripts. The full pipeline completed successfully after this change, which
is retained in the standard workflow.

GitHub contains shared code, configuration templates, references, and the
lockfile. Site data and generated results remain local.

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

Common complete-case samples are retained within mechanic families.
Model comparisons report AIC, apparent in-sample discrimination, event counts,
and computational diagnostics.

## Execution and dependencies

Use a fresh session opened through the repository `.Rproj`. The README provides
the full execution block: load `code/00_run_pipeline_v3.R`, run analysis, then
render reports. Loading definitions starts neither workflow; an analysis error
prevents reporting.

Scripts run sequentially in separate R processes. The parent calls `gc()` after
each successful return. The runner stops on failure, retains existing or partial
outputs, and does not retry automatically.

Restore the recorded `renv.lock` environment before execution. On Windows,
use a short checkout path to avoid package-header path problems. Restoration
downloads and installs packages; the analysis runner does not.

Analysis warnings are retained per script in a unique local directory and remain
visible. Logging failure stops the child. Suppressed warnings and writes
interrupted by hard crashes may not be retained.

Quarto and TeX render the PDFs from saved diagnostics without refitting models.
Reports show the pre-render snapshot. Each run retains staged files, results,
and backups. Publication and diagnostic synchronization are separate writes;
inspect retained evidence before recovering from a failure.

## Outputs

- `intermediate/`: working and analysis datasets.
- `diagnostics/cohort`, `diagnostics/table1`, and `diagnostics/models`: saved
  diagnostics; analysis-warning runs are under the model-diagnostic directory.
- `final/tables`, `final/figures`, and `final/reports`: analysis artifacts and PDFs.
- `final/optional_diagnostics`: reserved destination; no automatic RDS export.

Script 05 validates 19 analysis files: 14 CSVs and five PDFs. Reporting validates
a 20-file normalization report-inclusive manifest: 14 CSVs and six PDFs.
These are not whole-pipeline inventories.

Some datasets, diagnostics, and report tables contain identifiers; figures may
show individual points or small groups. Selected minimum-group filters are not
blanket cell suppression. Keep site data, active configuration, generated
outputs, and diagnostic RDS files out of GitHub.

Cross-site pooling, control cohorts, and the adult study's longitudinal
injury-stage workflows are not implemented.

## Next steps

The next generic test is another-site execution on Windows.

Growth-method and interpolation comparisons will be implemented as a separate
script after sites supply their data. This script is planned and is not part
of the current pipeline. The current reference CSVs match the official CDC
downloads; see [Growth references](reference-data/growth/README.md).

Study dates, participating sites, sample size, waivers, and access, retention,
and sharing policies are not specified here.

