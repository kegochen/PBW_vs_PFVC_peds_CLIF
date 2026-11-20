# Pediatric CLIF PBW versus PFVC — Executive Summary

## Study question

This retrospective pediatric study compares predicted body weight (PBW) and
predicted forced vital capacity (PFVC) normalization across demographic groups
and examines associations of normalized ventilation and respiratory-mechanic
metrics with clinical outcomes. PBW and PFVC are formula-derived predictions,
not measured lung volumes.

## Implementation status

Version 3a adds study-specific chronic-condition adjustment. Version 3 is
preserved at baseline commit `26e023228ebe6a27c2cb8364c97921b051f81420`.

The updated analysis completed on Windows. Separate reporting published all
three diagnostic PDFs, validated the normalization report-inclusive manifest,
and synchronized reporting status into the saved diagnostics. All three
reports were visually reviewed.

Local version 3a execution used the existing baseline project library through
an explicit override. Separate version 3a dependency restoration has not been
tested. The completed short-path Windows restoration test remains version 3
baseline evidence.

The published repository contains shared code, templates, references, and the
lockfile. Version 3a publication is pending. Site data and generated results
remain local.

## Cohort and scientific definitions

The cohort includes children aged 3 to under 18 and selects the latest eligible
hospitalization per patient.

Respiratory intubation anchors are identified independently of procedure-code
31500 evidence. Eligibility requires a positive set tidal volume, or positive
observed fallback, during confirmed invasive support within one hour of a
qualifying respiratory anchor, together with ICU eligibility.

Subsequent derivation requires a complete eligible invasive timepoint meeting
at least one hypoxemia criterion: OI >= 4, OSI >= 5, PF ratio <= 300, or
SF ratio <= 250. VT/PBW is not restricted to a 6–8 mL/kg eligibility band.

Index selection prefers the first qualifying timepoint with valid driving
pressure within the existing six-hour window; otherwise it selects the first
qualifying timepoint.

## Chronic-condition adjustment

Six overlapping diagnosis groups are customized from pediatric Complex Chronic
Conditions version 3, with asthma added separately:

- Mechanics-relevant neurologic/neuromuscular disease.
- Other neurologic disease.
- Chronic respiratory disease/malformations.
- Other respiratory disease.
- Asthma.
- Respiratory technology evidence.

Classification uses diagnoses present on admission during the selected index
hospitalization, with an explicit ICD10CM code-system label. All diagnosis
positions are eligible. Prior hospitalizations and procedure codes do not
contribute.

Observed qualifying evidence is coded present. Absence requires usable
present-on-admission diagnosis information without unresolved coding gaps;
otherwise the indicator remains missing. Respiratory technology evidence does
not specifically establish home invasive ventilation.

Fixed model-family adjustment sets retain the terms that were estimable during
development. All six indicators remain available descriptively and in the
applicable regression completeness requirements. See
[CCC adjustment policy](reference-data/ccc_adjustment_policy.md) for definitions,
missingness, and formula omissions.

## Analysis workflow

### Script 01 — Cohort identification

Loads CLIF tables, reconciles identifiers, identifies respiratory anchors and
one-hour tidal-volume eligibility, selects index hospitalizations, processes
respiratory support, extracts clinical inputs, classifies chronic-condition
diagnoses, and starts the hospitalization-level attrition log.

Height uses index-hospitalization vitals first, optional linked supplemental
height second, then an eligible cross-encounter donor. Donor eligibility uses
hospitalization admission timestamps and existing age-dependent windows.
Supplemental measurement dates remain unparsed provenance text.

### Script 02 — Quality checks

Applies retained laboratory, vital-sign, height, and respiratory-support
thresholds. Preserves availability, deduplication, provenance, and cleaning
diagnostics.

### Script 03 — Variable derivation

Calculates McLaren PBW, NCHS predicted weight, GLI-based PFVC, actual-weight
comparisons, normalized tidal volumes, respiratory mechanics, oxygenation,
pediatric SOFA, survival variables, and ventilator-free-day outcomes.

Links the chronic-condition indicators into analysis datasets, completes the
six-step hospitalization-level attrition log, and writes distribution summaries
using shared bins. Chronic-condition availability does not alter cohort
eligibility.

Selected tidal volumes are cleaned without substituting another source after
an invalid selected value becomes missing. PIP-derived driving-pressure
fallbacks remain distinguished from plateau-derived static driving pressure.

### Script 04 — Analysis

Produces Table 1, adjusted mortality and respiratory-mechanic models,
Fine–Gray liberation analyses with death as a competing event, Cox survival
models, demographic associations, PFVC-versus-weight comparisons, diagnostics,
figures, and unified effect-estimate tables.

Primary and minimally adjusted models share the applicable chronic-condition
adjustment while retaining their demographic-adjustment distinctions.
Family-specific omissions are documented. Kaplan–Meier curves and other raw
descriptive figures remain unadjusted.

Exploratory interactions evaluate asthma, chronic respiratory
disease/malformations, and mechanics-relevant neurologic disease separately
against four tidal-volume exposures for mortality and liberation. The 24
related interaction tests use unadjusted P values. Conditional slopes and
interaction estimates do not establish treatment benefit or a preferred
normalization method.

Iteration-only paired adjustment fits and duplicate exports were removed
after the final adjustment specification was selected.

### Script 05 — Normalization analysis

Examines normalization discordance, demographic and age patterns, and movement
between independently constructed, tie-preserving cohort-relative terciles.

Part 1 regression models include four chronic-condition terms. Adjusted
prediction curves use the stated sex/race references with those conditions
absent; condition-positive patients remain eligible for model fitting.
Raw distributions, mathematical identities, and tercile assignments are
unchanged.

Part 2 mortality models compare mechanic-only, PBW-locked, PFVC-locked, and
separate mechanic-plus-size specifications across seven alternative adjustment
terms. Demographic modification, age/PFVC ladders, balanced modifier models,
and encompassing size models retain the agreed chronic-condition adjustment.

Common complete-case samples are used within comparison families. Existing
computational suppression and curve-display criteria remain in place.
Unstable mechanical-power demographic interaction tests were retained as
suppressed results rather than interpreted estimates.

## Execution and dependencies

Use R 4.6.1 and the recorded lockfile environment. The README provides setup,
full-workflow, bounded-analysis, and reporting instructions using
`code/00_run_pipeline_v3a.R`.

Sourcing the runner defines functions only. Analysis scripts run sequentially
in separate R processes using the selected project library plus base R, with
child profiles disabled. The parent calls `gc()` after each successful script
return.

The runner stops at the first failure, retains partial outputs, and does not
retry automatically. Analysis warnings remain visible and are retained in
per-script CSVs. Hard crashes can interrupt logging.

Normal setup restores dependencies separately from execution. On Windows,
use a short checkout path. An optional development-library override selects
an existing compatible library without installing or updating packages.

Quarto and TeX render diagnostic PDFs from saved diagnostics without refitting
models. Nested R startup and library routing are controlled to prevent
unintended project bootstrapping. The reporting process calls `gc()` between
validated renders.

Reporting retains staged files, results, and backups. Final publication and
diagnostic synchronization are multiple non-atomic writes. Reports describe
the pre-render diagnostic snapshot.

## Outputs and interpretation

- `intermediate/`: working and analysis datasets.
- `diagnostics/`: cohort, Table 1, model, and warning artifacts.
- `final/tables`, `final/figures`, and `final/reports`: analysis deliverables.
- `final/optional_diagnostics`: reserved destination.

Script 05 validates 19 analysis files: 14 CSVs and five PDFs. Separate reporting
validates a 20-file normalization inventory: 14 CSVs and six PDFs. These are
not whole-pipeline inventories.

PBW/PFVC is a unit-bearing conversion factor. Terciles are cohort-relative
categories, not clinical injury thresholds. AIC and apparent in-sample
discrimination remain exploratory. Odds ratios are not risk ratios, and
Fine–Gray subdistribution hazard ratios are not changes in ventilator-free days.

Datasets, diagnostics, and report tables may contain identifiers; figures may
show individual points or small groups. Keep site data, active configuration,
generated outputs, and diagnostic RDS files out of GitHub.

Cross-site pooling, control cohorts, and the adult study's longitudinal
injury-stage workflows are not implemented.

## Remaining work

Complete version 3a publication and the focused neurologic evidence review
before buddy-site testing. Buddy-site testing remains deferred.

Pooled exploratory analyses and the separate growth-method/interpolation
comparison remain deferred until buddy-site testing is complete. Retained
pooled candidates are described in
[Deferred pooled analyses](reference-data/ccc_pooled_candidates.md).

The growth-reference CSVs match the official CDC downloads; see
[Growth references](reference-data/growth/README.md).

Study dates, participating sites, sample size, waivers, and access, retention,
and sharing policies are not specified here.

