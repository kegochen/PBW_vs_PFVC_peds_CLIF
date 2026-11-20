# CCC adjustment policy

The pediatric PBW-versus-PFVC pipeline uses study-specific diagnosis groups
scaffolded on pediatric Complex Chronic Conditions version 3. Asthma is added
separately. This is not an unchanged implementation of the published CCC
classification.

## Classification

Classification uses diagnoses from the selected index hospitalization only.
A diagnosis contributes when `poa_present == 1` and
`diagnosis_code_format == "ICD10CM"`. All diagnosis positions are eligible.
Procedure codes and prior hospitalizations do not contribute.

The six indicators are:

- `ccc_neuro_mechanics`: cerebral palsy, muscular dystrophies and myopathies,
  selected movement disorders, and selected scoliosis families.
- `ccc_neuro_other`: remaining neurologic diagnoses in the selected framework.
- `ccc_resp_chronic_malformation`: chronic respiratory diseases and respiratory
  malformations.
- `ccc_resp_other`: cystic fibrosis and respiratory transplant-related diagnoses.
- `ccc_asthma`: asthma.
- `ccc_resp_technology`: respiratory technology status, care, or complication
  evidence. This does not specifically establish home invasive ventilation.

Groups may overlap. The mapping utility contains the code definitions and
reference provenance.

Observed qualifying evidence gives an indicator value of 1, including when
other diagnosis rows have coding gaps. Absence is coded 0 only when usable
present-on-admission diagnosis information is available without unresolved
coding gaps. Otherwise the indicator is missing.

Code-format and structural checks are not validation against an annual
ICD-10-CM code dictionary.

## Fixed model-family term sets

Use the same term sets at each site:

- Mortality: neurologic mechanics, other neurologic, chronic respiratory
  disease/malformations, and asthma.
- Non-mechanical-power continuous outcomes: all indicators except other
  respiratory disease.
- Mechanical-power continuous outcomes: all indicators except other respiratory
  disease and respiratory technology.
- Fine–Gray liberation: all indicators except other respiratory disease.
- Cox mortality: all indicators except respiratory technology.
- Prediction-formula regressions and Script 05 regressions: neurologic mechanics,
  other neurologic, chronic respiratory disease/malformations, and asthma.

The omissions originated from zero prevalence, singleton leverage, or
sparse-event estimation problems during development. They are fixed study
specifications, not site-specific automatic selection rules.

Primary and minimally adjusted models share the applicable condition
adjustment. Their existing demographic-adjustment distinctions remain intact.

## Samples and diagnostics

All six indicators remain in the applicable regression completeness
requirements, including indicators omitted from a formula. Missing indicators
are not recoded as absence. Descriptive cohorts and cohort eligibility are
unchanged.

Retain the existing common-sample comparison policies, outcome definitions,
inference, and model diagnostics. Do not automatically remove a term or patient,
change estimators, or suppress warnings to resolve a local estimation failure.
Review actual compatibility or estimation findings before making a separately
agreed change.

Table 1 describes all six groups, including those omitted from some models.
Part 1 adjusted prediction curves use the stated sex/race references with the
four regression adjustment groups absent. This profile does not exclude
condition-positive patients from model fitting.

## Exploratory and deferred analyses

Individual-modifier interactions are exploratory. Their related exposure
specifications use unadjusted P values.

Potential restoration of omitted terms and co-occurrence interactions in
pooled cohorts is described in `ccc_pooled_candidates.md`. These candidates
are not enabled for automatic execution and remain deferred until buddy-site
testing is complete.