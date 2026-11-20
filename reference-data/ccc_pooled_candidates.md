# Deferred pooled CCC analyses

The study uses customized diagnosis groups scaffolded on pediatric Complex
Chronic Conditions version 3, with asthma added separately. Classification
uses present-on-admission diagnoses from the selected index hospitalization.

## Co-occurrence interactions

Retain these candidates for exploratory evaluation in pooled cohorts:

- Asthma, chronic respiratory disease/malformations, and mechanics-relevant
  neurologic disease all present.
- Asthma and chronic respiratory disease/malformations both present.
- Mechanics-relevant neurologic disease and chronic respiratory
  disease/malformations both present.

Local triple-positive information was insufficient for useful adjusted
interaction estimation. No triple-positive interaction was fitted. Pairwise
follow-ups were deferred under the planned gate, not individually established
as nonviable.

These are co-occurrence hypotheses, not a validated increasing-condition-count
score. Define the exposure, comparison group, interaction parameterization, and
criteria for further testing before fitting pooled models.

## Locally omitted adjustment terms

Retain the following terms for possible restoration when pooled data support
their estimation:

- Mortality: other respiratory disease and respiratory technology evidence.
- Non-mechanical-power continuous outcomes: other respiratory disease.
- Mechanical-power outcomes: other respiratory disease and respiratory
  technology evidence.
- Fine–Gray liberation: other respiratory disease.
- Cox mortality: respiratory technology evidence.

Local omissions addressed zero prevalence, singleton leverage with undefined
HC3 covariance, or sparse-event coefficient instability. Patients and condition
indicators were retained; these were formula omissions, not cohort exclusions.

## Pooled implementation

These candidates are not enabled for automatic execution. Before implementation,
agree on pooled model scope and assess code compatibility, condition availability,
group overlap, outcome events, exposure variation, and estimability. Larger
samples do not by themselves establish reliable estimation.

Preserve the agreed endpoints, inference, common-sample comparison policies,
and disclosure handling unless a separate change is approved. Retain these
analyses as exploratory and report the related testing scope.