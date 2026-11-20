# Custom respiratory-mechanics diagnosis groups scaffolded on the official CCC v3 table.
# This is not the unchanged published CCC classification or its technology algorithm.
# Diagnosis-only mapping; procedure codes do not contribute.
# Asthma and reassignment of selected scoliosis families are deliberate study adaptations.
#
# Source: Children's Hospital Association, Appendix Table: Categories of CCC v3
# and the Corresponding ICD-9 and ICD-10 Diagnosis and Procedure Codes.
# No local source download date is asserted.
#
# Entries denote the listed code and its descendants. Matching uses a separate
# normalized copy of each diagnosis code; source diagnosis text remains unchanged.

resp_mechanics_ccc_mapping <- function() {
  source_url <- paste0(
    "https://www.childrenshospitals.org/p/_bBqgXQiQby8yFkH2ztwpg/",
    "ep2Cgl9GSHieIn_Ym1eoOw-final-table-v3-with-rev-ccc-list-02282023.pdf"
  )
  
  prefixes <- list(
    ccc_neuro_mechanics = c(
      "G80",
      "G71", "G72",
      "G10", "G20", "G21", "G23", "G24.02", "G24.1", "G24.2", "G24.8",
      "G25.5", "G25.8", "G25.9", "G26",
      "M41.0", "M41.4", "M41.5"
    ),
    ccc_neuro_other = c(
      "G90.1", "Q00", "Q01", "Q02", "Q03", "Q04", "Q05", "Q06", "Q07",
      "F71", "F72", "F73",
      "F84.2", "G11", "G12", "G31.0", "G31.8", "G31.9", "G32.89", "G90.9",
      "G91.0", "G91.1", "G91.3", "G91.4", "G91.8", "G91.9",
      "G93.8", "G93.9", "G94", "G95.19", "G95.8", "Q85.1",
      "G40.011", "G40.019", "G40.111", "G40.119", "G40.211", "G40.219",
      "G40.311", "G40.319", "G40.411", "G40.419", "G40.42",
      "G40.803", "G40.804", "G40.813", "G40.814", "G40.823", "G40.824",
      "G40.833", "G40.834", "G40.911", "G40.919",
      "G40.A11", "G40.A19", "G40.B11", "G40.B19",
      "G36.1", "G36.8", "G36.9", "G37", "G81", "G82.5",
      "G83.5", "G83.9", "G93.1", "G93.5", "R40.3",
      "I63"
    ),
    ccc_resp_chronic_malformation = c(
      "Q30", "Q31", "Q32", "Q33", "Q34",
      "G47.35", "I27.0", "I27.82", "J84.1", "J84.2", "J84.8", "J84.9",
      "J96.1", "J96.2", "Z90.2"
    ),
    ccc_resp_other = c(
      "E84", "T86.3", "T86.81", "Z94.2", "Z94.3"
    ),
    ccc_asthma = "J45",
    ccc_resp_technology = c(
      "J95.0", "J95.850", "Z43.0", "Z93.0", "Z99.0", "Z99.11", "Z99.12"
    )
  )
  
  labels <- c(
    ccc_neuro_mechanics = "Mechanics-relevant neurologic/neuromuscular disease",
    ccc_neuro_other = "Other neurologic disease",
    ccc_resp_chronic_malformation = "Chronic respiratory disease/malformations",
    ccc_resp_other = "Other respiratory disease",
    ccc_asthma = "Asthma",
    ccc_resp_technology = "POA respiratory technology evidence"
  )
  
  list(
    mapping_version = "custom_resp_mechanics_v3a_1",
    classification = "Custom diagnosis groups scaffolded on the CCC v3 appendix table",
    source_url = source_url,
    code_system = "ICD-10-CM",
    matching = "Listed codes and descendants; uppercase and remove decimal points for matching",
    prefixes = prefixes,
    labels = labels,
    adaptations = c(
      "Asthma added as a separate J45 group",
      "Selected scoliosis families reassigned to the mechanics-relevant neurologic group",
      "Chronic respiratory diseases and respiratory malformations combined",
      "Remaining respiratory diagnosis subcategories combined",
      "Respiratory technology evidence independent of other condition groups",
      "Diagnosis codes only; published procedure-code contributions omitted"
    ),
    limitations = c(
      "Published Other respiratory subcategory has no diagnosis-code contribution",
      "Not an implementation of the unchanged CCC v3 algorithm",
      "No claim that technology evidence specifically identifies home invasive ventilation",
      "No annual ICD-10-CM validity crosswalk is applied by this mapping"
    )
  )
}