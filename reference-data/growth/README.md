# Growth references

These shared reference files are used by
`code/03_variable_derivation_v3.R`:

- `statage.csv`: stature-for-age reference.
- `wtage.csv`: weight-for-age reference.

## Reported source

Kevin Gochenour identified the following CDC page as the source of both files:

[CDC Growth Charts Data Files](https://www.cdc.gov/growthcharts/cdc-data-files.htm)

The CDC documentation describes the relevant older-child references as part
of the 2000 CDC Growth Charts and states that the data remain unchanged from
the initial May 30, 2000 release.

The original download date was not supplied and remains unknown.

## Complete file comparison

On November 19, 2026, both repository CSVs were compared with newly retrieved
official CDC files:

- [Stature-for-age CSV](https://www.cdc.gov/growthcharts/data/zscore/statage.csv)
- [Weight-for-age CSV](https://www.cdc.gov/growthcharts/data/zscore/wtage.csv)

Both comparisons passed: each file contained 436 rows, and byte-level MD5
checksums, column names, dimensions, row order, and all parsed column values
were identical. The comparison used the pipeline's existing base R
`read.csv()` parsing behavior and no numerical tolerance.

The downloaded files, retrieval timestamps, checksums, and complete comparison
results were retained in a local audit directory outside the repository.
Repository reference files were not replaced or modified.

This establishes equivalence to the CDC files retrieved for that comparison,
not the original download date, scientific suitability of the prediction
methods, or an independently approved interpolation policy.

## Fields used by this pipeline

Script 03 requires:

- `Sex`: sex code used to select the reference subset.
- `Agemos`: reference age in months.
- `M`: median reference measurement.

Both files also include LMS parameters and percentile columns. The complete
comparison confirmed that `M` equals `P50` throughout each repository file
and its corresponding CDC download.

## Scientific use

The existing McLaren PBW calculation:

1. Interpolates the child's height against median stature to obtain a
   corresponding reference age.
2. Interpolates that age against median weight.

The existing NCHS-labelled predicted-weight calculation interpolates
chronological age in months against the median weight-for-age reference.

Both calculations retain R's `approx(..., rule = 2)` behavior, including
endpoint values outside the reference range. This document records the
implementation; it does not endorse a new extrapolation policy.

Do not replace these files, change their values, or alter interpolation rules
without an explicit scientific review. Any future full-file comparison should
report differences before making replacements.
