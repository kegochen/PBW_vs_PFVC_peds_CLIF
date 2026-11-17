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

The original download date was not supplied. Complete equivalence between
the repository files and the official source files has not yet been verified.
The source citation is not a checksum or proof of file identity.

## Fields used by this pipeline

Script 03 requires:

- `Sex`: sex code used to select the reference subset.
- `Agemos`: reference age in months.
- `M`: median reference measurement.

The supplied file excerpts also include LMS parameters and percentile columns.
In the reviewed excerpts, `M` equals `P50`. That excerpt review does not
validate every row in either file.

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